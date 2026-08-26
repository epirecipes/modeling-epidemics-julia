### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 88a15edd-ce29-48dc-aeb0-6fb64429b7ab
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ dc547ed6-be9c-4923-ae58-02fcf7074a73
begin
    using Agents
    using ConcurrentSim
    using ResumableFunctions
    using Random
    using Distributions
    using OrdinaryDiffEq
    using CairoMakie
    using PlutoUI
end

# ╔═╡ 9102da18-896f-4640-8798-93accbf4121e
md"""
# Discrete-event simulation

Companion notebook for `a03-discrete-event-simulation` (see the
[book chapter](../../appendices/03-discrete-event-simulation.qmd) for the
full derivation and both implementations).

Adjust the event-process parameters below and watch a fresh
`EventQueueABM` realisation re-run reactively against the deterministic SIR
ODE. The **fixed-delay recovery** toggle switches the recovery
event's `timing` from the default exponential sampler to the constant
$\tau = 1/\gamma$ used in Exercise 1, connecting this unit back to the fixed
delay of the DDE appendix without any history function.

Both engines below draw partners the same way: a susceptible agent picks one
individual uniformly from all $N$, itself included. A draw that lands on itself,
or on a susceptible or recovered agent, passes without transmission. Counting
those null draws gives a force of infection of $\beta c I / N$, which is the
deterministic curve the stochastic paths are plotted against.
"""

# ╔═╡ 68d3bbd7-d0ab-4cc2-bd82-3d3d9c6e4aa8
md"""
β (transmission probability per contact): $(@bind β PlutoUI.Slider(0.01:0.01:0.15; default = 0.05, show_value = true))

c (contact rate): $(@bind c PlutoUI.Slider(2.0:1.0:20.0; default = 10.0, show_value = true))

γ (recovery rate, mean infectious period 1/γ): $(@bind γ PlutoUI.Slider(0.1:0.05:0.5; default = 0.25, show_value = true))

seed (rng seed for the event queue): $(@bind seed PlutoUI.Slider(1:20; default = 1, show_value = true))

Use a **fixed** recovery delay τ = 1/γ instead of exponential recovery:
$(@bind fixed_delay CheckBox(default = false))
"""

# ╔═╡ 5a85d311-b65c-492a-9604-d5258cbba059
begin
    @agent struct Person(NoSpaceAgent)
        status::Symbol
    end

    function transmit!(agent, model)
        alter = random_agent(model)
        if alter.status == :I && (rand(abmrng(model)) ≤ model.β)
            agent.status = :I
        end
        return nothing
    end
    recover!(agent, model) = (agent.status = :R)

    transmit_propensity(agent, model) = agent.status == :S ? model.c : 0.0
    recovery_propensity(agent, model) = agent.status == :I ? model.γ : 0.0
    fixed_delay_timing(agent, model, propensity) = agent.status == :I ? model.τ : Inf

    susceptible(x) = count(==(:S), x)
    infected(x) = count(==(:I), x)
    recovered(x) = count(==(:R), x)
    nothing
end

# ╔═╡ fee2aca7-cd0e-4b2b-b0a9-5565b0a36efa
function init_eventqueue_model(β, c, γ, N, I0, rng; fixed_delay = false)
    properties = Dict(:β => β, :c => c, :γ => γ, :τ => 1 / γ)
    transmit_event = AgentEvent(action! = transmit!, propensity = transmit_propensity)
    recovery_event = fixed_delay ?
        AgentEvent(action! = recover!, propensity = recovery_propensity, timing = fixed_delay_timing) :
        AgentEvent(action! = recover!, propensity = recovery_propensity)
    events = (transmit_event, recovery_event)
    model = EventQueueABM(Person, events; properties, rng)
    for i in 1:N
        status = i <= I0 ? :I : :S
        add_agent!(Person(; id = i, status = status), model)
    end
    return model
end

# ╔═╡ 74b95ca9-2c37-4ba9-a768-7ca20e4ede4f
begin
    N = 1000
    I0 = 10
    tf = 40.0
    sample_dt = 0.1

    rng = Xoshiro(seed)
    eq_model = init_eventqueue_model(β, c, γ, N, I0, rng; fixed_delay = fixed_delay)
    to_collect = [(:status, f) for f in (susceptible, infected, recovered)]
    eq_data, _ = run!(eq_model, tf; adata = to_collect, when = sample_dt)

    # exponential recovery, so this is the mean field only when fixed_delay is off
    function sir_ode!(du, u, p, t)
        S, I, R = u
        βp, cp, γp = p
        Ntot = S + I + R
        λ = βp * cp * I / Ntot
        du[1] = -λ * S
        du[2] = λ * S - γp * I
        du[3] = γp * I
        return nothing
    end
    sol_ode = solve(ODEProblem(sir_ode!, [Float64(N - I0), Float64(I0), 0.0], (0.0, tf), (β, c, γ)), Tsit5())

    # findmax over solver knots lands between steps, so read the interpolant
    # on a fine grid instead.
    ode_grid = range(0.0, tf; length = 4001)
    ode_I = [sol_ode(t)[2] for t in ode_grid]
    ode_peak, ode_peak_idx = findmax(ode_I)
    ode_peak_t = ode_grid[ode_peak_idx]
    nothing
end

# ╔═╡ 7d7f1fef-98f5-4170-8a97-65800e8d3a17
begin
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "infectious count",
              title = "EventQueueABM ($(fixed_delay ? "fixed" : "exponential") recovery) vs SIR ODE")
    stairs!(ax, eq_data.time, eq_data.infected_status; color = :firebrick, linewidth = 2, label = "I (EventQueueABM)")
    lines!(ax, ode_grid, ode_I; color = :black, linestyle = :dash, linewidth = 2, label = "I (SIR ODE, exponential recovery)")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 44c0a9f7-7a5e-4201-a743-1d46feac505d
begin
    peak_eq, idx_eq = findmax(eq_data.infected_status)
    conserved_eq = all(r -> r.susceptible_status + r.infected_status + r.recovered_status == N, eachrow(eq_data))
    recovery_note = fixed_delay ?
        "The toggle is **on**, so agents recover exactly τ = $(round(1/γ; digits = 2)) after infection. The dashed curve keeps exponential recovery with the same average duration, so it is a comparison at equal means rather than the mean field of this process. Matching a fixed infectious period needs an infection-age, renewal or delay formulation, which is what the DDE appendix builds." :
        "The toggle is **off**, so recovery is exponential in both, and the dashed curve is the mean field the stochastic path fluctuates around."
    recovery_label = fixed_delay ? "fixed τ = $(round(1/γ; digits = 2))" : "exponential"

    Markdown.parse("""
    **Current run** (β = $(β), c = $(c), γ = $(γ), seed = $(seed), recovery = $(recovery_label)):

    - Highest *recorded* EventQueueABM count: **$(peak_eq)** at t ≈ $(round(eq_data.time[idx_eq]; digits = 2)). Agents.jl records once the model has advanced by at least $(sample_dt) time units, at an event time, rather than after every event, so the true path maximum can fall between recorded observations and be higher.
    - SIR ODE peak: **$(round(ode_peak; digits = 1))** at t ≈ $(round(ode_peak_t; digits = 2)), read off the solver's interpolant on a fine grid rather than its adaptive knots.
    - Conservation of N holds in every recorded row: **$(conserved_eq)**

    $(recovery_note)

    With the fixed-delay toggle on, the recovery time is carried by the event the queue
    has already scheduled for that agent, not by a field on `Person`, which stores only
    `status`. That is how an individual-level model reaches a fixed infectious period
    with no history function.
    """)
end

# ╔═╡ b5a7713e-1ad6-4f56-80fa-94fedcd0e333
md"""
## Alternative engine: `ConcurrentSim.jl` with a limited resource

`ConcurrentSim` ships a `Resource`: a counted pool with a request queue. Nothing
in `EventQueueABM` plays that role, so a queue there is yours to build. The beds
below are **triaged**: each infectious individual is either severe or mild, and
severe cases are served first while beds are scarce.

Exercise 2 builds the plain first-come, first-served version. Try it before
reading on, since the code below already queues.

Three assumptions this model makes:

- **Every** infectious individual needs a bed, mild cases included, and the
  recovery clock starts only once a bed is held. Time spent queueing therefore
  lengthens the infectious period and, when a queue forms, can amplify
  transmission relative to an otherwise identical unlimited-bed process.
- Triage is **non-preemptive**. A severe case jumps ahead of mild cases *waiting*
  in the queue, but never displaces a mild case already occupying a bed.
- The bed is requested with `lock(m.beds; priority = ...)`, where lower numbers
  are served first, so severe cases carry 0 and mild cases 1. The chapter writes
  the same call as `request`, which is ConcurrentSim's other name for `lock`.
"""

# ╔═╡ 1cabfe77-0395-4843-b785-bfc19fd24bf4
md"""
Number of hospital beds: $(@bind num_beds PlutoUI.Slider(10:10:1000; default = 100, show_value = true))

Probability that a case is severe: $(@bind severe_fraction PlutoUI.Slider(0.0:0.05:0.5; default = 0.2, show_value = true))
"""

# ╔═╡ 556260cf-77ba-4529-8e54-c2377478bc4a
module DES
    using ConcurrentSim, ResumableFunctions, Distributions, Random

    mutable struct SIRPerson
        id::Int
        status::Symbol
        severe::Bool
    end

    mutable struct SIRModel
        sim::Simulation
        β::Float64
        c::Float64
        γ::Float64
        beds::Resource
        ta::Vector{Float64}
        Ia::Vector{Int}
        waits_severe::Vector{Float64}
        waits_mild::Vector{Float64}
        requests_severe::Int
        requests_mild::Int
        people::Vector{SIRPerson}
    end

    function record_i!(m::SIRModel)
        push!(m.ta, now(m.sim))
        push!(m.Ia, count(p -> p.status == :I, m.people))
        return nothing
    end

    # "n/a" rather than 0.0: a group with no admissions has no mean wait
    mean_wait(v) = isempty(v) ? "n/a" : round(sum(v) / length(v); digits = 2)

    # tested on the raw waits, since mean_wait and waited_share_admitted round for display
    any_waited(m::SIRModel) =
        any(>(1e-9), m.waits_severe) || any(>(1e-9), m.waits_mild)

    # denominator is completed admissions; those still queueing are in still_queued
    waited_share_admitted(m::SIRModel) = begin
        all_waits = vcat(m.waits_severe, m.waits_mild)
        isempty(all_waits) ? 0.0 :
            round(count(>(1e-9), all_waits) / length(all_waits); digits = 2)
    end

    still_queued(m::SIRModel) =
        (m.requests_severe - length(m.waits_severe), m.requests_mild - length(m.waits_mild))

    @resumable function live(sim::Simulation, person::SIRPerson, m::SIRModel, rng)
        while person.status == :S
            @yield timeout(sim, rand(rng, Exponential(1 / m.c)))
            alter = rand(rng, m.people)      # uniform over all N, person included
            if alter.status == :I && rand(rng) < m.β
                person.status = :I
                record_i!(m)
            end
        end
        if person.status == :I
            # counted at request time, before the queue is entered
            if person.severe
                m.requests_severe += 1
            else
                m.requests_mild += 1
            end
            request_time = now(sim)
            @yield lock(m.beds; priority = person.severe ? 0 : 1)
            waited = now(sim) - request_time
            push!(person.severe ? m.waits_severe : m.waits_mild, waited)
            @yield timeout(sim, rand(rng, Exponential(1 / m.γ)))
            person.status = :R
            @yield unlock(m.beds)
            record_i!(m)
        end
    end

    function init_concurrentsim_model(S, I, R, β, c, γ, nbeds, severe_fraction, rng)
        sim = Simulation()
        beds = Resource(sim, nbeds)
        people = SIRPerson[]
        for i in 1:(S + I + R)
            status = i <= S ? :S : (i <= S + I ? :I : :R)
            push!(people, SIRPerson(i, status, rand(rng) < severe_fraction))
        end
        SIRModel(sim, β, c, γ, beds, [0.0], [I], Float64[], Float64[], 0, 0, people)
    end
end

# ╔═╡ d3117a11-7b9f-4d57-b3d0-b24738b05b9e
begin
    cs_rng = Xoshiro(seed)
    cs_model = DES.init_concurrentsim_model(N - I0, I0, 0, β, c, γ, num_beds,
                                            severe_fraction, cs_rng)
    for person in cs_model.people
        @process DES.live(cs_model.sim, person, cs_model, cs_rng)
    end
    run(cs_model.sim, tf)
    nothing
end

# ╔═╡ e3c3507c-52a5-4475-800e-0a629195348d
begin
    fig2 = Figure(size = (720, 380))
    ax2 = Axis(fig2[1, 1]; xlabel = "time", ylabel = "infectious count",
               title = "ConcurrentSim: $(num_beds) triaged beds")
    stairs!(ax2, cs_model.ta, cs_model.Ia; color = :steelblue, linewidth = 2, label = "I (ConcurrentSim)")
    lines!(ax2, ode_grid, ode_I; color = :black, linestyle = :dash, linewidth = 2, label = "I (SIR ODE, unlimited beds)")
    axislegend(ax2; position = :rt)
    fig2
end

# ╔═╡ 48b8d3ee-63cc-4439-919e-1180aaac01fd
begin
    cs_peak = maximum(cs_model.Ia)
    q_sev, q_mild = DES.still_queued(cs_model)
    adm_sev = length(cs_model.waits_severe)
    adm_mild = length(cs_model.waits_mild)

    # branch on this run's own numbers
    verdict = if !DES.any_waited(cs_model) && q_sev + q_mild == 0
        "At $(num_beds) beds nobody queued in this run, so the bed limit never bound and the blue curve is the unconstrained process."
    elseif cs_peak > 1.05 * ode_peak
        "Queueing occurred here, and this run's peak of $(cs_peak) exceeds the unlimited-bed ODE peak of $(round(ode_peak; digits = 1)). That is consistent with queue-induced amplification, though part of any single run's gap is ordinary stochastic variation."
    else
        "Here some individuals queued, but the peak of $(cs_peak) did not rise clearly above the unlimited-bed ODE peak of $(round(ode_peak; digits = 1))."
    end

    # interpolate into a plain string, then parse: md"..." mis-pairs the **
    # markers when a line carries several values
    Markdown.parse("""
    **Bed-constrained run** ($(num_beds) beds, severity probability $(round(Int, 100 * severe_fraction))%, unconstrained SIR R₀ = $(round(β * c / γ; digits = 2))):

    - Peak infectious count: **$(cs_peak)**
    - Severe cases: **$(adm_sev)** admitted of $(cs_model.requests_severe), mean wait **$(DES.mean_wait(cs_model.waits_severe))**
    - Mild cases: **$(adm_mild)** admitted of $(cs_model.requests_mild), mean wait **$(DES.mean_wait(cs_model.waits_mild))**
    - Still queueing when the run stopped: severe **$(q_sev)**, mild **$(q_mild)**
    - Share of *completed admissions* that waited at all: **$(DES.waited_share_admitted(cs_model))**

    The wait figures cover only admissions completed by t = $(tf). Anyone still queueing
    at that point appears on the "still queueing" line and is absent from the means, so
    under heavy congestion these means understate the true waiting, most of all for mild
    cases at the back of the queue.

    $(verdict)

    When both groups are represented and beds are scarce, compare their mean waits.
    Change the seed a few times: severe cases wait less each time, and that is the
    triage priority. The size of the gap in any one run also carries randomness. How
    many beds it takes before the queue stops binding depends on β, c and γ as well, so
    it shifts as you change them.
    """)
end

# ╔═╡ ec7779dc-56a8-4b33-83cb-f3269e2b0c50
(
    course_unit = "a03-discrete-event-simulation",
    status = "complete",
    β = β,
    c = c,
    γ = γ,
    seed = seed,
    fixed_delay = fixed_delay,
    num_beds = num_beds,
    severe_fraction = severe_fraction,
    eq_peak = peak_eq,
    cs_peak = cs_peak,
)

# ╔═╡ Cell order:
# ╠═88a15edd-ce29-48dc-aeb0-6fb64429b7ab
# ╠═dc547ed6-be9c-4923-ae58-02fcf7074a73
# ╟─9102da18-896f-4640-8798-93accbf4121e
# ╟─68d3bbd7-d0ab-4cc2-bd82-3d3d9c6e4aa8
# ╠═5a85d311-b65c-492a-9604-d5258cbba059
# ╠═fee2aca7-cd0e-4b2b-b0a9-5565b0a36efa
# ╠═74b95ca9-2c37-4ba9-a768-7ca20e4ede4f
# ╟─7d7f1fef-98f5-4170-8a97-65800e8d3a17
# ╟─44c0a9f7-7a5e-4201-a743-1d46feac505d
# ╟─b5a7713e-1ad6-4f56-80fa-94fedcd0e333
# ╟─1cabfe77-0395-4843-b785-bfc19fd24bf4
# ╠═556260cf-77ba-4529-8e54-c2377478bc4a
# ╠═d3117a11-7b9f-4d57-b3d0-b24738b05b9e
# ╟─e3c3507c-52a5-4475-800e-0a629195348d
# ╟─48b8d3ee-63cc-4439-919e-1180aaac01fd
# ╠═ec7779dc-56a8-4b33-83cb-f3269e2b0c50
