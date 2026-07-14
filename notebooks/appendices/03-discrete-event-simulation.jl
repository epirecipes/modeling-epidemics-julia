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
`EventQueueABM` realization re-run reactively against the deterministic ODE
mean field. The **fixed-delay recovery** toggle switches the recovery
event's `timing` from the default exponential sampler to the constant
$\tau = 1/\gamma$ used in Exercise 1, connecting this unit back to the fixed
delay of the DDE appendix without any history function.
"""

# ╔═╡ 68d3bbd7-d0ab-4cc2-bd82-3d3d9c6e4aa8
md"""
β (transmission probability per contact): $(@bind β Slider(0.01:0.01:0.15; default = 0.05, show_value = true))

c (contact rate): $(@bind c Slider(2.0:1.0:20.0; default = 10.0, show_value = true))

γ (recovery rate, mean infectious period 1/γ): $(@bind γ Slider(0.1:0.05:0.5; default = 0.25, show_value = true))

seed (rng seed for the event queue): $(@bind seed Slider(1:20; default = 1, show_value = true))

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

    rng = Xoshiro(seed)
    eq_model = init_eventqueue_model(β, c, γ, N, I0, rng; fixed_delay = fixed_delay)
    to_collect = [(:status, f) for f in (susceptible, infected, recovered)]
    eq_data, _ = run!(eq_model, tf; adata = to_collect)

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
    nothing
end

# ╔═╡ 7d7f1fef-98f5-4170-8a97-65800e8d3a17
begin
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "infectious count",
              title = "EventQueueABM ($(fixed_delay ? "fixed" : "exponential") recovery) vs ODE mean field")
    stairs!(ax, eq_data.time, eq_data.infected_status; color = :firebrick, linewidth = 2, label = "I (EventQueueABM)")
    lines!(ax, sol_ode.t, sol_ode[2, :]; color = :black, linestyle = :dash, linewidth = 2, label = "I (ODE mean field)")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 44c0a9f7-7a5e-4201-a743-1d46feac505d
begin
    peak_eq, idx_eq = findmax(eq_data.infected_status)
    peak_ode2, idx_ode2 = findmax(sol_ode[2, :])
    conserved_eq = all(r -> r.susceptible_status + r.infected_status + r.recovered_status == N, eachrow(eq_data))
    md"""
    **Current run** (β = $(β), c = $(c), γ = $(γ), seed = $(seed), recovery = $(fixed_delay ? "fixed τ = $(round(1/γ; digits=2))" : "exponential")):

    - EventQueueABM peak infectious count: **$(peak_eq)** at t ≈ $(round(eq_data.time[idx_eq]; digits = 2))
    - ODE mean-field peak: **$(round(peak_ode2; digits = 1))** at t ≈ $(round(sol_ode.t[idx_ode2]; digits = 2))
    - Conservation of N holds in every recorded row: **$(conserved_eq)**

    With the fixed-delay toggle on, every infectious agent carries its own
    infection time and recovers exactly τ time units later — no history
    function needed, unlike the DDE appendix, because the individual
    retains the information a population-level state discards.
    """
end

# ╔═╡ b5a7713e-1ad6-4f56-80fa-94fedcd0e333
md"""
## Alternative engine: `ConcurrentSim.jl` with a limited resource

Set the number of hospital beds below to see how resource contention
(Exercise 2) changes the infectious curve: an infectious individual must
successfully `request` a bed before its recovery `timeout` starts, and
`release`s it on recovery. With few beds, individuals stay counted as
infectious for longer (waiting **and** occupying a bed), which — unlike in
`EventQueueABM`, which has no resource abstraction — measurably prolongs and
amplifies the epidemic.
"""

# ╔═╡ 1cabfe77-0395-4843-b785-bfc19fd24bf4
md"""
Number of hospital beds: $(@bind num_beds Slider(10:10:1000; default = 1000, show_value = true))
"""

# ╔═╡ 556260cf-77ba-4529-8e54-c2377478bc4a
begin
    mutable struct SIRPerson
        id::Int
        status::Symbol
    end

    mutable struct SIRModel
        sim::Simulation
        β::Float64
        c::Float64
        γ::Float64
        beds::Resource
        ta::Vector{Float64}
        Ia::Vector{Int}
        wait_times::Vector{Float64}
        people::Vector{SIRPerson}
    end

    function record_i!(m::SIRModel)
        push!(m.ta, now(m.sim))
        push!(m.Ia, count(p -> p.status == :I, m.people))
        return nothing
    end

    @resumable function live(sim::Simulation, person::SIRPerson, m::SIRModel, rng)
        while person.status == :S
            @yield timeout(sim, rand(rng, Exponential(1 / m.c)))
            alter = person
            while alter === person
                alter = rand(rng, m.people)
            end
            if alter.status == :I && rand(rng) < m.β
                person.status = :I
                record_i!(m)
            end
        end
        if person.status == :I
            request_time = now(sim)
            @yield request(m.beds)
            push!(m.wait_times, now(sim) - request_time)
            @yield timeout(sim, rand(rng, Exponential(1 / m.γ)))
            person.status = :R
            @yield release(m.beds)
            record_i!(m)
        end
    end

    function init_concurrentsim_model(S, I, R, β, c, γ, nbeds)
        sim = Simulation()
        beds = Resource(sim, nbeds)
        people = SIRPerson[]
        for i in 1:S
            push!(people, SIRPerson(i, :S))
        end
        for i in 1:I
            push!(people, SIRPerson(S + i, :I))
        end
        for i in 1:R
            push!(people, SIRPerson(S + I + i, :R))
        end
        SIRModel(sim, β, c, γ, beds, [0.0], [I], Float64[], people)
    end
    nothing
end

# ╔═╡ d3117a11-7b9f-4d57-b3d0-b24738b05b9e
begin
    cs_rng = Xoshiro(seed)
    cs_model = init_concurrentsim_model(N - I0, I0, 0, β, c, γ, num_beds)
    for person in cs_model.people
        @process live(cs_model.sim, person, cs_model, cs_rng)
    end
    run(cs_model.sim, tf)
    nothing
end

# ╔═╡ e3c3507c-52a5-4475-800e-0a629195348d
begin
    fig2 = Figure(size = (720, 380))
    ax2 = Axis(fig2[1, 1]; xlabel = "time", ylabel = "infectious count",
               title = "ConcurrentSim with $(num_beds) beds")
    stairs!(ax2, cs_model.ta, cs_model.Ia; color = :steelblue, linewidth = 2, label = "I (ConcurrentSim)")
    lines!(ax2, sol_ode.t, sol_ode[2, :]; color = :black, linestyle = :dash, linewidth = 2, label = "I (ODE mean field)")
    axislegend(ax2; position = :rt)
    fig2
end

# ╔═╡ 48b8d3ee-63cc-4439-919e-1180aaac01fd
md"""
**Bed-constrained run** ($(num_beds) beds):

- Peak infectious count: **$(maximum(cs_model.Ia))**
- Mean wait time for a bed: **$(isempty(cs_model.wait_times) ? 0.0 : round(sum(cs_model.wait_times) / length(cs_model.wait_times); digits = 2))**
- Fraction of infectious individuals who waited at all: **$(isempty(cs_model.wait_times) ? 0.0 : round(count(>(1e-9), cs_model.wait_times) / length(cs_model.wait_times); digits = 2))**

Reducing the bed count below 1000 makes waiting — and hence the effective
infectious period — longer, which is why the peak can grow well beyond the
unconstrained model even though β, c, and γ are unchanged. `EventQueueABM`
has no built-in analogue for this kind of shared, capacity-limited resource.
"""

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
    eq_peak = peak_eq,
    cs_peak = maximum(cs_model.Ia),
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
