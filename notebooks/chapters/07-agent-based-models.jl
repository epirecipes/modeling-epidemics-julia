### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ e3c4ef43-1852-4e2a-a531-04e201cb2745
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ e18a72cb-31b8-48a5-8cba-7f5ea4b2280e
begin
    using Agents
    using CairoMakie
    using Distributions
    using PlutoUI
    using StableRNGs
end

# ╔═╡ 70e6c414-e8c4-4643-bb51-527bff4b047f
md"""
# Agent-based models

Each `Person` is a stateful individual.  The two-phase update first reads all
old statuses and stores `next_status`; the model step then commits them.  This
prevents activation order from changing who can transmit.

This notebook starts from 500 people with 5 infected, half the 1000 and 10 of
the chapter.  Halving both keeps the same starting proportion infected, so the
epidemic reaches about the same final size, and the smaller population keeps
each slider change quick enough to explore.  Fewer agents also means more
noise: the height and timing of the peak move further from seed to seed than
they do in the chapter.
"""

# ╔═╡ 238accad-5aa0-4cb6-8b2b-a403da8430d0
@bind population PlutoUI.Slider(100:50:500; default = 500, show_value = true)

# ╔═╡ 8e396acc-ad78-4c44-afb6-82c7b6db6c66
@bind initial_infected PlutoUI.Slider(1:10; default = 5, show_value = true)

# ╔═╡ 966785e5-2ac5-48fb-a81e-7f5832e26c26
@bind contact_rate PlutoUI.Slider(4.0:1.0:14.0; default = 10.0, show_value = true)

# ╔═╡ bb5f836a-1549-4976-91f1-4ccf5b3a40fe
@bind horizon PlutoUI.Slider(10:5:40; default = 30, show_value = true)

# ╔═╡ 59a8f381-2b02-44ec-84b9-d6e1a32f3d65
begin
    β = 0.05
    γ = 0.25
    Δt = 0.1

    @agent struct Person(NoSpaceAgent)
        status::Symbol
        next_status::Symbol
    end

    function agent_step!(person, model)
        person.next_status = person.status
        if person.status == :S
            contacts = rand(abmrng(model), Poisson(model.c * model.Δt))
            for _ in 1:contacts
                contact = random_agent(model)
                if contact.status == :I && rand(abmrng(model)) < model.β
                    person.next_status = :I
                    break
                end
            end
        elseif person.status == :I &&
               rand(abmrng(model)) < 1 - exp(-model.γ * model.Δt)
            person.next_status = :R
        end
        return nothing
    end

    function model_step!(model)
        for person in allagents(model)
            person.status = person.next_status
        end
        return nothing
    end

    susceptible_count(model) = count(p.status == :S for p in allagents(model))
    infected_count(model) = count(p.status == :I for p in allagents(model))
    recovered_count(model) = count(p.status == :R for p in allagents(model))
end

# ╔═╡ 7765cfb5-cb4e-450b-a19f-dd8f88067365
function initialize_model(seed)
    model = StandardABM(
        Person;
        agent_step!,
        model_step!,
        properties = (β = β, c = contact_rate, γ = γ, Δt = Δt),
        rng = StableRNG(seed),
    )
    for id in 1:population
        status = id <= initial_infected ? :I : :S
        add_agent!(Person(; id, status, next_status = status), model)
    end
    return model
end

# ╔═╡ 956e7e93-a51d-41c1-8eef-008f4d27542a
begin
    nsteps = round(Int, horizon / Δt)
    model = initialize_model(20250711)
    _, model_data = run!(
        model,
        nsteps;
        mdata = [susceptible_count, infected_count, recovered_count],
        showprogress = false,
    )
    abm_states = Matrix(model_data[:, [
        :susceptible_count, :infected_count, :recovered_count,
    ]])'
    abm_times = Float64.(model_data.time) .* Δt
end

# ╔═╡ ed9bba52-531a-49bf-b14b-7e38dd15c100
begin
    total_error = maximum(abs.(vec(sum(abm_states; dims = 1)) .- population))
    peak_infected = maximum(abm_states[2, :])
    @assert total_error == 0
    (; total_error, peak_infected, nsteps)
end

# ╔═╡ 3acb0f23-74a0-42a2-b22b-5bb409e1d22d
begin
    fig = Figure(size = (760, 420))
    ax = Axis(fig[1, 1], xlabel = "time", ylabel = "people",
        title = "Well-mixed individual-based SIR")
    for (j, (name, color)) in enumerate(zip(["S", "I", "R"],
                                             [:steelblue, :firebrick, :seagreen]))
        stairs!(ax, abm_times, abm_states[j, :]; color, linewidth = 2,
            label = name)
    end
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ a358036d-b6e9-4a8f-a6c1-0e12b86d2fd6
md"""
The controls change an individual-level realization, not just a deterministic
curve.  Population conservation remains exact because every agent has one and
only one status at every recorded step.
"""

# ╔═╡ 9a1527c1-f05b-46cb-959a-b6502f0f0aef
(
    course_unit = "ch07-agent-based-models",
    status = "complete",
    controls = (; population, initial_infected, contact_rate, horizon),
)

# ╔═╡ Cell order:
# ╠═e3c4ef43-1852-4e2a-a531-04e201cb2745
# ╠═e18a72cb-31b8-48a5-8cba-7f5ea4b2280e
# ╟─70e6c414-e8c4-4643-bb51-527bff4b047f
# ╠═238accad-5aa0-4cb6-8b2b-a403da8430d0
# ╠═8e396acc-ad78-4c44-afb6-82c7b6db6c66
# ╠═966785e5-2ac5-48fb-a81e-7f5832e26c26
# ╠═bb5f836a-1549-4976-91f1-4ccf5b3a40fe
# ╠═59a8f381-2b02-44ec-84b9-d6e1a32f3d65
# ╠═7765cfb5-cb4e-450b-a19f-dd8f88067365
# ╠═956e7e93-a51d-41c1-8eef-008f4d27542a
# ╠═ed9bba52-531a-49bf-b14b-7e38dd15c100
# ╠═3acb0f23-74a0-42a2-b22b-5bb409e1d22d
# ╟─a358036d-b6e9-4a8f-a6c1-0e12b86d2fd6
# ╠═9a1527c1-f05b-46cb-959a-b6502f0f0aef
