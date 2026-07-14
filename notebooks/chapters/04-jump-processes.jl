### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 7e6f32ef-1e69-4c5a-a9b4-9aa8d4e91f5a
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 0b6bb2cf-f4db-41dd-9afc-20cba6de6f4d
begin
    using CairoMakie
    using PlutoUI
    using Random
    using StableRNGs
end

# ╔═╡ 8e28d52c-d88d-4de8-9e8c-2d4e8d60999b
md"""
# Continuous-time Markov chains and jump processes

This notebook makes the Gillespie direct algorithm explicit.  Move the
controls to change the epidemic while the seeded simulation remains
reproducible.
"""

# ╔═╡ 2505a9fd-c001-4b95-b59a-bfe3bd6ecbf2
begin
    @bind transmission PlutoUI.Slider(0.02:0.01:0.10; default = 0.05, show_value = true)
    @bind initial_infected PlutoUI.Slider(1:20; default = 10, show_value = true)
    @bind horizon PlutoUI.Slider(10:5:60; default = 40, show_value = true)
end

# ╔═╡ 6d17e08c-7db3-43ae-b8bc-d3a199c912ba
md"""
The infection propensity is \(a_1=\beta cSI/N\), and the recovery
propensity is \(a_2=\gamma I\).  Each event changes an integer state by one.
"""

# ╔═╡ 69e9f825-a8bc-4cf1-bd22-d99edc15bb1c
function simulate_ctmc(N, I₀, β, c, γ, tmax, rng)
    t = 0.0
    state = [N - I₀, I₀, 0]
    event_times = [t]
    states = [copy(state)]
    while t < tmax && state[2] > 0
        infection = β * c * state[1] * state[2] / N
        recovery = γ * state[2]
        total = infection + recovery
        total == 0 && break
        t += randexp(rng) / total
        t > tmax && break
        if rand(rng) * total < infection
            state[1] -= 1
            state[2] += 1
        else
            state[2] -= 1
            state[3] += 1
        end
        push!(event_times, t)
        push!(states, copy(state))
    end
    push!(event_times, tmax)
    push!(states, copy(state))
    return (; event_times, states)
end

# ╔═╡ d4207f34-b1dd-4f70-81f8-873ddc33e0fd
path = simulate_ctmc(
    1000,
    initial_infected,
    transmission,
    10.0,
    0.25,
    Float64(horizon),
    StableRNG(20250711),
)

# ╔═╡ 660f51a2-0e99-4eb2-90dd-3a21d18ddb0f
begin
    event_states = reduce(hcat, path.states)
    mass_error = maximum(abs.(sum(event_states; dims = 1) .- 1000))
    peak_infected = maximum(event_states[2, :])
    @assert mass_error == 0
    summary = (; events = length(path.event_times) - 2, peak_infected, mass_error)
end

# ╔═╡ 2c1e2f42-a213-4b0a-8ca2-07a36d4bc681
begin
    fig = Figure(size = (760, 420))
    ax = Axis(fig[1, 1], xlabel = "time", ylabel = "people",
        title = "Seeded Gillespie SIR path")
    names = ["S", "I", "R"]
    colors = [:steelblue, :firebrick, :seagreen]
    for j in 1:3
        stairs!(ax, path.event_times, event_states[j, :]; color = colors[j],
            linewidth = 2, label = names[j])
    end
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 9f0c2d82-0b9a-4a4a-8852-4ab9b0a6d6ea
md"""
**Check:** `summary` reports exact conservation and the number of events.
The path is piecewise constant; `saveat`-style plotting would only sample
this event history, not create the events.
"""

# ╔═╡ 4a4f7f4f-8d8d-4f50-bec6-3a40cce5bb43
(
    course_unit = "ch04-jump-processes",
    status = "complete",
    controls = (; transmission, initial_infected, horizon),
)

# ╔═╡ Cell order:
# ╠═7e6f32ef-1e69-4c5a-a9b4-9aa8d4e91f5a
# ╠═0b6bb2cf-f4db-41dd-9afc-20cba6de6f4d
# ╟─8e28d52c-d88d-4de8-9e8c-2d4e8d60999b
# ╠═2505a9fd-c001-4b95-b59a-bfe3bd6ecbf2
# ╟─6d17e08c-7db3-43ae-b8bc-d3a199c912ba
# ╠═69e9f825-a8bc-4cf1-bd22-d99edc15bb1c
# ╠═d4207f34-b1dd-4f70-81f8-873ddc33e0fd
# ╠═660f51a2-0e99-4eb2-90dd-3a21d18ddb0f
# ╠═2c1e2f42-a213-4b0a-8ca2-07a36d4bc681
# ╟─9f0c2d82-0b9a-4a4a-8852-4ab9b0a6d6ea
# ╠═4a4f7f4f-8d8d-4f50-bec6-3a40cce5bb43
