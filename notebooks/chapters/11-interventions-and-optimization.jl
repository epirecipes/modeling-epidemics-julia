### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 802552b4-06aa-4224-a8e4-c2a72ad45bc9
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 983c632c-dd6e-4257-a7dd-7a7317a4bbaa
begin
    using CairoMakie
    using DiffEqCallbacks
    using OrdinaryDiffEq
    using PlutoUI
end

# ╔═╡ 10e0ae67-3437-46c6-9b53-4208a186e62f
md"""
# Interventions and optimization

Move the intervention start, reduction, and duration sliders. The notebook
simulates a rectangular policy and compares its final cumulative incidence
with the no-policy counterfactual.
"""

# ╔═╡ 6616b485-a374-4e54-81e7-fcbe9f29742d
@bind intervention_start11 PlutoUI.Slider(0:1:25; default = 8, show_value = true)

# ╔═╡ d6f1872d-fc78-4402-bee8-40cc90b2e606
@bind intervention_reduction11 PlutoUI.Slider(0.0:0.05:0.75;
    default = 0.5, show_value = true)

# ╔═╡ 9eaf01e1-0461-4530-a38d-c3ceb7e61b21
@bind intervention_duration11 PlutoUI.Slider(2:1:15; default = 10, show_value = true)

# ╔═╡ 2d9a5577-c9db-4a8d-bc1a-1b32d2e0d8e0
begin
    β11, c11, γ11, N11 = 0.05, 10.0, 0.25, 1000.0
    u011 = [990.0, 10.0, 0.0, 0.0]
    time11 = collect(0.0:0.5:45.0)
    function policy11!(du, u, p, t)
        β, c, γ, reduction, N = p
        incidence = β * c * (1 - reduction) * u[1] * u[2] / N
        du[1] = -incidence
        du[2] = incidence - γ * u[2]
        du[3] = γ * u[2]
        du[4] = incidence
    end

    function scenario11(start_time; duration = 10.0, reduction = 0.5)
        stop_time = start_time + duration
        parameters = [β11, c11, γ11, 0.0, N11]
        function affect!(integrator)
            integrator.p[4] = integrator.t < stop_time ? reduction : 0.0
        end
        callback = PresetTimeCallback([start_time, stop_time], affect!)
        solve(ODEProblem(policy11!, u011, (0.0, last(time11)), parameters),
            Tsit5(); callback = callback, saveat = time11,
            abstol = 1e-8, reltol = 1e-8)
    end

    no_policy11 = solve(ODEProblem(policy11!, u011, (0.0, last(time11)),
        [β11, c11, γ11, 0.0, N11]), Tsit5(); saveat = time11,
        abstol = 1e-8, reltol = 1e-8)
    selected11 = scenario11(intervention_start11;
        duration = intervention_duration11, reduction = intervention_reduction11)
end

# ╔═╡ 0f48a36b-7dd0-4f5a-a6d3-82b5e6edb3dc
begin
    no_array11 = Array(no_policy11)
    selected_array11 = Array(selected11)
    fig11 = Figure(size = (850, 620))
    ax11a = Axis(fig11[1, 1], xlabel = "time", ylabel = "infectious")
    lines!(ax11a, time11, no_array11[2, :], label = "none")
    lines!(ax11a, selected11.t, selected_array11[2, :], label = "selected policy")
    vspan!(ax11a, intervention_start11,
        intervention_start11 + intervention_duration11, color = (:steelblue, 0.15))
    axislegend(ax11a)
    ax11b = Axis(fig11[2, 1], xlabel = "time", ylabel = "cumulative incidence")
    lines!(ax11b, time11, no_array11[4, :], label = "none")
    lines!(ax11b, selected11.t, selected_array11[4, :], label = "selected policy")
    vspan!(ax11b, intervention_start11,
        intervention_start11 + intervention_duration11, color = (:steelblue, 0.15))
    fig11
end

# ╔═╡ e0c456bc-0153-4fbb-bb4d-e6dca00e4e6e
begin
    candidate_starts11 = collect(0.0:1.0:25.0)
    candidate_objectives11 = [
        Array(scenario11(start_time; duration = intervention_duration11,
            reduction = intervention_reduction11))[4, end]
        for start_time in candidate_starts11
    ]
    grid_best11 = candidate_starts11[argmin(candidate_objectives11)]
    conservation_error11 = maximum(
        abs.(vec(sum(selected_array11[1:3, :]; dims = 1)) .- N11),
    )
    cumulative_monotone11 = all(diff(selected_array11[4, :]) .>= -1e-10)
    @assert conservation_error11 < 1e-6
    @assert cumulative_monotone11
    (
        course_unit = "ch11-interventions-and-optimization",
        status = "complete",
        selected_start = intervention_start11,
        selected_reduction = intervention_reduction11,
        selected_duration = intervention_duration11,
        final_incidence_without_policy = no_array11[4, end],
        final_incidence_selected = selected_array11[4, end],
        grid_best_start = grid_best11,
        conservation_error = conservation_error11,
    )
end

# ╔═╡ Cell order:
# ╠═802552b4-06aa-4224-a8e4-c2a72ad45bc9
# ╠═983c632c-dd6e-4257-a7dd-7a7317a4bbaa
# ╟─10e0ae67-3437-46c6-9b53-4208a186e62f
# ╠═6616b485-a374-4e54-81e7-fcbe9f29742d
# ╠═d6f1872d-fc78-4402-bee8-40cc90b2e606
# ╠═9eaf01e1-0461-4530-a38d-c3ceb7e61b21
# ╠═2d9a5577-c9db-4a8d-bc1a-1b32d2e0d8e0
# ╠═0f48a36b-7dd0-4f5a-a6d3-82b5e6edb3dc
# ╠═e0c456bc-0153-4fbb-bb4d-e6dca00e4e6e
