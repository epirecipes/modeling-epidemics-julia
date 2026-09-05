### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 9ae06a3f-c63a-4fd3-8648-aab037c12fc1
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ a5b80b48-2fc3-4bd2-9130-f80f4c1f6af9
begin
    using CairoMakie
    using CSV
    using DataFrames
    using DiffEqCallbacks
    using OrdinaryDiffEq
    using PlutoUI
    using SHA
    using TOML
end

# ╔═╡ 07f973e7-1dce-43a5-be0e-ff7ae40028dc
md"""
# Data-to-decision capstone

This notebook uses only the committed Eyam CSV and its provenance metadata.
If the CSV is unavailable, a deterministic SIR trajectory with the same schema
is used; no network request is made. Sliders explore a conditional policy.
"""

# ╔═╡ de216038-3dd7-4865-bebe-146184a12995
@bind policy_start12 PlutoUI.Slider(0.0:0.5:4.0; default = 1.0, show_value = true)

# ╔═╡ a97fd0d7-245f-40db-bf9a-5f6d04b7c5fa
@bind policy_reduction12 PlutoUI.Slider(0.0:0.1:0.8; default = 0.5, show_value = true)

# ╔═╡ 5cb729e9-4f4f-44cf-b88a-c87f48ee72db
begin
    root12 = normpath(joinpath(@__DIR__, "..", ".."))
    data_path12 = joinpath(root12, "data", "eyam_plague_1666.csv")
    metadata12 = TOML.parsefile(joinpath(root12, "data", "eyam_plague_1666.toml"))
    population12 = Float64(metadata12["population"])

    function fallback12(N)
        times = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
        u0 = [N - 7.0, 7.0, 0.0]
        rhs!(du, u, p, t) = begin
            λ, γ = p
            incidence = λ * u[1] * u[2] / N
            du[1] = -incidence
            du[2] = incidence - γ * u[2]
            du[3] = γ * u[2]
        end
        states = Array(solve(ODEProblem(rhs!, u0, (0.0, 4.0), [1.0, 0.5]),
            Tsit5(); saveat = times, abstol = 1e-9, reltol = 1e-9))
        DataFrame(time_months = times, S = states[1, :],
            I = states[2, :], R = states[3, :])
    end

    committed12 = isfile(data_path12)
    data12 = committed12 ? CSV.read(data_path12, DataFrame) : fallback12(population12)
    data_source12 = committed12 ? "committed Eyam CSV" : "deterministic fallback"
    if committed12
        @assert bytes2hex(sha256(read(data_path12))) == metadata12["csv_sha256"]
    end
    @assert names(data12) == ["time_months", "S", "I", "R"]
    @assert nrow(data12) == 8
    @assert all(isapprox.(data12.S .+ data12.I .+ data12.R, population12;
        atol = 1e-8))
end

# ╔═╡ 5dbbd9c0-8145-4a2d-9cc3-e57dfcca08cf
begin
    times12 = Float64.(data12.time_months)
    initial12 = Float64.(collect(data12[1, [:S, :I, :R]]))
    observed12 = permutedims(Matrix{Float64}(data12[:, [:S, :I, :R]]))
    function fit_rhs12!(du, u, p, t)
        λ, γ = p
        incidence = λ * u[1] * u[2] / population12
        du[1] = -incidence
        du[2] = incidence - γ * u[2]
        du[3] = γ * u[2]
    end
    function predict12(p)
        Array(solve(ODEProblem(fit_rhs12!, initial12,
            (first(times12), last(times12)), p), Tsit5();
            saveat = times12, abstol = 1e-8, reltol = 1e-8))
    end
    candidates12 = [
        (sum(((predict12([λ, γ]) .- observed12) ./ population12).^2), λ, γ)
        for λ in range(0.5, 6.0; length = 12)
        for γ in range(0.5, 4.0; length = 12)
    ]
    best12 = candidates12[argmin(first.(candidates12))]
    fitted12 = [best12[2], best12[3]]
    prediction12 = predict12(fitted12)
    residuals12 = prediction12 .- observed12
end

# ╔═╡ 1f90a1aa-076e-46ca-9c54-3c5f1cf3cd82
begin
    fig_fit12 = Figure(size = (850, 640))
    ax_fit12 = Axis(fig_fit12[1, 1], xlabel = "months since 18 June 1666",
        ylabel = "people", title = "Fit: $(data_source12)")
    for (j, label) in enumerate(["S", "I", "R"])
        lines!(ax_fit12, times12, prediction12[j, :], label = "$label model")
        scatter!(ax_fit12, times12, observed12[j, :], label = "$label data")
    end
    axislegend(ax_fit12, nbanks = 2)
    ax_res12 = Axis(fig_fit12[2, 1], xlabel = "time", ylabel = "model − data")
    for j in 1:3
        lines!(ax_res12, times12, residuals12[j, :], label = ["S", "I", "R"][j])
    end
    hlines!(ax_res12, [0.0], linestyle = :dash)
    fig_fit12
end

# ╔═╡ e7b6c2d8-8cd9-4b1e-93be-7d8a6c5a9b7b
begin
    function decision_rhs12!(du, u, p, t)
        λ, γ, reduction, N = p
        incidence = (1 - reduction) * λ * u[1] * u[2] / N
        du[1] = -incidence
        du[2] = incidence - γ * u[2]
        du[3] = γ * u[2]
        du[4] = incidence
    end
    function decision12(start; reduction = 0.5, duration = 1.0)
        stop = start + duration
        p = [fitted12[1], fitted12[2], 0.0, population12]
        problem = ODEProblem(decision_rhs12!, [initial12..., 0.0],
            (0.0, 8.0), p)
        if reduction == 0.0
            return solve(problem, Tsit5(); saveat = 0.1)
        end
        affect!(integrator) = (integrator.p[3] =
            integrator.t < stop ? reduction : 0.0)
        callback = PresetTimeCallback([start, stop], affect!)
        solve(problem, Tsit5(); callback = callback, saveat = 0.1)
    end
    no_policy12 = decision12(0.0; reduction = 0.0)
    selected12 = decision12(policy_start12; reduction = policy_reduction12)
    no_array12, selected_array12 = Array(no_policy12), Array(selected12)
end

# ╔═╡ 3a76057d-5acd-4f4e-a319-f4bfb6bcf9c2
begin
    fig_decision12 = Figure(size = (850, 430))
    ax_decision12 = Axis(fig_decision12[1, 1], xlabel = "months",
        ylabel = "cumulative incidence")
    lines!(ax_decision12, no_policy12.t, no_array12[4, :], label = "none")
    lines!(ax_decision12, selected12.t, selected_array12[4, :],
        label = "selected policy")
    vspan!(ax_decision12, policy_start12, policy_start12 + 1.0,
        color = (:steelblue, 0.15))
    axislegend(ax_decision12)
    fig_decision12
end

# ╔═╡ 829dcb3f-112a-496b-ae99-94322920c2e4
begin
    conservation_error12 = maximum(
        abs.(vec(sum(prediction12; dims = 1)) .- population12),
    )
    policy_reduction12_value = no_array12[4, end] - selected_array12[4, end]
    @assert conservation_error12 < 1e-6
    @assert all(isfinite, prediction12)
    @assert policy_reduction12_value >= -1e-8
    (
        course_unit = "ch15-data-to-decision-capstone",
        status = "complete",
        data_source = data_source12,
        fitted_effective_rate = fitted12[1],
        fitted_removal_rate = fitted12[2],
        policy_start = policy_start12,
        policy_reduction = policy_reduction12,
        final_incidence_reduction = policy_reduction12_value,
        conservation_error = conservation_error12,
        interpretation = "counterfactual under a limited homogeneous SIR model",
    )
end

# ╔═╡ Cell order:
# ╠═9ae06a3f-c63a-4fd3-8648-aab037c12fc1
# ╠═a5b80b48-2fc3-4bd2-9130-f80f4c1f6af9
# ╟─07f973e7-1dce-43a5-be0e-ff7ae40028dc
# ╠═de216038-3dd7-4865-bebe-146184a12995
# ╠═a97fd0d7-245f-40db-bf9a-5f6d04b7c5fa
# ╠═5cb729e9-4f4f-44cf-b88a-c87f48ee72db
# ╠═5dbbd9c0-8145-4a2d-9cc3-e57dfcca08cf
# ╠═1f90a1aa-076e-46ca-9c54-3c5f1cf3cd82
# ╠═e7b6c2d8-8cd9-4b1e-93be-7d8a6c5a9b7b
# ╠═3a76057d-5acd-4f4e-a319-f4bfb6bcf9c2
# ╠═829dcb3f-112a-496b-ae99-94322920c2e4
