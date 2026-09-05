### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 06c06160-80a6-41fe-93f1-cc82361f57a7
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ b30f92a1-d60a-43f3-b5c8-7b87897d4aac
begin
    using CairoMakie
    using OrdinaryDiffEq
    using PlutoUI
end

# ╔═╡ 74ba9e9c-9596-4e3e-a6cd-2c74ce088117
md"""
# Sensitivity and identifiability

This notebook explores how the infectious peak changes over a parameter range.
The finite-difference derivative is local; the heat map is a small global
range experiment. Both use deterministic ODE solves.

Raw derivatives are plotted beside dimensionless elasticities: the fractional
change in I per fractional change in a parameter. Wherever they are nonzero the
raw β and c curves are a factor c/β = 200 apart, which is a units effect rather
than a difference in influence, since the model uses only the product βc.
Rescaling removes it and the two curves then coincide exactly.
"""

# ╔═╡ 09c00763-7c1b-4af4-8488-3679ccbd8806
@bind sensitivity_span PlutoUI.Slider(5:5:30; default = 15, show_value = true)

# ╔═╡ c2df82d0-6b65-4b7c-86dd-c10d399fc86f
begin
    using EpiModelingCourse

    const N10 = BASELINE_PARAMETERS.N
    const u010 = [BASELINE_INITIAL_STATE.S, BASELINE_INITIAL_STATE.I,
        BASELINE_INITIAL_STATE.R]
    function sir10!(du, u, p, t)
        β, c, γ = p
        incidence = β * c * u[1] * u[2] / N10
        du[1] = -incidence
        du[2] = incidence - γ * u[2]
        du[3] = γ * u[2]
    end
    times10 = collect(0.0:1.0:40.0)
    p010 = [BASELINE_PARAMETERS.beta, BASELINE_PARAMETERS.contact_rate,
        BASELINE_PARAMETERS.gamma]
end

# ╔═╡ 1f6f76a7-1d1c-4d9a-8d80-69f9502e188d
begin
    function trajectory10(p)
        Array(solve(ODEProblem(sir10!, u010, (0.0, 40.0), p), Tsit5();
            saveat = times10, abstol = 1e-8, reltol = 1e-8))
    end

    # Dense interpolation on a sub-daily grid, since the peak falls between days.
    function peak_infectious10(p)
        solution = solve(ODEProblem(sir10!, u010, (0.0, 40.0), p), Tsit5();
            abstol = 1e-8, reltol = 1e-8)
        maximum(solution(t)[2] for t in 0.0:0.05:40.0)
    end

    base10 = trajectory10(p010)
    local_sensitivity10 = Matrix{Float64}(undef, length(times10), 3)
    for j in 1:3
        δ = 1e-5 * p010[j]
        plus = copy(p010)
        minus = copy(p010)
        plus[j] += δ
        minus[j] -= δ
        local_sensitivity10[:, j] =
            (trajectory10(plus)[2, :] - trajectory10(minus)[2, :]) ./ (2δ)
    end
    # Elasticities: fractional change in I per fractional change in a parameter.
    elasticity10 = hcat(
        (p010[j] .* local_sensitivity10[:, j] ./ base10[2, :] for j in 1:3)...,
    )
end

# ╔═╡ 8c4d7c0c-bf3f-4d4d-94bb-335cf994a6cb
begin
    fig_local10 = Figure(size = (960, 400))
    ax_local10 = Axis(fig_local10[1, 1], xlabel = "time",
        ylabel = "∂I/∂θⱼ", title = "Raw derivatives (unit-dependent)")
    ax_elastic10 = Axis(fig_local10[1, 2], xlabel = "time",
        ylabel = "(θⱼ/I) ∂I/∂θⱼ", title = "Elasticities (dimensionless)")
    for (j, label) in enumerate(["β", "c", "γ"])
        style = j == 2 ? :dash : :solid
        lines!(ax_local10, times10, local_sensitivity10[:, j], label = label,
            linestyle = style)
        lines!(ax_elastic10, times10, elasticity10[:, j], label = label,
            linestyle = style)
    end
    axislegend(ax_elastic10)
    fig_local10
end

# ╔═╡ 6c2a0b89-c9a5-4de5-8ad2-2cce6ec86b8e
begin
    beta_range10 = range(p010[1] * (1 - sensitivity_span / 100),
        p010[1] * (1 + sensitivity_span / 100); length = 7)
    gamma_range10 = range(p010[3] * (1 - sensitivity_span / 100),
        p010[3] * (1 + sensitivity_span / 100); length = 7)
    peak_surface10 = [
        peak_infectious10([β, p010[2], γ])
        for γ in gamma_range10, β in beta_range10
    ]
    fig_surface10 = Figure(size = (800, 440))
    ax_surface10 = Axis(fig_surface10[1, 1], xlabel = "β", ylabel = "γ",
        title = "Peak I over ±$(sensitivity_span)% ranges")
    heatmap!(ax_surface10, collect(beta_range10), collect(gamma_range10),
        peak_surface10)
    fig_surface10
end

# ╔═╡ a736e6cf-7d16-4d4c-9df3-a980f60f4040
begin
    conservation_error10 = maximum(abs.(vec(sum(base10; dims = 1)) .- N10))
    ε10 = 1e-5
    p_plus10 = copy(p010)
    p_plus10[1] += ε10
    finite_difference10 = (trajectory10(p_plus10)[2, end] - base10[2, end]) / ε10
    derivative_error10 = abs(finite_difference10 -
        local_sensitivity10[end, 1]) / (1 + abs(local_sensitivity10[end, 1]))
    elasticity_gap10 = maximum(abs.(elasticity10[:, 1] .- elasticity10[:, 2]))
    @assert conservation_error10 < 1e-6
    @assert derivative_error10 < 1e-3
    @assert elasticity_gap10 < 1e-6
    (
        course_unit = "ch12-sensitivity-and-identifiability",
        status = "complete",
        sensitivity_span_percent = sensitivity_span,
        conservation_error = conservation_error10,
        finite_difference_relative_error = derivative_error10,
        elasticity_gap = elasticity_gap10,
        identifiability_note = "incidence identifies β*c, not β and c separately",
    )
end

# ╔═╡ Cell order:
# ╠═06c06160-80a6-41fe-93f1-cc82361f57a7
# ╠═b30f92a1-d60a-43f3-b5c8-7b87897d4aac
# ╟─74ba9e9c-9596-4e3e-a6cd-2c74ce088117
# ╠═09c00763-7c1b-4af4-8488-3679ccbd8806
# ╠═c2df82d0-6b65-4b7c-86dd-c10d399fc86f
# ╠═1f6f76a7-1d1c-4d9a-8d80-69f9502e188d
# ╠═8c4d7c0c-bf3f-4d4d-94bb-335cf994a6cb
# ╠═6c2a0b89-c9a5-4de5-8ad2-2cce6ec86b8e
# ╠═a736e6cf-7d16-4d4c-9df3-a980f60f4040
