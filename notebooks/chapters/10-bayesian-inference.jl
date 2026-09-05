### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ b705b5b7-4e29-433c-9703-58867199cacf
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ cbda1dd5-a97f-44a7-9283-99ae64918792
begin
    using OrdinaryDiffEq
    using SciMLSensitivity
    using Turing
    using MCMCChains
    using Distributions
    using StableRNGs
    using CairoMakie
    using PlutoUI
end

# ╔═╡ 107678bc-a791-40e9-aac8-6c3da0e891cf
md"""
# Bayesian inference

Companion notebook for the `ch10-bayesian-inference` unit. NUTS is run **once**
with lightweight, reproducible settings (300 samples × 2 chains). The sliders
then explore the posterior *predictive* distribution without re-sampling: change
the credible level and the number of posterior draws used to build the band.
"""

# ╔═╡ b433198b-1367-441d-a3a7-e0f7e5339da9
function sir_ode!(du, u, p, t)
    S, I, R, C = u
    β, c, γ = p
    N = S + I + R
    infection = β * c * I / N * S
    du[1] = -infection
    du[2] = infection - γ * I
    du[3] = γ * I
    du[4] = infection
end

# ╔═╡ 49deb993-b7ec-4a07-91d7-333e3f93c32a
begin
    tmax = 40.0
    tspan = (0.0, tmax)
    obstimes = 1.0:1.0:tmax
    u0 = [990.0, 10.0, 0.0, 0.0]
    prob_ode = ODEProblem(sir_ode!, u0, tspan, [0.05, 10.0, 0.25])
    sol_true = solve(prob_ode, Tsit5(); saveat = 1.0)
    incidence_true = Array(sol_true)[4, 2:end] .- Array(sol_true)[4, 1:(end - 1)]
    data = rand.(StableRNG(20240607), Poisson.(incidence_true))
end;

# ╔═╡ cc175dda-ea4f-4a39-9e4f-7cffda832a49
@model function bayes_sir(y, prob)
    l = length(y)
    i₀ ~ Beta(1.0, 20.0)
    β ~ Beta(1.0, 10.0)
    I = i₀ * 1000.0
    prob = remake(prob; u0 = [1000.0 - I, I, 0.0, 0.0], p = [β, 10.0, 0.25],
        tspan = (0.0, float(l)))
    sol = solve(prob, Tsit5(); saveat = 1.0)
    if !SciMLBase.successful_retcode(sol)
        Turing.@addlogprob! -Inf
        return nothing
    end
    C = Array(sol)[4, :]
    X = C[2:end] .- C[1:(end - 1)]
    for i in 1:l
        y[i] ~ Poisson(max(X[i], 0.0) + 1e-6)
    end
    return nothing
end

# ╔═╡ 92531747-9e27-4a25-94ae-715563b32126
chain = sample(StableRNG(2), bayes_sir(data, prob_ode), NUTS(0.8),
    MCMCSerial(), 300, 2; chain_type = MCMCChains.Chains, progress = false);

# ╔═╡ 6365306f-c354-4891-aee4-da501b2d1ca1
begin
    i0s = vec(Array(chain[:, :i₀, :]))
    βs = vec(Array(chain[:, :β, :]))
    diagnostics = (
        max_rhat = round(maximum(rhat(chain).nt.rhat); digits = 3),
        min_ess = round(minimum(ess(chain).nt.ess); digits = 1),
        posterior_mean = (i₀ = round(mean(i0s); digits = 4),
            β = round(mean(βs); digits = 4)),
        truth = (i₀ = 0.01, β = 0.05),
    )
end

# ╔═╡ 38dbba5b-8bb8-4d5c-90a3-7d78dc17c132
md"""**Credible level for the posterior predictive band:**"""

# ╔═╡ 5d8cdd6b-f346-4fd5-ba6b-c7987d98ea38
@bind cred_level PlutoUI.Slider(0.50:0.05:0.99; default = 0.95, show_value = true)

# ╔═╡ 7819a12f-124e-46ac-a1e1-4560cfc5ca4c
md"""**Number of posterior draws used for the predictive band:**"""

# ╔═╡ 6351b00e-eeec-4fc7-b2e9-0c7fea55dc38
@bind n_draws PlutoUI.Slider(50:25:200; default = 150, show_value = true)

# ╔═╡ 9ddab425-1a67-4706-9c44-cefccb3806a0
begin
    function predict_incidence(i₀, β, l)
        I = i₀ * 1000.0
        solp = solve(remake(prob_ode; u0 = [1000.0 - I, I, 0.0, 0.0],
            p = [β, 10.0, 0.25], tspan = (0.0, float(l))), Tsit5(); saveat = 1.0)
        C = Array(solp)[4, :]
        C[2:end] .- C[1:(end - 1)]
    end
    l = length(data)
    idx = rand(StableRNG(3), 1:length(βs), n_draws)
    pred_curves = reduce(hcat, [predict_incidence(i0s[k], βs[k], l) for k in idx])
end;

# ╔═╡ 17f28ab9-b574-44c1-b7bc-ab6bc917e78b
let
    α = (1 - cred_level) / 2
    lo = [quantile(pred_curves[t, :], α) for t in 1:length(data)]
    hi = [quantile(pred_curves[t, :], 1 - α) for t in 1:length(data)]
    md = [quantile(pred_curves[t, :], 0.5) for t in 1:length(data)]
    fig = Figure(size = (760, 420))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "new cases per day",
        title = "posterior predictive ($(Int(round(100cred_level)))% band, $(n_draws) draws)")
    band!(ax, collect(obstimes), lo, hi; color = (:steelblue, 0.3),
        label = "credible band")
    lines!(ax, collect(obstimes), md; color = :steelblue, linewidth = 3,
        label = "posterior median")
    scatter!(ax, collect(obstimes), data; color = :grey20, label = "data")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 64655d77-054b-4eae-9f06-c5901685f9bb
md"""
Widening the credible level or increasing the number of draws smooths and widens
the band. Because the data are informative, even a modest number of draws
produces a band that tracks the observations. The book chapter adds trace plots
and marginal densities as convergence diagnostics.
"""

# ╔═╡ ff1c100b-b3dc-4075-ba86-7a12462b1a85
(
    course_unit = "ch10-bayesian-inference",
    status = "complete",
    cred_level = cred_level,
    n_draws = n_draws,
    diagnostics = diagnostics,
)

# ╔═╡ Cell order:
# ╠═b705b5b7-4e29-433c-9703-58867199cacf
# ╠═cbda1dd5-a97f-44a7-9283-99ae64918792
# ╟─107678bc-a791-40e9-aac8-6c3da0e891cf
# ╠═b433198b-1367-441d-a3a7-e0f7e5339da9
# ╠═49deb993-b7ec-4a07-91d7-333e3f93c32a
# ╠═cc175dda-ea4f-4a39-9e4f-7cffda832a49
# ╠═92531747-9e27-4a25-94ae-715563b32126
# ╠═6365306f-c354-4891-aee4-da501b2d1ca1
# ╟─38dbba5b-8bb8-4d5c-90a3-7d78dc17c132
# ╠═5d8cdd6b-f346-4fd5-ba6b-c7987d98ea38
# ╟─7819a12f-124e-46ac-a1e1-4560cfc5ca4c
# ╠═6351b00e-eeec-4fc7-b2e9-0c7fea55dc38
# ╠═9ddab425-1a67-4706-9c44-cefccb3806a0
# ╠═17f28ab9-b574-44c1-b7bc-ab6bc917e78b
# ╟─64655d77-054b-4eae-9f06-c5901685f9bb
# ╠═ff1c100b-b3dc-4075-ba86-7a12462b1a85
