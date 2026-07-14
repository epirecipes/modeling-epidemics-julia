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

# ╔═╡ e29bd1d0-fd18-47ca-9d85-d3b31ca4eaf5
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "environments", "ch08"));
        io = devnull)
end

# ╔═╡ d64c31c9-3f09-4a12-8ab9-1575172f2790
begin
    using OrdinaryDiffEq
    using ProfileLikelihood
    using Optimization
    using OptimizationOptimJL
    using Distributions
    using StableRNGs
    using CairoMakie
    using PlutoUI
end

# ╔═╡ de64ee5f-424e-4793-ad03-612fad9f3ef9
md"""
# Frequentist parameter estimation

Companion notebook for the `ch08-parameter-estimation` unit. It runs in the
isolated `environments/ch08` project. The maximum-likelihood estimate (MLE) is
computed once from simulated Poisson case reports; the sliders let you resimulate
the data and slide a **candidate** transmission probability `β` to feel the shape
of the likelihood around the optimum.
"""

# ╔═╡ 1834055d-be1a-4ee3-a18f-712c76e823ba
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

# ╔═╡ 1df98c0f-a7c4-4e33-b5ed-966bec8a9ec3
begin
    δt = 1.0
    tmax = 40.0
    tspan = (0.0, tmax)
    u0 = [990.0, 10.0, 0.0, 0.0]
    p_true = [0.05, 10.0, 0.25]      # β, c, γ; truth i₀ = 0.01, β = 0.05
    obstimes = δt:δt:tmax
    prob_ode = ODEProblem(sir_ode!, u0, tspan, p_true)
end;

# ╔═╡ 67e8b590-f765-4789-80a0-f46288d1ddd2
md"""**Data seed** (resimulate the Poisson case reports):"""

# ╔═╡ bd3d1bb2-9cde-4b25-a274-d82c2c6b07bc
@bind seed PlutoUI.Slider(1:50; default = 7, show_value = true)

# ╔═╡ 1b8830e4-a096-4023-9195-bc217eeafeb7
begin
    sol_true = solve(prob_ode, Tsit5(); saveat = δt)
    incidence_true = Array(sol_true)[4, 2:end] .- Array(sol_true)[4, 1:(end - 1)]
    data = rand.(StableRNG(seed), Poisson.(incidence_true))
end;

# ╔═╡ 7cfb281b-ea49-4db2-9b58-b6dc4226cd1e
begin
    function incidence_at(θ)
        i₀, β = θ
        I = i₀ * 1000.0
        solf = solve(remake(prob_ode; u0 = [1000.0 - I, I, 0.0, 0.0],
            p = [β, 10.0, 0.25]), Tsit5(); saveat = δt)
        C = Array(solf)[4, :]
        C[2:end] .- C[1:(end - 1)]
    end
    function poisson_loglik(θ)
        X = incidence_at(θ)
        any(X .<= 0) && return -Inf
        sum(logpdf.(Poisson.(X), data))
    end
end

# ╔═╡ 97b0178b-4fd3-471a-84e2-7f90d4cc4b40
begin
    function loglik(θ, data, integrator)
        i₀, β = θ
        integrator.p[1] = β; integrator.p[2] = 10.0; integrator.p[3] = 0.25
        I = i₀ * 1000.0
        reinit!(integrator, [1000.0 - I, I, 0.0, 0.0])
        solve!(integrator)
        C = Array(integrator.sol)[4, :]
        X = C[2:end] .- C[1:(end - 1)]
        any(X .<= 0) && return -Inf
        sum(logpdf.(Poisson.(X), data))
    end
    prob = LikelihoodProblem(loglik, [0.01, 0.1], sir_ode!, u0, tspan;
        syms = [:i₀, :β], data = data, ode_parameters = p_true,
        ode_kwargs = (verbose = false, saveat = δt),
        f_kwargs = (adtype = Optimization.AutoFiniteDiff(),),
        prob_kwargs = (lb = [0.0, 0.0], ub = [1.0, 1.0]), ode_alg = Tsit5())
    θ̂ = get_mle(mle(prob, NelderMead()))
end

# ╔═╡ b2674dd3-ee99-42d5-8451-6e549fd5d259
md"""**Candidate β** (compare against the data and the MLE fit):"""

# ╔═╡ 75e6fbfa-2b1e-4f13-9c05-564443614ec0
@bind β_candidate PlutoUI.Slider(0.03:0.001:0.07; default = 0.05, show_value = true)

# ╔═╡ 4a762770-f352-4320-9b3d-fd338aaf5dd5
let
    cand = [θ̂[1], β_candidate]
    ll = round(poisson_loglik(cand); digits = 1)
    fig = Figure(size = (760, 420))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "new cases per day",
        title = "candidate β = $(β_candidate), log-likelihood = $(ll)")
    barplot!(ax, collect(obstimes), data; color = (:grey30, 0.35), label = "data")
    lines!(ax, collect(obstimes), incidence_at([0.01, 0.05]); color = :firebrick,
        linewidth = 3, label = "truth")
    lines!(ax, collect(obstimes), incidence_at(θ̂); color = :steelblue,
        linewidth = 3, linestyle = :dash, label = "MLE")
    lines!(ax, collect(obstimes), incidence_at(cand); color = :darkorange,
        linewidth = 2, label = "candidate")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ b3b1cf9a-de7e-46e2-b043-5c788b6172fa
md"""
Slide the candidate `β` away from the MLE and the log-likelihood drops: the
orange curve pulls away from both the data and the blue MLE fit. The likelihood
is maximised near the truth `β = 0.05`. The book chapter turns this curvature
into profile-likelihood confidence intervals via Wilks' theorem.
"""

# ╔═╡ 4af6682f-0032-4756-8d18-bce2a03b7290
(
    course_unit = "ch08-parameter-estimation",
    status = "complete",
    seed = seed,
    mle = (i₀ = round(θ̂[1]; digits = 4), β = round(θ̂[2]; digits = 4)),
    candidate_β = β_candidate,
    candidate_loglik = round(poisson_loglik([θ̂[1], β_candidate]); digits = 1),
    mle_loglik = round(poisson_loglik(θ̂); digits = 1),
)

# ╔═╡ Cell order:
# ╠═e29bd1d0-fd18-47ca-9d85-d3b31ca4eaf5
# ╠═d64c31c9-3f09-4a12-8ab9-1575172f2790
# ╟─de64ee5f-424e-4793-ad03-612fad9f3ef9
# ╠═1834055d-be1a-4ee3-a18f-712c76e823ba
# ╠═1df98c0f-a7c4-4e33-b5ed-966bec8a9ec3
# ╟─67e8b590-f765-4789-80a0-f46288d1ddd2
# ╠═bd3d1bb2-9cde-4b25-a274-d82c2c6b07bc
# ╠═1b8830e4-a096-4023-9195-bc217eeafeb7
# ╠═7cfb281b-ea49-4db2-9b58-b6dc4226cd1e
# ╠═97b0178b-4fd3-471a-84e2-7f90d4cc4b40
# ╟─b2674dd3-ee99-42d5-8451-6e549fd5d259
# ╠═75e6fbfa-2b1e-4f13-9c05-564443614ec0
# ╠═4a762770-f352-4320-9b3d-fd338aaf5dd5
# ╟─b3b1cf9a-de7e-46e2-b043-5c788b6172fa
# ╠═4af6682f-0032-4756-8d18-bce2a03b7290
