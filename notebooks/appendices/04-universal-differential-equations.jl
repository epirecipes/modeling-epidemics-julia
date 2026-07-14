### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 0f4a28cd-91e9-4ebb-9c89-000d7c384f2c
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ bb428d14-7f87-47d2-8d95-20a4568e132e
begin
    using OrdinaryDiffEq
    using Distributions
    using StableRNGs
    using Random
    using Lux
    using DiffEqFlux
    using SciMLSensitivity
    using Optimization
    using OptimizationOptimisers
    using OptimizationOptimJL
    using CairoMakie
    using PlutoUI
end

# ╔═╡ f4fa7a60-9969-4ba1-9e16-43592d421cec
md"""
# Universal differential equations

Companion notebook for the `a04-universal-differential-equations` unit. We keep
the SIR bookkeeping and replace only the **force of infection** with a small
neural network, then fit the hybrid model by differentiating through the ODE
solver.

The network is trained **once** (a deliberately small, seeded run so the
notebook stays responsive). The sliders below are then *instant*: they explore
the learned closure and, crucially, its **extrapolation region** beyond the
range of infectious proportions the training data actually visited.
"""

# ╔═╡ f240f960-27e5-485b-bba2-87fb63a9509f
md"""
The data-generating model closes the force of infection with mass action,
``\lambda^\star(I) = \beta c\, I = 0.5\,I`` (course baseline in proportions),
with ``\gamma = 0.25``. Noisy daily incidence is generated from a fixed
`StableRNG`.
"""

# ╔═╡ 626a5341-43c7-430b-a5ba-422cb755d684
begin
    N = 1000.0
    βc = 0.5
    γ = 0.25
    u0 = [0.99, 0.01, 0.0]        # S, I, C as proportions
    tspan = (0.0, 40.0)
    δt = 1.0
    train_time = 30.0

    function sir_true!(du, u, p, t)
        S, I, C = u
        λ, ν = p
        du[1] = -λ * S * I
        du[2] = λ * S * I - ν * I
        du[3] = λ * S * I
    end

    sol_true = solve(ODEProblem(sir_true!, u0, tspan, [βc, γ]), Tsit5(); saveat = δt)
    ts = 0.0:δt:train_time
    Ctrue = reduce(hcat, sol_true(ts).u)[3, :]
    data_rng = StableRNG(20240501)
    noisy = [rand(data_rng, Poisson(max(c, 0.0))) for c in N .* diff(Ctrue)]
    Imax_train = maximum(reduce(hcat, sol_true(ts).u)[2, :])
    (observations = length(noisy), Imax_train = round(Imax_train; digits = 4))
end

# ╔═╡ 2d28507d-45e5-4b68-be0c-87ee89dc24a8
md"""
A `Lux.Chain` with one hidden layer of five `tanh` units and a `softplus`
output (so ``\lambda_\theta \ge 0`` by construction) — **16 parameters**. We
flatten/unflatten the parameters by hand for this fixed architecture.
"""

# ╔═╡ e44c9c39-d477-42a8-a208-06182c0195a2
begin
    foi = Lux.Chain(Lux.Dense(1, 5, tanh), Lux.Dense(5, 1, softplus))
    _, foi_state = Lux.setup(StableRNG(1), foi)

    function params_from_vector(θ)
        return (layer_1 = (weight = reshape(θ[1:5], 5, 1), bias = θ[6:10]),
                layer_2 = (weight = reshape(θ[11:15], 1, 5), bias = θ[16:16]))
    end

    force_of_infection(I, θ) = foi([I], params_from_vector(θ), foi_state)[1][1]
    p_init = 0.1 .* randn(StableRNG(2), Float64, 16)
    length(p_init)
end

# ╔═╡ a88fdb66-1bd4-4e4b-ae44-d450b6a65981
md"""
### Train once (fixed, not tied to the sliders)

`Adam` makes fast early progress; `BFGS` polishes. Gradients flow through the
ODE solve via an `InterpolatingAdjoint` with a reverse-mode vector–Jacobian
product. This cell runs a single time; the sliders further down do **not**
retrain it.
"""

# ╔═╡ 9ac30325-540e-46fb-bd30-1fee6a58c692
θ_hat = let
    function sir_ude!(du, u, θ, t)
        S, I, C = u
        λ = force_of_infection(I, θ)
        du[1] = -λ * S
        du[2] = λ * S - γ * I
        du[3] = λ * S
    end
    prob_ude = ODEProblem(sir_ude!, u0, (0.0, train_time), p_init)
    predict(θ) = Array(solve(prob_ude, Tsit5(); p = θ, saveat = δt,
        sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP())))
    loss(θ) = (μ = max.(N .* diff(predict(θ)[3, :]), 1e-6); sum(μ .- noisy .* log.(μ)))

    optf = OptimizationFunction((θ, _) -> loss(θ), Optimization.AutoZygote())
    res_adam = solve(OptimizationProblem(optf, p_init),
        OptimizationOptimisers.Adam(0.05); maxiters = 120)
    res = solve(OptimizationProblem(optf, res_adam.u),
        OptimizationOptimJL.BFGS(); maxiters = 40)
    res.u
end

# ╔═╡ c7d7815b-fad2-4be9-8382-898ec19f98c9
md"""
### Explore the learned closure

**Display range for the infectious proportion ``I``** (how far past the
training data to extrapolate):
"""

# ╔═╡ 15e11195-072e-4b69-a420-9df087722857
@bind view_Imax PlutoUI.Slider(0.1:0.05:0.5; default = 0.35, show_value = true)

# ╔═╡ 2dc6d6c2-b4d5-49c1-864b-13a32fd1e142
md"""**Evaluate the network at a single infectious proportion ``I``:**"""

# ╔═╡ 90269c98-cc1e-4b33-801e-61409b6fe612
@bind eval_I PlutoUI.Slider(0.0:0.005:0.5; default = 0.1, show_value = true)

# ╔═╡ b484ae43-15b8-4866-b862-8f6be4f8a0d1
begin
    λ_nn_here = force_of_infection(eval_I, θ_hat)
    λ_true_here = βc * eval_I
    in_training_range = eval_I <= Imax_train
    (I = eval_I, λ_network = round(λ_nn_here; digits = 4),
     λ_true = round(λ_true_here; digits = 4),
     abs_error = round(abs(λ_nn_here - λ_true_here); sigdigits = 3),
     within_training_range = in_training_range)
end

# ╔═╡ 35beb4d3-cb27-471c-b949-5999c7e6967d
begin
    Igrid = 0.0:0.005:view_Imax
    λ_nn = [force_of_infection(I, θ_hat) for I in Igrid]
    λ_true = βc .* Igrid
    inr = Igrid .<= Imax_train
    max_abs_error = maximum(abs.(λ_nn[inr] .- λ_true[inr]))

    fig = Figure(size = (760, 430))
    ax = Axis(fig[1, 1]; xlabel = "infectious proportion I",
        ylabel = "force of infection λ",
        title = "Learned closure (max error in training range ≈ " *
                string(round(max_abs_error; sigdigits = 2)) * ")")
    view_Imax > Imax_train && vspan!(ax, Imax_train, view_Imax; color = (:gray, 0.12))
    lines!(ax, Igrid, λ_true; color = :black, linewidth = 2, label = "true βc·I")
    lines!(ax, Igrid, λ_nn; color = :firebrick, linestyle = :dash, linewidth = 2,
        label = "neural network")
    vlines!(ax, [Imax_train]; color = :gray, linestyle = :dot)
    scatter!(ax, [eval_I], [force_of_infection(eval_I, θ_hat)]; color = :firebrick,
        markersize = 12)
    axislegend(ax; position = :lt)
    fig
end

# ╔═╡ d8763c86-1845-4a45-8e5e-1fd0af76a955
md"""
**What to notice.** The network matches ``\lambda^\star(I) = \beta c\,I`` closely
to the *left* of the dotted line (the largest infectious proportion in the
training data, ``I_{\max}\approx`` $(round(Imax_train; digits=3))). In the shaded
region to the *right* it is unconstrained — push the display slider past the
dotted line and the dashed curve drifts away from the truth. Identifiability is
local to the data.
"""

# ╔═╡ b3410ca9-29b2-46cf-83a6-a2da4e0ed309
begin
    @assert max_abs_error < 0.05 "closure should match the truth within the training range"
    (
        course_unit = "a04-universal-differential-equations",
        status = "complete",
        controls = (; view_Imax, eval_I),
        max_abs_error_in_range = round(max_abs_error; sigdigits = 3),
    )
end

# ╔═╡ Cell order:
# ╠═0f4a28cd-91e9-4ebb-9c89-000d7c384f2c
# ╠═bb428d14-7f87-47d2-8d95-20a4568e132e
# ╟─f4fa7a60-9969-4ba1-9e16-43592d421cec
# ╟─f240f960-27e5-485b-bba2-87fb63a9509f
# ╠═626a5341-43c7-430b-a5ba-422cb755d684
# ╟─2d28507d-45e5-4b68-be0c-87ee89dc24a8
# ╠═e44c9c39-d477-42a8-a208-06182c0195a2
# ╟─a88fdb66-1bd4-4e4b-ae44-d450b6a65981
# ╠═9ac30325-540e-46fb-bd30-1fee6a58c692
# ╟─c7d7815b-fad2-4be9-8382-898ec19f98c9
# ╠═15e11195-072e-4b69-a420-9df087722857
# ╟─2dc6d6c2-b4d5-49c1-864b-13a32fd1e142
# ╠═90269c98-cc1e-4b33-801e-61409b6fe612
# ╠═b484ae43-15b8-4866-b862-8f6be4f8a0d1
# ╠═35beb4d3-cb27-471c-b949-5999c7e6967d
# ╟─d8763c86-1845-4a45-8e5e-1fd0af76a955
# ╠═b3410ca9-29b2-46cf-83a6-a2da4e0ed309
