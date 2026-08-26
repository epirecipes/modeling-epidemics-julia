### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 9b0c4228-97aa-4779-9a5b-f71cb693e2c1
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 56da20cc-06d7-4f5d-9b44-1d53070b90b0
begin
    using DelayDiffEq
    using OrdinaryDiffEq
    using StochasticDiffEq
    using SparseArrays
    using Random
    using CairoMakie
    using PlutoUI
end

# ╔═╡ 894e7c1d-7d7a-4c38-b130-d4ef01caafb8
md"""
# Stochastic delay differential equations

Companion notebook for `a02-stochastic-delay-differential-equations` (see the
[book chapter](../../appendices/02-stochastic-delay-differential-equations.qmd)
for the full derivation).

**This unit's software combination is experimental**, and this notebook is
built to make that concrete rather than to hide it: change the seed below and
watch the minimum state swing from safely nonnegative to substantially
negative, with conservation holding regardless. Toggle the step size to see
that finer steps do **not** monotonically improve nonnegativity the way they
would for a convergent numerical method. Toggle the post-hoc clamp to see
that "fixing" negative states by clamping the *current* value breaks exact
conservation, because the delayed recovery term keeps depending on the
uncorrected history.
"""

# ╔═╡ a395554d-f61f-4255-834c-3b8714e6befe
md"""
Random seed (passed to the **global** `Random.seed!`, not `rng=`, per the
book chapter's caveat): $(@bind seed PlutoUI.Slider(1:10; default = 4, show_value = true))

Fixed step size `dt` (adaptive stepping fails at the history discontinuity,
so this is fixed, not tuned for accuracy): $(@bind dt Select([0.05, 0.1, 0.2, 0.4]; default = 0.1))

Apply a post-hoc `DiscreteCallback` that clamps I to zero when negative:
$(@bind use_clamp CheckBox(default = false))
"""

# ╔═╡ 9c7f18a1-431b-4510-9f7d-dea2c6eedd2b
function sir_dde!(du, u, h, p, t)
    S, I, R = u
    β, c, τ = p
    N = S + I + R
    λ = β * c * I / N
    Sd, Id, Rd = h(p, t - τ)
    Nd = Sd + Id + Rd
    λd = β * c * Id / Nd
    du[1] = -λ * S
    du[2] = λ * S - λd * Sd
    du[3] = λd * Sd
    return nothing
end

# ╔═╡ 09e1ada8-832a-44a9-91cd-d4f14b2faabd
function sir_delayed_noise!(G, u, h, p, t)
    S, I, R = u
    β, c, τ = p
    N = S + I + R
    infection = max(β * c * I / N * S, 0.0)
    Sd, Id, Rd = h(p, t - τ)
    Nd = Sd + Id + Rd
    recovery = max(β * c * Id / Nd * Sd, 0.0)
    G[1, 1] = -sqrt(infection)
    G[2, 1] = sqrt(infection)
    G[2, 2] = -sqrt(recovery)
    G[3, 2] = sqrt(recovery)
    return nothing
end

# ╔═╡ 2e2d2023-7d79-4612-bfe1-37e26cd5639f
begin
    sir_history(p, t) = [1000.0, 0.0, 0.0]
    noise_shape() = sparse([1.0 0.0; 1.0 1.0; 0.0 1.0])
    u0 = [990.0, 10.0, 0.0]
    γ = 0.25
    p_sdde = (0.05, 10.0, 1 / γ)
    tspan = (0.0, 40.0)
    nothing
end

# ╔═╡ a8129575-b19a-4c0c-b378-c328a6a8fec7
begin
    fire_count = Ref(0)
    function clamp_condition(u, t, integrator)
        return u[2] < 0
    end
    function clamp_affect!(integrator)
        fire_count[] += 1
        integrator.u[2] = max(integrator.u[2], 0.0)
        return nothing
    end
    clamp_cb = DiscreteCallback(clamp_condition, clamp_affect!)

    Random.seed!(seed)   # global RNG seeding: the documented workaround for this combination
    prob_sdde = SDDEProblem(sir_dde!, sir_delayed_noise!, u0, sir_history, tspan, p_sdde;
                             noise_rate_prototype = noise_shape())
    sol_sdde = if use_clamp
        solve(prob_sdde, LambaEM(); dt = dt, adaptive = false, callback = clamp_cb)
    else
        solve(prob_sdde, LambaEM(); dt = dt, adaptive = false)
    end

    sol_dde = solve(
        DDEProblem(sir_dde!, u0, sir_history, tspan, p_sdde; constant_lags = [p_sdde[3]]),
        MethodOfSteps(Tsit5());
        reltol = 1e-8, abstol = 1e-8,
    )
    states_sdde = Array(sol_sdde)
    conservation_error = maximum(abs.(vec(sum(states_sdde; dims = 1)) .- 1000.0))
    min_state = minimum(states_sdde)
    nothing
end

# ╔═╡ 29bacdbd-451f-4f71-91d6-54ae785d1202
begin
    grid = tspan[1]:dt:tspan[2]
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "number",
              title = "SDDE realization (seed $(seed), dt = $(dt)) vs deterministic DDE drift")
    lines!(ax, grid, [sol_sdde(t)[1] for t in grid]; color = :steelblue, linewidth = 3, label = "S (SDDE)")
    lines!(ax, grid, [sol_sdde(t)[2] for t in grid]; color = :firebrick, linewidth = 3, label = "I (SDDE)")
    lines!(ax, grid, [sol_sdde(t)[3] for t in grid]; color = :seagreen, linewidth = 3, label = "R (SDDE)")
    lines!(ax, grid, [sol_dde(t)[2] for t in grid]; color = :firebrick, linewidth = 2, linestyle = :dash,
           label = "I (DDE drift)")
    hlines!(ax, [0.0]; color = :black, linestyle = :dot)
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ a78a0dea-de26-44c3-8522-fdde4e48bd0d
md"""
**Current run** (seed = $(seed), dt = $(dt), clamp callback $(use_clamp ? "ON" : "OFF")):

- Maximum conservation error over the full trajectory: **$(round(conservation_error; sigdigits = 3))**
- Minimum state value reached (any compartment): **$(round(min_state; digits = 2))**
- Nonnegative throughout: **$(min_state >= -1e-6)**
$(use_clamp ? "- Clamp callback fired **$(fire_count[])** times" : "")

$(min_state < -1 && !use_clamp ? "This seed/step combination drives a compartment **substantially negative** — exactly the reproducible caveat the book chapter discusses: conservation is structural (built into the noise matrix), but nonnegativity is not." : "")

$(use_clamp && conservation_error > 1e-4 ? "With the clamp callback on, conservation is now **violated** (error ≫ solver tolerance): clamping the *current* state cannot repair the *delayed* recovery term, which keeps depending on the pre-clamp history for up to τ more time units. This is the concrete demonstration behind Exercise 2 in the book chapter." : "")
"""

# ╔═╡ 0dd40cbb-74e3-49a2-bd24-793efa4fb6ea
md"""
## Why this is documented as experimental

1. **Adaptive stepping is currently unusable** for `SDDEProblem` at this
   history discontinuity — the step-size controller collapses to machine
   epsilon almost immediately, so every run above uses a fixed step.
2. **The `rng=` keyword is unreliable** for `SDDEProblem`/`MethodOfSteps`;
   reproducibility currently requires seeding the **global** RNG with
   `Random.seed!` before each `solve` call, unlike the rest of the course.
3. **Nonnegativity is not guaranteed**, and — as the clamp toggle above
   demonstrates — a per-step correction is a much more incomplete fix here
   than for a non-delayed chemical Langevin SDE, because it cannot repair
   the delayed drift/diffusion terms that keep referencing the uncorrected
   history.

None of this invalidates the chemical Langevin approximation itself; it
means this particular package combination is, at present, best treated as a
teaching device rather than a production tool.
"""

# ╔═╡ 8c8e3b0e-9281-421d-ad65-9a875ec3416b
(
    course_unit = "a02-stochastic-delay-differential-equations",
    status = "complete",
    seed = seed,
    dt = dt,
    use_clamp = use_clamp,
    conservation_error = round(conservation_error; sigdigits = 3),
    min_state = round(min_state; digits = 3),
)

# ╔═╡ Cell order:
# ╠═9b0c4228-97aa-4779-9a5b-f71cb693e2c1
# ╠═56da20cc-06d7-4f5d-9b44-1d53070b90b0
# ╟─894e7c1d-7d7a-4c38-b130-d4ef01caafb8
# ╟─a395554d-f61f-4255-834c-3b8714e6befe
# ╠═9c7f18a1-431b-4510-9f7d-dea2c6eedd2b
# ╠═09e1ada8-832a-44a9-91cd-d4f14b2faabd
# ╠═2e2d2023-7d79-4612-bfe1-37e26cd5639f
# ╠═a8129575-b19a-4c0c-b378-c328a6a8fec7
# ╟─29bacdbd-451f-4f71-91d6-54ae785d1202
# ╟─a78a0dea-de26-44c3-8522-fdde4e48bd0d
# ╟─0dd40cbb-74e3-49a2-bd24-793efa4fb6ea
# ╠═8c8e3b0e-9281-421d-ad65-9a875ec3416b
