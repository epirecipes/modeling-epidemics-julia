### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c7e994ee-e85a-4cc8-b97c-7975f60443ee
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ f8ade818-176d-42b0-8f42-22d4eff62e7c
begin
    using DelayDiffEq
    using OrdinaryDiffEq
    using CairoMakie
    using PlutoUI
end

# ╔═╡ 37099f8e-34af-49b5-86fa-19ca64d55581
md"""
# Delay differential equations

Companion notebook for `a01-delay-differential-equations`
(see the [book chapter](../../appendices/01-delay-differential-equations.qmd)
for the full derivation).

Move the sliders below to change the fixed recovery delay $\tau$ and the
number of initially infectious individuals $I(0)$, and watch how the
**method of steps** solution and its comparison with the exponential-recovery
ODE respond. Everything below re-solves reactively — nothing is
pre-computed or pasted in.
"""

# ╔═╡ a48ce48b-5d69-44ac-b766-09dbe1c0e6b4
md"""
τ (fixed recovery delay): $(@bind τ PlutoUI.Slider(1.0:0.5:10.0; default = 4.0, show_value = true))

I(0) (initial infections, seeded at t = 0): $(@bind I0 PlutoUI.Slider(1:1:50; default = 10, show_value = true))
"""

# ╔═╡ 85eb8c0c-6bde-4f83-a9ad-68873a3c9966
function sir_dde!(du, u, h, p, t)
    S, I, R = u
    β, c, τp = p
    N = S + I + R
    λ = β * c * I / N
    Sd, Id, Rd = h(p, t - τp)
    Nd = Sd + Id + Rd
    λd = β * c * Id / Nd
    du[1] = -λ * S
    du[2] = λ * S - λd * Sd
    du[3] = λd * Sd
    return nothing
end

# ╔═╡ d746a046-b605-4ef3-ac52-9c6e2d82011b
function sir_ode!(du, u, p, t)
    S, I, R = u
    β, c, γrate = p
    N = S + I + R
    λ = β * c * I / N
    du[1] = -λ * S
    du[2] = λ * S - γrate * I
    du[3] = γrate * I
    return nothing
end

# ╔═╡ 9cf9c141-ae6b-490b-861c-004e16363cd2
begin
    sir_history(p, t) = [1000.0, 0.0, 0.0]

    u0 = [1000.0 - I0, Float64(I0), 0.0]
    β_fixed = 0.05
    c_fixed = 10.0
    tspan = (0.0, 40.0)
    p_dde = (β_fixed, c_fixed, τ)
    γ_matched = 1 / τ

    sol_dde = solve(
        DDEProblem(sir_dde!, u0, sir_history, tspan, p_dde; constant_lags = [τ]),
        MethodOfSteps(Tsit5());
        reltol = 1e-8, abstol = 1e-8,
    )
    sol_ode = solve(
        ODEProblem(sir_ode!, u0, tspan, (β_fixed, c_fixed, γ_matched)),
        Tsit5();
        reltol = 1e-8, abstol = 1e-8,
    )
    nothing
end

# ╔═╡ fb581411-3ba8-4057-b92f-ece68a0387f1
begin
    tt = 0.0:0.2:40.0
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "number",
              title = "Fixed-delay DDE (τ = $(τ)) vs exponential ODE (γ = 1/τ)")
    lines!(ax, tt, [sol_dde(t)[1] for t in tt]; color = :steelblue, linewidth = 3, label = "S (DDE)")
    lines!(ax, tt, [sol_dde(t)[2] for t in tt]; color = :firebrick, linewidth = 3, label = "I (DDE)")
    lines!(ax, tt, [sol_dde(t)[3] for t in tt]; color = :seagreen, linewidth = 3, label = "R (DDE)")
    lines!(ax, tt, [sol_ode(t)[2] for t in tt]; color = :firebrick, linewidth = 2,
           linestyle = :dash, label = "I (ODE)")
    vlines!(ax, [τ]; color = :black, linestyle = :dot, label = "τ")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 58bdef23-9cea-4495-9177-feda24cb3119
begin
    peak_dde, idx_dde = findmax(sol_dde[2, :])
    peak_ode, idx_ode = findmax(sol_ode[2, :])
    conserved = all(t -> isapprox(sum(sol_dde(t)), 1000.0; atol = 1e-4), 0.0:0.5:40.0)
    nonnegative = all(t -> all(>=(-1e-8), sol_dde(t)), 0.0:0.5:40.0)
    md"""
    **Current run** (τ = $(τ), I(0) = $(I0), matched γ = $(round(γ_matched; digits = 4))):

    - DDE peak infectious count: **$(round(peak_dde; digits = 1))** at t ≈ $(round(sol_dde.t[idx_dde]; digits = 2))
    - Matched-mean ODE peak: **$(round(peak_ode; digits = 1))** at t ≈ $(round(sol_ode.t[idx_ode]; digits = 2))
    - DDE peak exceeds ODE peak: **$(peak_dde > peak_ode)**
    - Conservation of N holds on a dense grid: **$(conserved)**
    - All DDE states stay nonnegative on the same grid: **$(nonnegative)**

    As τ grows (holding β, c, N fixed), the gap between the fixed-delay and
    exponential-recovery peaks widens — the fixed delay lets the initial
    cohort transmit at full strength for the *entire* interval $[0, τ)$
    instead of thinning out early, exactly as Exercise 1 in the book chapter
    asks you to quantify.
    """
end

# ╔═╡ 03028e91-748d-447c-a5f3-820b2bcc0022
md"""
## About the method of steps

`MethodOfSteps` solves the DDE as a sequence of ODE problems: on $[0, τ]$ the
delayed term is fully determined by the history function, so the segment is
an ordinary ODE; that segment's dense solution then supplies the history for
$[τ, 2τ]$, and so on. Declaring `constant_lags = [τ]` tells the solver where
the resulting derivative discontinuities are, instead of forcing it to
rediscover them through failed error control.
"""

# ╔═╡ c24f139f-195f-4397-9a94-5cf2fe04283e
(
    course_unit = "a01-delay-differential-equations",
    status = "complete",
    τ = τ,
    I0 = I0,
    dde_peak = round(peak_dde; digits = 2),
    ode_peak = round(peak_ode; digits = 2),
)

# ╔═╡ Cell order:
# ╠═c7e994ee-e85a-4cc8-b97c-7975f60443ee
# ╠═f8ade818-176d-42b0-8f42-22d4eff62e7c
# ╟─37099f8e-34af-49b5-86fa-19ca64d55581
# ╟─a48ce48b-5d69-44ac-b766-09dbe1c0e6b4
# ╠═85eb8c0c-6bde-4f83-a9ad-68873a3c9966
# ╠═d746a046-b605-4ef3-ac52-9c6e2d82011b
# ╠═9cf9c141-ae6b-490b-861c-004e16363cd2
# ╟─fb581411-3ba8-4057-b92f-ece68a0387f1
# ╟─58bdef23-9cea-4495-9177-feda24cb3119
# ╟─03028e91-748d-447c-a5f3-820b2bcc0022
# ╠═c24f139f-195f-4397-9a94-5cf2fe04283e
