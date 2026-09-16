### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 8ca2287b-252c-47d6-827d-a16293f7765e
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 4cdea50e-8056-41ee-a90a-cf29da2d9b36
begin
    using CairoMakie
    using NonlinearSolve
    using OrdinaryDiffEq
    using PlutoUI
end

# ╔═╡ 6b793091-dde4-45bc-8ed3-7d3ac8fcb876
md"""
# Foundations: SIR theory and the Julia workflow

*Course unit `ch01-foundations`.* Turn the parameter dials and watch the whole
epidemic respond: the reproduction number ``R_0``, the trajectory, the threshold
line ``S/N = 1/R_0``, and the analytic final size. The
[book chapter](../../chapters/01-foundations.html) gives the derivations.
"""

# ╔═╡ 9d6d0ac1-4cfd-43a7-9f47-46945b2abf5a
md"""
## The model

With ``\lambda = \beta c I / N``,

```math
\frac{dS}{dt} = -\lambda S, \qquad \frac{dI}{dt} = \lambda S - \gamma I, \qquad \frac{dR}{dt} = \gamma I,
```

```math
R_0 = \frac{\beta c}{\gamma}, \qquad s_\infty = s_0\, e^{-R_0 (1 - s_\infty)} .
```
"""

# ╔═╡ b1e000a2-6224-48e8-b67d-a76532f7fd5c
function sir_ode!(du, u, p, t)
    S, I, R = u
    λ = p.β * p.c * I / p.N
    du[1] = -λ * S
    du[2] = λ * S - p.γ * I
    du[3] = p.γ * I
    return nothing
end

# ╔═╡ f41acd32-4bda-40b1-947b-0b6ca6d23873
md"""
## Explore the parameters

Baseline: ``\beta=0.05``, ``c=10``, ``\gamma=0.25`` (so ``R_0=2``). Drop ``R_0``
below 1 and the epidemic fails to take off.
"""

# ╔═╡ f92c10f8-c80a-4869-aba6-0e6ffba9c78a
@bind β_ui PlutoUI.Slider(0.01:0.005:0.10; default = 0.05, show_value = true)

# ╔═╡ b1d9dda9-3aa1-41b6-96c1-61e66935c92e
@bind c_ui PlutoUI.Slider(2.0:1.0:16.0; default = 10.0, show_value = true)

# ╔═╡ d729cd77-b14b-4b9c-9ff7-c5c75b1b55e0
@bind γ_ui PlutoUI.Slider(0.1:0.05:0.6; default = 0.25, show_value = true)

# ╔═╡ 6b88430c-6cc7-4f82-80d4-67abfe108148
begin
    u0_ch01 = [990.0, 10.0, 0.0]
    p_ch01 = (β = β_ui, c = c_ui, γ = γ_ui, N = 1000.0)
    sol_ch01 = solve(ODEProblem(sir_ode!, u0_ch01, (0.0, 120.0), p_ch01), Tsit5(); saveat = 0.5)
    R0_ch01 = p_ch01.β * p_ch01.c / p_ch01.γ
    s0_ch01 = u0_ch01[1] / p_ch01.N
    s_inf_ch01 = solve(
        NonlinearProblem((s, R0) -> s - s0_ch01 * exp(-R0 * (1 - s)), 0.2, R0_ch01),
        NewtonRaphson()).u
    (R0 = round(R0_ch01; digits = 3),
     s_infinity = round(s_inf_ch01; digits = 3),
     attack_rate = round(1 - s_inf_ch01; digits = 3))
end

# ╔═╡ 1055bec3-709e-4b28-851f-fe39caa2637d
let
    t = sol_ch01.t
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "individuals",
              title = "SIR dynamics, R₀ = $(round(R0_ch01; digits = 2))")
    lines!(ax, t, sol_ch01[1, :]; color = :steelblue, linewidth = 3, label = "S")
    lines!(ax, t, sol_ch01[2, :]; color = :firebrick, linewidth = 3, label = "I")
    lines!(ax, t, sol_ch01[3, :]; color = :seagreen,  linewidth = 3, label = "R")
    if R0_ch01 > 1
        hlines!(ax, [1000.0 / R0_ch01]; color = :black, linestyle = :dash, label = "S = N/R₀")
    end
    axislegend(ax; position = :rc)
    fig
end

# ╔═╡ f2ec3f55-e159-4ded-9b32-f7772add4053
md"""
## Validation

Population is conserved, and when the epidemic grows at first
(``R_0 S(0)/N > 1``) the infectious peak sits on the threshold ``S/N = 1/R_0``.
When it does not, ``I`` falls from the start and its largest value is ``I(0)``,
so the check below is skipped.
"""

# ╔═╡ be6f1bf3-767f-471d-aa43-210e173e33ac
let
    cons = all(u -> isapprox(sum(u), 1000.0; atol = 1e-6), sol_ch01.u)
    I = sol_ch01[2, :]
    peak_ok = if R0_ch01 > 1.05
        isapprox(sol_ch01[1, argmax(I)] / 1000.0, 1 / R0_ch01; atol = 0.02)
    else
        true
    end
    (conservation_ok = cons, threshold_ok = peak_ok)
end

# ╔═╡ cb1e9984-0c0f-47d3-ad48-f7e4d6e5e1f4
md"""
## Takeaways

- ``R_0`` and the initial state set the final size; the recovery rate sets the timing.
- The epidemic grows only while ``S/N > 1/R_0``, so a peak at a positive time sits
  exactly at equality.
- `NonlinearSolve` roots the final-size relation, which cannot be rearranged to
  give the surviving susceptible fraction on its own.
- The deterministic ODE is the large-population limit, not the exact average of a
  finite epidemic; later chapters add stochasticity.
"""

# ╔═╡ 51e8bca4-dc6e-49af-a0f1-3aaabc4199b1
(course_unit = "ch01-foundations", status = "complete")

# ╔═╡ Cell order:
# ╠═8ca2287b-252c-47d6-827d-a16293f7765e
# ╠═4cdea50e-8056-41ee-a90a-cf29da2d9b36
# ╟─6b793091-dde4-45bc-8ed3-7d3ac8fcb876
# ╟─9d6d0ac1-4cfd-43a7-9f47-46945b2abf5a
# ╠═b1e000a2-6224-48e8-b67d-a76532f7fd5c
# ╟─f41acd32-4bda-40b1-947b-0b6ca6d23873
# ╠═f92c10f8-c80a-4869-aba6-0e6ffba9c78a
# ╠═b1d9dda9-3aa1-41b6-96c1-61e66935c92e
# ╠═d729cd77-b14b-4b9c-9ff7-c5c75b1b55e0
# ╠═6b88430c-6cc7-4f82-80d4-67abfe108148
# ╠═1055bec3-709e-4b28-851f-fe39caa2637d
# ╟─f2ec3f55-e159-4ded-9b32-f7772add4053
# ╠═be6f1bf3-767f-471d-aa43-210e173e33ac
# ╟─cb1e9984-0c0f-47d3-ad48-f7e4d6e5e1f4
# ╠═51e8bca4-dc6e-49af-a0f1-3aaabc4199b1
