### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 515e9f35-fe32-4a27-9c47-31342a750cf6
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ ae91bcaf-6f7b-4ea6-ba39-f64a0c5195c4
begin
    using CairoMakie
    using DiffEqCallbacks
    using OrdinaryDiffEq
    using PlutoUI
end

# ╔═╡ 1681e83a-417d-47dc-b7dd-f3ac609aa3f6
md"""
# Ordinary differential equations

*Course unit `ch03-ordinary-differential-equations`.* Design a lockdown
interactively: choose its strength and window and watch the infectious curve
flatten while a `ContinuousCallback` reports the peak. See the
[book chapter](../../chapters/03-ordinary-differential-equations.html).
"""

# ╔═╡ 38e13e78-2172-4614-8398-b948617cefcc
md"""
## Model and intervention

Base SIR with a time-dependent, piecewise-constant ``\beta(t)``:

```math
\lambda(t) = \frac{\beta(t)\, c\, I}{N}, \qquad
\beta(t) = \begin{cases} \beta_\text{lock}, & t_1 \le t < t_2, \\ \beta, & \text{otherwise.} \end{cases}
```

The peak solves ``dI/dt = \lambda(t) S - \gamma I = 0``.
"""

# ╔═╡ 36abb1f0-ab30-4b9c-9aee-7c9bf7199149
function sir_ode!(du, u, p, t)
    S, I, R = u
    β, c, γ, N = p
    λ = β * c * I / N
    du[1] = -λ * S
    du[2] = λ * S - γ * I
    du[3] = γ * I
    return nothing
end

# ╔═╡ 21994dfc-c43d-442d-a867-670d45bbc517
md"""
## Configure the lockdown

`β_lock` is the reduced transmission probability during the window `[t₁, t₂]`.
The baseline (no intervention) curve is shown for comparison.
"""

# ╔═╡ 3bce08c3-5d3e-421c-bba1-9f148a66769b
@bind βlock_ui Slider(0.010:0.005:0.050; default = 0.025, show_value = true)

# ╔═╡ 194bb316-5374-473a-9ac9-9a9a4e153851
@bind t1_ui Slider(0.0:2.0:30.0; default = 10.0, show_value = true)

# ╔═╡ 803c1a74-252c-4216-a048-1cbcfaa58d0b
@bind window_ui Slider(5.0:5.0:40.0; default = 20.0, show_value = true)

# ╔═╡ 7a104907-389f-4367-bf1e-fa812c66d358
begin
    u0_ch03 = [990.0, 10.0, 0.0]
    p_ch03 = (0.05, 10.0, 0.25, 1000.0)
    tspan_ch03 = (0.0, 80.0)
    t2_ui = t1_ui + window_ui

    sol_base = solve(ODEProblem(sir_ode!, u0_ch03, tspan_ch03, collect(p_ch03)), Tsit5();
                     reltol = 1e-8, abstol = 1e-8)

    lockdown!(integ) = (integ.p[1] = integ.t < t2_ui ? βlock_ui : 0.05)
    cb = PresetTimeCallback([t1_ui, t2_ui], lockdown!)
    sol_lock = solve(ODEProblem(sir_ode!, u0_ch03, tspan_ch03, [0.05, 10.0, 0.25, 1000.0]),
                     Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-8)

    peak_t = Float64[]
    dIdt(u, p) = p[1] * p[2] * u[2] / p[4] * u[1] - p[3] * u[2]
    peak_cb = ContinuousCallback((u, t, i) -> dIdt(u, i.p), i -> push!(peak_t, i.t))
    solve(ODEProblem(sir_ode!, u0_ch03, tspan_ch03, collect(p_ch03)), Tsit5();
          callback = peak_cb, reltol = 1e-8, abstol = 1e-8)

    (window = (t1_ui, t2_ui),
     baseline_peak = round(maximum(sol_base[2, :]); digits = 1),
     lockdown_peak = round(maximum(sol_lock[2, :]); digits = 1),
     detected_peak_time = round.(peak_t; digits = 2))
end

# ╔═╡ 91c476ca-5c92-4ae8-8a1d-df2a7147a2b7
let
    tt = 0.0:0.25:80.0
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "infectious I(t)",
              title = "Lockdown: β → $(βlock_ui) over [$(t1_ui), $(t2_ui)]")
    vspan!(ax, t1_ui, t2_ui; color = (:grey, 0.15))
    lines!(ax, tt, [sol_base(t)[2] for t in tt]; color = :firebrick, linewidth = 3, label = "no intervention")
    lines!(ax, tt, [sol_lock(t)[2] for t in tt]; color = :steelblue, linewidth = 3, linestyle = :dash, label = "lockdown")
    vlines!(ax, peak_t; color = :black, linestyle = :dot, label = "detected peak")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 3781b693-99a9-4ada-ac1e-8e28e0437c4a
md"""
## Validation

The adaptive solver conserves ``N``, its dense interpolant is self-consistent,
and the callback-detected peak matches a brute-force scan of the baseline.
"""

# ╔═╡ 458be579-3533-447c-968c-32bd76c17b3a
let
    cons = all(u -> isapprox(sum(u), 1000.0; atol = 1e-6), sol_base.u)
    interp = isapprox(sol_base(sol_base.t[10]), sol_base.u[10]; atol = 1e-8)
    dense_t = 0.0:0.01:80.0
    brute = dense_t[argmax([sol_base(t)[2] for t in dense_t])]
    peak_ok = isempty(peak_t) ? false : isapprox(peak_t[1], brute; atol = 0.05)
    (conservation_ok = cons, interpolation_ok = interp, peak_matches = peak_ok)
end

# ╔═╡ c1a9ffb9-b396-4e82-91d1-04a2a8f8b409
md"""
## Takeaways

- `PresetTimeCallback` mutates a parameter mid-solve; `ContinuousCallback`
  root-finds an event exactly.
- Dense output means `sol(t)` is a high-order interpolant at *any* time.
- Earlier and stronger lockdowns flatten the peak more but can shift it later.
- The same `ODEProblem` scales to stiff solvers, sensitivity, and fitting.
"""

# ╔═╡ c184c8c2-28e9-4d56-b09d-e12aa2f21f14
(course_unit = "ch03-ordinary-differential-equations", status = "complete")
# ╔═╡ Cell order:
# ╠═515e9f35-fe32-4a27-9c47-31342a750cf6
# ╠═ae91bcaf-6f7b-4ea6-ba39-f64a0c5195c4
# ╟─1681e83a-417d-47dc-b7dd-f3ac609aa3f6
# ╟─38e13e78-2172-4614-8398-b948617cefcc
# ╠═36abb1f0-ab30-4b9c-9aee-7c9bf7199149
# ╟─21994dfc-c43d-442d-a867-670d45bbc517
# ╠═3bce08c3-5d3e-421c-bba1-9f148a66769b
# ╠═194bb316-5374-473a-9ac9-9a9a4e153851
# ╠═803c1a74-252c-4216-a048-1cbcfaa58d0b
# ╠═7a104907-389f-4367-bf1e-fa812c66d358
# ╠═91c476ca-5c92-4ae8-8a1d-df2a7147a2b7
# ╟─3781b693-99a9-4ada-ac1e-8e28e0437c4a
# ╠═458be579-3533-447c-968c-32bd76c17b3a
# ╟─c1a9ffb9-b396-4e82-91d1-04a2a8f8b409
# ╠═c184c8c2-28e9-4d56-b09d-e12aa2f21f14
