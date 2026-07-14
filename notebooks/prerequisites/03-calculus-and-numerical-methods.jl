### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ fa144ee6-6cb3-4a3e-b2bf-af8d2dda1505
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ bf981898-c2fd-40d5-ba7d-20b9befc73f1
begin
    using CairoMakie
    using PlutoUI
end

# ╔═╡ e5a22a06-f02f-459d-baf1-5bd8bf8a51ef
md"""
# Calculus and numerical methods

*Course unit `p03-calculus-and-numerical-methods`.* This notebook lets you feel
the accuracy/cost trade-off of numerical integration. Pick a step size and a
method (forward Euler or midpoint) and compare against a fine reference; then
watch the empirical order of convergence. See the
[book chapter](../../prerequisites/03-calculus-and-numerical-methods.html).
"""

# ╔═╡ 8830ee6f-2b9c-466d-b017-ab3cd8b0e272
md"""
## Forward Euler and midpoint

Dropping the quadratic Taylor remainder gives forward Euler,

```math
u_{n+1} = u_n + \delta t\, f(u_n),
```

with global error ``O(\delta t)``. The midpoint (2nd-order Runge–Kutta) update

```math
u_{n+1} = u_n + \delta t\, f\!\left(u_n + \tfrac{\delta t}{2} f(u_n)\right)
```

has global error ``O(\delta t^2)`` — one order higher.
"""

# ╔═╡ f1de4368-ec10-47e2-a972-b20941a35d68
function sir_ode!(du, u, p, t)
    S, I, R = u
    β, c, γ, N = p
    λ = β * c * I / N
    du[1] = -λ * S
    du[2] = λ * S - γ * I
    du[3] = γ * I
    return nothing
end

# ╔═╡ bce0d344-8b1c-40f9-936a-227617649f31
function integrate(f, u0, p, δt, tmax; method = :euler)
    ts = 0.0:δt:tmax
    us = Matrix{Float64}(undef, length(u0), length(ts))
    us[:, 1] .= u0
    k1 = similar(u0); k2 = similar(u0)
    for n in 1:length(ts) - 1
        u = @view us[:, n]
        f(k1, u, p, ts[n])
        if method === :euler
            us[:, n + 1] .= u .+ δt .* k1
        else
            f(k2, u .+ (δt / 2) .* k1, p, ts[n] + δt / 2)
            us[:, n + 1] .= u .+ δt .* k2
        end
    end
    return collect(ts), us
end

# ╔═╡ fc798f26-bd1d-4e6c-88e4-5d5f3299a3eb
md"""
## Choose a step and a method

The reference is a very fine forward-Euler run, so we measure convergence to the
integrator's own exact solution at the final time ``t=40``.
"""

# ╔═╡ ae063eee-953c-44ad-95d6-27a45bfa6e28
@bind method_ui Select(["euler" => "Forward Euler", "midpoint" => "Midpoint (RK2)"]; default = "euler")

# ╔═╡ 5a3f1dce-8efb-42c8-96ee-fccc8606e557
@bind logstep_ui Slider(3:9; default = 6, show_value = true)

# ╔═╡ d4306274-d0ab-4050-b271-be3737b543e7
begin
    u0_p03 = [990.0, 10.0, 0.0]
    p_p03 = (0.05, 10.0, 0.25, 1000.0)
    tmax_p03 = 40.0
    method_sym = Symbol(method_ui)
    δt_sel = tmax_p03 / 2^logstep_ui
    _, us_ref_p03 = integrate(sir_ode!, u0_p03, p_p03, tmax_p03 / 2^16, tmax_p03)
    u_ref_p03 = us_ref_p03[:, end]
    ts_sel, us_sel = integrate(sir_ode!, u0_p03, p_p03, δt_sel, tmax_p03; method = method_sym)
    (method = method_ui, δt = round(δt_sel; sigdigits = 3),
     final_state = round.(us_sel[:, end]; digits = 3),
     final_error = round(maximum(abs.(us_sel[:, end] .- u_ref_p03)); sigdigits = 3))
end

# ╔═╡ ef4c7bf6-dc1c-4a47-b0f6-72f4b75c38ea
let
    fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "individuals",
              title = "$(method_ui) at δt = $(round(δt_sel; sigdigits = 3))")
    lines!(ax, ts_sel, us_sel[1, :]; color = :steelblue, linewidth = 3, label = "S")
    lines!(ax, ts_sel, us_sel[2, :]; color = :firebrick, linewidth = 3, label = "I")
    lines!(ax, ts_sel, us_sel[3, :]; color = :seagreen,  linewidth = 3, label = "R")
    axislegend(ax; position = :rc)
    fig
end

# ╔═╡ 067ceba6-4675-4869-8eac-3a5825dd6bdb
md"""
## Empirical order of convergence

The slope of ``\log(\text{error})`` against ``\log(\delta t)`` should be ``\approx 1``
for Euler and ``\approx 2`` for midpoint.
"""

# ╔═╡ 93c7e95b-3cc1-465c-91b9-1d83023d9259
let
    δts = [tmax_p03 / 2^k for k in 5:11]
    errs = map(δts) do dt
        _, us = integrate(sir_ode!, u0_p03, p_p03, dt, tmax_p03; method = method_sym)
        maximum(abs.(us[:, end] .- u_ref_p03))
    end
    orders = [log2(errs[k] / errs[k + 1]) for k in 1:length(errs) - 1]
    expected = method_sym === :euler ? 1.0 : 2.0
    mean_order = sum(orders) / length(orders)
    (mean_order = round(mean_order; digits = 3), expected = expected,
     matches = isapprox(mean_order, expected; atol = 0.15))
end

# ╔═╡ 26670a67-0ce5-45c5-a9a3-e6942396442c
md"""
## Takeaways

- Forward Euler is first order: halving the step only halves the error.
- Midpoint is second order — far more accuracy per unit work.
- Too-large steps break stability; adaptive solvers (`Tsit5`) automate the choice.
"""

# ╔═╡ 3e2a2d44-22b6-4636-b210-265106dff41c
(course_unit = "p03-calculus-and-numerical-methods", status = "complete")
# ╔═╡ Cell order:
# ╠═fa144ee6-6cb3-4a3e-b2bf-af8d2dda1505
# ╠═bf981898-c2fd-40d5-ba7d-20b9befc73f1
# ╟─e5a22a06-f02f-459d-baf1-5bd8bf8a51ef
# ╟─8830ee6f-2b9c-466d-b017-ab3cd8b0e272
# ╠═f1de4368-ec10-47e2-a972-b20941a35d68
# ╠═bce0d344-8b1c-40f9-936a-227617649f31
# ╟─fc798f26-bd1d-4e6c-88e4-5d5f3299a3eb
# ╠═ae063eee-953c-44ad-95d6-27a45bfa6e28
# ╠═5a3f1dce-8efb-42c8-96ee-fccc8606e557
# ╠═d4306274-d0ab-4050-b271-be3737b543e7
# ╠═ef4c7bf6-dc1c-4a47-b0f6-72f4b75c38ea
# ╟─067ceba6-4675-4869-8eac-3a5825dd6bdb
# ╠═93c7e95b-3cc1-465c-91b9-1d83023d9259
# ╟─26670a67-0ce5-45c5-a9a3-e6942396442c
# ╠═3e2a2d44-22b6-4636-b210-265106dff41c
