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

# ╔═╡ 4260b035-3d72-4cd6-b8a4-9d544d3a7373
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 1758b121-b621-4fe6-af8b-aef1949cd17e
begin
    using OrdinaryDiffEq
    using SparseArrays
    using LinearAlgebra
    using StableRNGs
    using CairoMakie
    using PlutoUI
end

# ╔═╡ ca6995de-560c-4584-bd6d-afbdb2b1ba5d
md"""
# Finite-state projection

Companion notebook for the `a05-finite-state-projection` unit. Instead of
sampling trajectories, we solve directly for the **probability distribution**
over ``(S, I)`` states by integrating a truncated chemical master equation
``d\mathbf{p}/dt = \mathbf{Q}\,\mathbf{p}``.

We build the sparse generator ``\mathbf{Q}`` **transparently** — we do *not* use
`FiniteStateProjection.jl`, whose resolved release cannot be loaded without
mutating a third-party module. The sliders let you (i) watch the distribution
evolve in time and (ii) deliberately **truncate** the retained infectious count
and see probability leak out of the box.
"""

# ╔═╡ a7b771dc-c379-42bd-8fb5-a5009fcf962d
md"""
The reachable states form a triangle ``\{(s,i): s+i \le N\}``. The kernel below
enumerates it, indexes each reachable ``(s,i)``, and fills the sparse CTMC
generator from the infection propensity ``\beta s i`` and recovery propensity
``\gamma i``. An optional `Imax` caps the retained infectious count so we can
study truncation error.
"""

# ╔═╡ 110a9e6e-e9df-4e7e-8462-5898dd717344
function sir_fsp_generator(N, β, γ; Smax = N, Imax = N)
    index = Dict{Tuple{Int,Int},Int}()
    states = Tuple{Int,Int}[]
    for s in 0:Smax, i in 0:Imax
        if s + i <= N
            push!(states, (s, i))
            index[(s, i)] = length(states)
        end
    end
    rows = Int[]; cols = Int[]; vals = Float64[]
    add!(r, c, v) = (push!(rows, r); push!(cols, c); push!(vals, v))
    for (col, (s, i)) in enumerate(states)
        a_inf = β * s * i
        a_rec = γ * i
        add!(col, col, -(a_inf + a_rec))
        a_inf > 0 && haskey(index, (s - 1, i + 1)) && add!(index[(s - 1, i + 1)], col, a_inf)
        a_rec > 0 && haskey(index, (s, i - 1)) && add!(index[(s, i - 1)], col, a_rec)
    end
    return sparse(rows, cols, vals, length(states), length(states)), states, index
end

# ╔═╡ 75acedd5-2367-408c-adc5-a06f50522368
begin
    Ncap = 100
    β = 0.005
    γ = 0.25
    tspan = (0.0, 40.0)
    counts = 0:Ncap
end

# ╔═╡ 41eaa5cd-4927-45e9-a588-7a5927d619e6
md"""
### Explore the projected master equation

**Retained infectious count ``I_{\max}``** (smaller ⇒ tighter truncation ⇒ more
leakage):
"""

# ╔═╡ 472b54aa-a63d-4931-a508-d7bbaa6fe135
@bind Imax_cap PlutoUI.Slider(10:5:100; default = 100, show_value = true)

# ╔═╡ e93a4d57-d6aa-4c8b-beda-db1ff8485e81
md"""**Report the marginal distribution at time ``t``:**"""

# ╔═╡ b627fb0a-3f02-42be-aace-fdbb45834a42
@bind t_report PlutoUI.Slider(0.0:1.0:40.0; default = 20.0, show_value = true)

# ╔═╡ cbd0d043-12a5-4931-9ce9-c4cc866d84be
begin
    Q, states, index = sir_fsp_generator(Ncap, β, γ; Imax = Imax_cap)
    p0 = zeros(length(states))
    p0[index[(99, 1)]] = 1.0
    sol_fsp = solve(ODEProblem((dp, p, Q, t) -> mul!(dp, Q, p), p0, tspan, Q),
        Vern7(); abstol = 1e-10, reltol = 1e-8, saveat = 1.0)
    (n_states = length(states), nnz_Q = nnz(Q))
end

# ╔═╡ b865398d-362a-4843-9c1f-85ed0d3495d6
begin
    marginalI(p) = (m = zeros(Ncap + 1);
        for (k, (s, i)) in enumerate(states); m[i + 1] += p[k]; end; m)
    p_t = sol_fsp(t_report)
    PI = marginalI(p_t)
    mean_I = sum(counts .* PI)
    total_prob_t = sum(p_t)
    leakage = 1 - minimum(sum(u) for u in sol_fsp.u)
    (mean_I = round(mean_I; digits = 3),
     total_prob_at_t = round(total_prob_t; sigdigits = 6),
     max_leakage = round(leakage; sigdigits = 3))
end

# ╔═╡ 75f38e45-7f13-4036-987d-1f8899f0c1ab
begin
    total_prob_curve = [sum(u) for u in sol_fsp.u]

    fig = Figure(size = (900, 400))
    ax1 = Axis(fig[1, 1]; xlabel = "infectious count I",
        ylabel = "probability",
        title = "P(I) at t = $(t_report)  (mean ≈ $(round(mean_I; digits = 2)))")
    barplot!(ax1, counts, PI; color = (:firebrick, 0.75))
    vlines!(ax1, [mean_I]; color = :black, linestyle = :dash)
    xlims!(ax1, -1, 45)

    ax2 = Axis(fig[1, 2]; xlabel = "time (days)", ylabel = "total probability Σp",
        title = "Conservation (leak ≈ $(round(leakage; sigdigits = 2)))")
    lines!(ax2, sol_fsp.t, total_prob_curve; color = :steelblue, linewidth = 2)
    hlines!(ax2, [1.0]; color = :gray, linestyle = :dot)
    ylims!(ax2, min(0.9, minimum(total_prob_curve) - 0.02), 1.01)
    fig
end

# ╔═╡ 6b1c2ad5-a23b-4cd2-9efd-8fb68571c826
md"""
**What to notice.** At ``I_{\max} = 100`` (the full triangle) no probability can
escape: the right panel sits flat at 1 and the leakage is at solver-tolerance
level. Drag ``I_{\max}`` down and, once the epidemic pushes ``I`` toward the cap,
probability leaks out — the total drops below 1. That shortfall is the
finite-state-projection **truncation error**, and it is the built-in diagnostic
for whether the projection is trustworthy.
"""

# ╔═╡ 620c3243-4d5d-4c94-be6b-7999b6f44bb7
begin
    @assert isapprox(total_prob_t, sum(p_t); atol = 1e-12)
    @assert all(PI .>= -1e-8) "marginal probabilities must be nonnegative (up to solver tolerance)"
    (
        course_unit = "a05-finite-state-projection",
        status = "complete",
        controls = (; Imax_cap, t_report),
        n_states = length(states),
        max_leakage = round(leakage; sigdigits = 3),
    )
end

# ╔═╡ Cell order:
# ╠═4260b035-3d72-4cd6-b8a4-9d544d3a7373
# ╠═1758b121-b621-4fe6-af8b-aef1949cd17e
# ╟─ca6995de-560c-4584-bd6d-afbdb2b1ba5d
# ╟─a7b771dc-c379-42bd-8fb5-a5009fcf962d
# ╠═110a9e6e-e9df-4e7e-8462-5898dd717344
# ╠═75acedd5-2367-408c-adc5-a06f50522368
# ╟─41eaa5cd-4927-45e9-a588-7a5927d619e6
# ╠═472b54aa-a63d-4931-a508-d7bbaa6fe135
# ╟─e93a4d57-d6aa-4c8b-beda-db1ff8485e81
# ╠═b627fb0a-3f02-42be-aace-fdbb45834a42
# ╠═cbd0d043-12a5-4931-9ce9-c4cc866d84be
# ╠═b865398d-362a-4843-9c1f-85ed0d3495d6
# ╠═75f38e45-7f13-4036-987d-1f8899f0c1ab
# ╟─6b1c2ad5-a23b-4cd2-9efd-8fb68571c826
# ╠═620c3243-4d5d-4c94-be6b-7999b6f44bb7
