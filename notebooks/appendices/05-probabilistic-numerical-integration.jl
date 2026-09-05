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

# ╔═╡ 85e02015-f31a-4440-a7af-5586bfedebff
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ bfa34744-6f39-4419-b601-e26ff0f68a33
begin
    using EpiModelingCourse
    using OrdinaryDiffEq
    using CairoMakie
    using PlutoUI
end

# ╔═╡ 53effeeb-9b02-4dad-a1cc-eb994ecd0822
md"""
# Probabilistic numerical integration

Companion notebook for the `a05-probabilistic-numerical-integration` unit. A
probabilistic ODE solver returns a *distribution* over trajectories whose spread
reflects the discretisation error. Here we make that error **concrete and
executable**: we measure the true error of a working solver against a
machine-precision reference.

> **Environment note.** The pinned `ProbNumDiffEq.jl` v0.17 does not precompile
> against the resolved `OrdinaryDiffEq` stack (`UndefVarError: standardtag`), so
> its solvers are **not** loaded here. This notebook runs the deterministic
> error-quantification core; the book chapter documents the gated probabilistic
> API. Slide the controls to see how the *measured* error responds to the
> solver's tolerance and output step.
"""

# ╔═╡ a7f2300a-2aec-4c9c-820e-a06b609842e6
md"""
Baseline SIR (``S_0=990,\,I_0=10,\,N=1000,\,c=10,\,\beta=0.05,\,\gamma=0.25``).
The reference is a high-order `Vern9` solve at tolerance ``10^{-12}``; the
working solver is `Tsit5` at an adjustable tolerance.
"""

# ╔═╡ 22b27fe7-6ac7-40d6-bb42-8ff4d8c62c18
begin
    p = BASELINE_PARAMETERS
    s0 = BASELINE_INITIAL_STATE
    u0 = [s0.S, s0.I, s0.R]
    tspan = (0.0, 40.0)
    prob = ODEProblem(sir_ode!, u0, tspan, p)
end

# ╔═╡ a8a4d078-86cd-485f-97f0-bceeed25c51a
md"""
### Explore the numerical error

**Working-solver tolerance** (log₁₀; smaller ⇒ more accurate, more work):
"""

# ╔═╡ cbebfb14-4529-4016-88ad-d1fee0abc19e
@bind log_tol PlutoUI.Slider(-8:1:-1; default = -1, show_value = true)

# ╔═╡ 780b2999-59ff-4538-9088-80d46ff90d57
md"""**Output step ``\delta t``** (dense-output spacing at which error is sampled):"""

# ╔═╡ 58e0d932-6c83-45fc-bab8-805e66f67ee2
@bind δt PlutoUI.Slider([0.1, 0.25, 0.5, 1.0]; default = 0.5, show_value = true)

# ╔═╡ 4be34aa3-c1cc-47d5-879b-fbb8233a5c1d
begin
    tol = 10.0^log_tol
    ref = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = δt)
    working = solve(prob, Tsit5(); abstol = tol, reltol = tol, saveat = δt)
    refm = reduce(hcat, ref.u)
    err = abs.(reduce(hcat, working.u) .- refm)
    max_error = (S = maximum(err[1, :]), I = maximum(err[2, :]), R = maximum(err[3, :]))
    pop_drift = maximum(abs.([sum(u) for u in ref.u] .- p.N))
    (tolerance = tol, δt = δt, max_error = map(x -> round(x; sigdigits = 3), max_error),
     reference_pop_drift = round(pop_drift; sigdigits = 3))
end

# ╔═╡ 82c58d32-9384-48a2-8a8e-899c91f12256
begin
    fig = Figure(size = (900, 400))
    ax1 = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "absolute error",
        yscale = log10, title = "Error vs. reference (tol = $(tol))")
    for (i, (lab, col)) in enumerate(zip(("S", "I", "R"),
            (:steelblue, :firebrick, :seagreen)))
        lines!(ax1, ref.t, max.(err[i, :], 1e-12); color = col, linewidth = 2, label = lab)
    end
    axislegend(ax1; position = :rb)

    # Refinement sweep at the current δt: error should fall as tolerance tightens.
    tols = [1e-1, 1e-3, 1e-5, 1e-7]
    Ierr = Float64[]
    for tt in tols
        s = solve(prob, Tsit5(); abstol = tt, reltol = tt, saveat = δt)
        push!(Ierr, maximum(abs.(reduce(hcat, s.u)[2, :] .- refm[2, :])))
    end
    ax2 = Axis(fig[1, 2]; xlabel = "solver tolerance", ylabel = "max |I error|",
        xscale = log10, yscale = log10, title = "Calibration: error shrinks with tolerance")
    scatterlines!(ax2, tols, max.(Ierr, 1e-12); color = :purple, linewidth = 2)
    fig
end

# ╔═╡ d03f62af-0110-44a0-898a-c02faaa085af
md"""
**What to notice.** The error is largest around the **epidemic peak**, where the
vector field changes fastest — exactly where a probabilistic solver would widen
its posterior. Tightening the tolerance (left slider) pushes the whole error
curve down, and the right panel shows that decrease monotonically: this is the
deterministic version of a well-*calibrated* uncertainty estimate. Changing the
output step ``\delta t`` barely moves the error, because `Tsit5` is *adaptive* —
``\delta t`` only sets where the dense output is sampled, not the internal steps.
"""

# ╔═╡ 5e8aaae5-47f9-410b-b26f-80f5e9c48cd9
begin
    @assert pop_drift < 1e-6 "reference must conserve the population to near machine precision"
    @assert issorted(Ierr; rev = true) "error must decrease as tolerance tightens"
    (
        course_unit = "a05-probabilistic-numerical-integration",
        status = "complete",
        controls = (; log_tol, δt),
        max_I_error = round(max_error.I; sigdigits = 3),
    )
end

# ╔═╡ Cell order:
# ╠═85e02015-f31a-4440-a7af-5586bfedebff
# ╠═bfa34744-6f39-4419-b601-e26ff0f68a33
# ╟─53effeeb-9b02-4dad-a1cc-eb994ecd0822
# ╟─a7f2300a-2aec-4c9c-820e-a06b609842e6
# ╠═22b27fe7-6ac7-40d6-bb42-8ff4d8c62c18
# ╟─a8a4d078-86cd-485f-97f0-bceeed25c51a
# ╠═cbebfb14-4529-4016-88ad-d1fee0abc19e
# ╟─780b2999-59ff-4538-9088-80d46ff90d57
# ╠═58e0d932-6c83-45fc-bab8-805e66f67ee2
# ╠═4be34aa3-c1cc-47d5-879b-fbb8233a5c1d
# ╠═82c58d32-9384-48a2-8a8e-899c91f12256
# ╟─d03f62af-0110-44a0-898a-c02faaa085af
# ╠═5e8aaae5-47f9-410b-b26f-80f5e9c48cd9
