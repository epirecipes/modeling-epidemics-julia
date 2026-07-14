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

# ╔═╡ 9e2cf155-5d8e-4c3e-b78d-bc7bb2d37ce4
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 7f08e5d0-48b0-4a18-bfbb-9fcd7cf1c30a
begin
    using CairoMakie
    using EpiModelingCourse
    using OrdinaryDiffEq
    using PlutoUI
    using Random
    using Statistics
    using Surrogates
end

# ╔═╡ 89e3c0ea-7437-4b82-9b12-8f014b9f2c69
md"""
# Surrogate and emulator models

This notebook emulates the finite-horizon recovered fraction of the baseline
SIR model. Change the training size and inspect the held-out error: a good
training fit is not evidence of good interpolation.

**Classification:** deterministic, continuous-time, continuous-state,
population-level.
"""

# ╔═╡ 74ba7e3a-7f90-4726-91a1-a910f72eb82f
@bind n_train Slider(6:2:24; default = 12, show_value = true)

# ╔═╡ 4362acc2-cf0b-480e-93c5-019857b3db56
@bind selected_beta Slider(0.02:0.001:0.08; default = 0.05, show_value = true)

# ╔═╡ 5af2d6de-6776-4c37-9398-cd4ad40e5b3b
begin
    const A07_LB = 0.02
    const A07_UB = 0.08
    const A07_TSPAN = (0.0, 40.0)

    function recovered_fraction(β)
        p = SIRParameters(beta = β, contact_rate = 10.0, gamma = 0.25, N = 1000.0)
        problem = ODEProblem(sir_ode!, [990.0, 10.0, 0.0], A07_TSPAN, p)
        solution = solve(problem, Tsit5(); saveat = [last(A07_TSPAN)],
                         abstol = 1e-9, reltol = 1e-9)
        return solution[3, end] / p.N
    end
end

# ╔═╡ 98f36776-2dc3-4fd4-af7d-0eb35fcc6e79
begin
    Random.seed!(20260711)
    β_train = sample(n_train, A07_LB, A07_UB, LatinHypercubeSample())
    y_train = recovered_fraction.(β_train)
    emulator = RadialBasis(
        β_train,
        y_train,
        A07_LB,
        A07_UB;
        rad = cubicRadial(),
    )
    β_grid = collect(range(A07_LB, A07_UB; length = 64))
    y_true = recovered_fraction.(β_grid)
    y_pred = emulator.(β_grid)
    errors = y_pred .- y_true
end

# ╔═╡ f46f2baa-e5e0-45a3-921c-44cf341cfa82
begin
    fig = Figure(size = (900, 380))
    ax1 = Axis(fig[1, 1], xlabel = "transmission probability β",
        ylabel = "recovered fraction at t = 40")
    lines!(ax1, β_grid, y_true; color = :black, linewidth = 3, label = "SIR")
    lines!(ax1, β_grid, y_pred; color = :darkorange, label = "surrogate")
    scatter!(ax1, β_train, y_train; color = :steelblue, label = "training")
    axislegend(ax1; position = :lt)
    ax2 = Axis(fig[1, 2], xlabel = "β", ylabel = "|surrogate − SIR|")
    scatter!(ax2, β_grid, abs.(errors); color = :firebrick)
    fig
end

# ╔═╡ 7b99d7de-fb55-40d2-b8a5-ea7a51d07346
begin
    selected_error = emulator(selected_beta) - recovered_fraction(selected_beta)
    validation = (
        course_unit = "a07-surrogate-models",
        status = "complete",
        training_points = n_train,
        selected_beta = selected_beta,
        selected_prediction = emulator(selected_beta),
        selected_error = selected_error,
        rmse = sqrt(mean(errors .^ 2)),
        max_abs_error = maximum(abs.(errors)),
    )
    validation
end

# ╔═╡ 0f4aa3b0-64f0-47f0-87e1-56cb91a74f9e
md"""
The radial-basis emulator is an interpolation device on
`[$(A07_LB), $(A07_UB)]`. Increase the training size until the held-out
diagnostic is adequate for the intended use, then re-check any decision with
the mechanistic simulator. `Surrogates.update!` can add a point after a local
error check, but it does not replace a held-out validation design.
"""

# ╔═╡ Cell order:
# ╠═9e2cf155-5d8e-4c3e-b78d-bc7bb2d37ce4
# ╠═7f08e5d0-48b0-4a18-bfbb-9fcd7cf1c30a
# ╟─89e3c0ea-7437-4b82-9b12-8f014b9f2c69
# ╠═74ba7e3a-7f90-4726-91a1-a910f72eb82f
# ╠═4362acc2-cf0b-480e-93c5-019857b3db56
# ╠═5af2d6de-6776-4c37-9398-cd4ad40e5b3b
# ╠═98f36776-2dc3-4fd4-af7d-0eb35fcc6e79
# ╠═f46f2baa-e5e0-45a3-921c-44cf341cfa82
# ╠═7b99d7de-fb55-40d2-b8a5-ea7a51d07346
# ╟─0f4aa3b0-64f0-47f0-87e1-56cb91a74f9e
