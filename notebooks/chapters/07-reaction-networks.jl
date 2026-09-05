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

# ╔═╡ 0770a691-39bf-4bc2-8a3f-c42328ede616
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ b6b44eb3-2ca6-4474-844c-bb5b9e68597d
begin
    using CairoMakie
    using Catalyst
    using JumpProcesses
    using MomentClosure
    using OrdinaryDiffEq
    using PlutoUI
    using StableRNGs
    using Statistics
end

# ╔═╡ 33eda5de-6217-42c2-9100-752c6b22608b
md"""
# Reaction networks and moment closure

The chapter fixed the population at ``N = 1000``, where a second-order normal
closure tracks a jump-process ensemble closely.  The interesting question is
where that stops being true.

Moment closure is an *uncontrolled* approximation: there is no step size to
refine.  Its accuracy depends on how close the distribution of the state is to
the shape the closure assumes, and that depends on population size.  Shrink the
population and the distribution of the infectious class becomes skewed, extinction
becomes common, and a symmetric closure has no way to represent either.

Move the sliders to find where the closure stops following the ensemble.
"""

# ╔═╡ dd184794-f18d-4192-8db3-3ad4df8f06f1
@bind population PlutoUI.Slider(50:50:1000; default = 1000, show_value = true)

# ╔═╡ 6625c919-1bda-46dc-bd5c-65e0d48e1365
@bind trajectories PlutoUI.Slider(50:50:400; default = 200, show_value = true)

# ╔═╡ cd789643-33b2-4437-b882-0e873b078367
@bind horizon PlutoUI.Slider(20:5:60; default = 40, show_value = true)

# ╔═╡ d7aea2c6-c457-4e0b-8a88-4e842b02065b
begin
    sir_network = @reaction_network begin
        κ, S + I --> 2I
        γ, I --> R
    end
    κ_value = 0.05 * 10.0 / population
    parameter_map = [:κ => κ_value, :γ => 0.25]
    tspan = (0.0, Float64(horizon))
    times = collect(0.0:0.5:Float64(horizon))
    S₀ = round(Int, 0.99 * population)
    I₀ = max(round(Int, 0.01 * population), 1)
end

# ╔═╡ 867eae9b-e073-4831-88b0-50291ea770eb
begin
    raw_equations = generate_raw_moment_eqs(sir_network, 2)
    closed_equations = moment_closure(raw_equations, "normal")
    closure_solution = solve(
        ODEProblem(closed_equations, [Float64(S₀), Float64(I₀), 0.0], tspan, parameter_map),
        Tsit5(); saveat = times)
    closure_names = string.(unknowns(closed_equations.odes))
    closure_mean = closure_solution[findfirst(==("μ₀₁₀(t)"), closure_names), :]
    closure_var = closure_solution[findfirst(==("μ₀₂₀(t)"), closure_names), :] .- closure_mean .^ 2
    closure_sd = sqrt.(max.(closure_var, 0.0))
end

# ╔═╡ d1daefa9-3350-450c-8d36-76e10e8033ea
begin
    function simulate_jump(seed)
        problem = JumpProblem(
            sir_network, [:S => S₀, :I => I₀, :R => 0], tspan, parameter_map;
            aggregator = Direct(), rng = StableRNG(seed))
        return Array(solve(problem, SSAStepper(); saveat = times)(times))
    end
    ensemble = cat((simulate_jump(9000 + k) for k in 1:trajectories)...; dims = 3)
    ensemble_mean = dropdims(mean(ensemble; dims = 3), dims = 3)[2, :]
    ensemble_sd = dropdims(std(ensemble; dims = 3), dims = 3)[2, :]
    extinct_fraction = mean(ensemble[2, end, :] .== 0)
    largest_gap = maximum(abs.(closure_mean .- ensemble_mean))
end

# ╔═╡ 681d82a9-99be-4aa7-80f5-4c91c7815dd8
let
    figure = Figure(size = (760, 380))
    axis = Axis(figure[1, 1], xlabel = "Time", ylabel = "Infectious",
                title = "N = $(population), $(trajectories) trajectories")
    band!(axis, times, ensemble_mean .- ensemble_sd, ensemble_mean .+ ensemble_sd;
          color = (:grey, 0.35))
    lines!(axis, times, ensemble_mean; color = :grey25, linewidth = 2,
           label = "Ensemble mean ± sd")
    lines!(axis, times, closure_mean; color = :firebrick, linewidth = 2,
           label = "Normal closure")
    lines!(axis, times, closure_mean .+ closure_sd; color = :firebrick,
           linewidth = 1, linestyle = :dot)
    lines!(axis, times, closure_mean .- closure_sd; color = :firebrick,
           linewidth = 1, linestyle = :dot)
    axislegend(axis; position = :rt, framevisible = false)
    figure
end

# ╔═╡ 30e1bb5e-d5fe-4ce4-b28c-d24f6741cc6e
(
    course_unit = "ch07-reaction-networks",
    status = "complete",
    controls = (; population, trajectories, horizon),
    largest_absolute_gap = round(largest_gap; digits = 3),
    gap_relative_to_population = round(largest_gap / population; digits = 4),
    fraction_extinct_by_end = round(extinct_fraction; digits = 3),
)

# ╔═╡ Cell order:
# ╠═0770a691-39bf-4bc2-8a3f-c42328ede616
# ╠═b6b44eb3-2ca6-4474-844c-bb5b9e68597d
# ╟─33eda5de-6217-42c2-9100-752c6b22608b
# ╠═dd184794-f18d-4192-8db3-3ad4df8f06f1
# ╠═6625c919-1bda-46dc-bd5c-65e0d48e1365
# ╠═cd789643-33b2-4437-b882-0e873b078367
# ╠═d7aea2c6-c457-4e0b-8a88-4e842b02065b
# ╠═867eae9b-e073-4831-88b0-50291ea770eb
# ╠═d1daefa9-3350-450c-8d36-76e10e8033ea
# ╠═681d82a9-99be-4aa7-80f5-4c91c7815dd8
# ╠═30e1bb5e-d5fe-4ce4-b28c-d24f6741cc6e
