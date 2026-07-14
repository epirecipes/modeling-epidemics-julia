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

# ╔═╡ 3d41b2cd-4f94-4eb0-a1ec-e682982d2b1b
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 8fa6ad0c-8a2e-47c0-bf7d-84d928d635d9
begin
    using CairoMakie
    using Catalyst
    using ModelingToolkit
    using OrdinaryDiffEq
    using PlutoUI
    using ModelingToolkit: t_nounits as t, D_nounits as D
end

# ╔═╡ 0dc671a0-ff7a-4ef7-840a-7f5605e45540
md"""
# Model composition and structured populations

This notebook has two experiments:

1. compose named ModelingToolkit `System` components and compile the
   connections;
2. vary cross-group mixing in a Catalyst reaction network.

The first is deterministic, continuous-time, continuous-state. The second is
population structured and can be interpreted as a deterministic ODE or as a
stochastic reaction model.
"""

# ╔═╡ 84e1caf3-e809-447d-8ca3-651202d3ed89
@bind cross_mixing Slider(0.0:0.1:1.0; default = 0.5, show_value = true)

# ╔═╡ 2d05f568-7257-4e96-a972-d55726bfc535
begin
    @variables S(t) I(t) R(t)
    @parameters β γ

    @named susceptible_component = System(
        [D(S) ~ -β * S * I / (S + I + R)],
        t,
    )
    @named infectious_component = System(
        [D(I) ~ β * S * I / (S + I + R) - γ * I],
        t,
    )
    @named recovered_component = System([D(R) ~ γ * I], t)

    sir_composed = compose(
        System(
            [
                S ~ infectious_component.S,
                I ~ susceptible_component.I,
                R ~ susceptible_component.R,
                infectious_component.S ~ susceptible_component.S,
                susceptible_component.I ~ infectious_component.I,
                susceptible_component.R ~ recovered_component.R,
                infectious_component.R ~ recovered_component.R,
                recovered_component.I ~ infectious_component.I,
            ],
            t,
            [S, I, R],
            [β, γ],
            bindings = [
                susceptible_component.β => β,
                infectious_component.β => β,
                infectious_component.γ => γ,
                recovered_component.γ => γ,
            ],
            name = :sir,
        ),
        susceptible_component,
        infectious_component,
        recovered_component,
    )
    compiled_sir = mtkcompile(sir_composed)
end

# ╔═╡ fcbf5748-7974-4c86-92cb-2d580d03fdfa
begin
    component_problem = ODEProblem(
        compiled_sir,
        Dict(
            susceptible_component.S => 990.0,
            infectious_component.I => 10.0,
            recovered_component.R => 0.0,
            β => 0.05,
            γ => 0.25,
        ),
        (0.0, 10.0),
    )
    component_solution = solve(component_problem, Tsit5(); saveat = 0.1)
end

# ╔═╡ 4c02e53c-6dd4-4dd1-bd48-5e76b99d7346
begin
    two_group_network = @reaction_network begin
        @parameters N
        (β * c1) / N, S1 + I1 --> 2I1
        (β * c1 * m12) / N, S1 + I2 --> I1 + I2
        (β * c2 * m21) / N, S2 + I1 --> I2 + I1
        (β * c2) / N, S2 + I2 --> 2I2
        γ, I1 --> R1
        γ, I2 --> R2
    end

    group_problem = ODEProblem(
        two_group_network,
        [
            :S1 => 495.0, :I1 => 5.0, :R1 => 0.0,
            :S2 => 495.0, :I2 => 5.0, :R2 => 0.0,
        ],
        (0.0, 4.0),
        [
            :N => 1000.0, :β => 0.05, :c1 => 10.0, :c2 => 10.0,
            :m12 => cross_mixing, :m21 => cross_mixing, :γ => 0.25,
        ],
    )
    group_solution = solve(group_problem, Tsit5(); saveat = 0.05)
end

# ╔═╡ 2f921c7a-0c3a-47d5-8b39-05f2b1de873d
begin
    fig = Figure(size = (900, 380))
    ax1 = Axis(fig[1, 1], xlabel = "time", ylabel = "I", title = "composed SIR")
    lines!(ax1, component_solution.t, component_solution[2, :]; color = :firebrick)
    ax2 = Axis(fig[1, 2], xlabel = "time", ylabel = "I₁ + I₂",
        title = "two-group Catalyst ODE")
    lines!(ax2, group_solution.t, group_solution[:I1] .+ group_solution[:I2],
        color = :steelblue)
    fig
end

# ╔═╡ 5cb30acd-5bcc-4db3-8da7-a9cd983e1972
begin
    population = group_solution[:S1] .+ group_solution[:I1] .+
        group_solution[:R1] .+ group_solution[:S2] .+
        group_solution[:I2] .+ group_solution[:R2]
    diagnostics = (
        course_unit = "a09-model-composition",
        status = "complete",
        cross_mixing = cross_mixing,
        compiled_equations = length(equations(compiled_sir)),
        component_population_error = maximum(abs.(
            component_solution[1, :] .+ component_solution[2, :] .+
            component_solution[3, :] .- 1000.0)),
        group_population_error = maximum(abs.(population .- 1000.0)),
    )
    diagnostics
end

# ╔═╡ Cell order:
# ╠═3d41b2cd-4f94-4eb0-a1ec-e682982d2b1b
# ╠═8fa6ad0c-8a2e-47c0-bf7d-84d928d635d9
# ╟─0dc671a0-ff7a-4ef7-840a-7f5605e45540
# ╠═84e1caf3-e809-447d-8ca3-651202d3ed89
# ╠═2d05f568-7257-4e96-a972-d55726bfc535
# ╠═fcbf5748-7974-4c86-92cb-2d580d03fdfa
# ╠═4c02e53c-6dd4-4dd1-bd48-5e76b99d7346
# ╠═2f921c7a-0c3a-47d5-8b39-05f2b1de873d
# ╠═5cb30acd-5bcc-4db3-8da7-a9cd983e1972
