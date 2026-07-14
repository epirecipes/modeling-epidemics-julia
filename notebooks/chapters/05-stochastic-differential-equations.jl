### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 29190e7d-56ee-4647-b7b1-fbb2f376f62d
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 9db38bed-cd2e-4947-b3b4-8cde0c84a3cb
begin
    using CairoMakie
    using OrdinaryDiffEq
    using PlutoUI
    using StableRNGs
    using StochasticDiffEq
end

# ╔═╡ 6d4136f6-94ce-4d85-8d95-e81e0b3e38d4
md"""
# Stochastic differential equations

The chemical-Langevin SIR model replaces each reaction count by a drift plus
Gaussian noise.  The two noise columns below correspond to infection and
recovery, so their stoichiometry preserves \(S+I+R\).
"""

# ╔═╡ 9d41926f-e1e8-4c26-b1f0-3911f7c027d6
begin
    @bind step_choice PlutoUI.Slider(1:3; default = 2, show_value = false)
    @bind transmission PlutoUI.Slider(0.02:0.01:0.10; default = 0.05, show_value = true)
    @bind horizon PlutoUI.Slider(10:5:40; default = 30, show_value = true)
    dt_options = (0.01, 0.02, 0.05)
    dt = dt_options[step_choice]
end

# ╔═╡ c4cb108c-38c4-48c1-a5c2-3cfa26bdb7e0
begin
    N = 1000.0
    c = 10.0
    γ = 0.25
    κ = transmission * c / N
    u₀ = [990.0, 10.0, 0.0]
    times = collect(0.0:0.1:Float64(horizon))

    function propensities(u, κ, γ)
        S, I = max(u[1], 0.0), max(u[2], 0.0)
        return κ * S * I, γ * I
    end

    function drift!(du, u, p, t)
        infection, recovery = propensities(u, p[1], p[2])
        du[1] = -infection
        du[2] = infection - recovery
        du[3] = recovery
        return nothing
    end

    function diffusion!(G, u, p, t)
        infection, recovery = propensities(u, p[1], p[2])
        a = sqrt(infection)
        b = sqrt(recovery)
        G[1, 1], G[2, 1], G[3, 1] = -a, a, 0.0
        G[1, 2], G[2, 2], G[3, 2] = 0.0, -b, b
        return nothing
    end
end

# ╔═╡ f96fe6ed-d5ee-471a-9cac-5b7d24fb4aec
begin
    sde_problem = SDEProblem(
        drift!,
        diffusion!,
        u₀,
        (0.0, Float64(horizon)),
        (κ, γ);
        noise_rate_prototype = zeros(3, 2),
    )
    sde_solution = solve(
        sde_problem,
        LambaEM();
        dt,
        adaptive = false,
        saveat = times,
        rng = StableRNG(20250711),
    )
    sde_states = Array(sde_solution)
end

# ╔═╡ 3e4cffab-56fc-4be7-bb54-5fad730f5de9
begin
    ode_solution = solve(
        ODEProblem(drift!, u₀, (0.0, Float64(horizon)), (κ, γ)),
        Tsit5();
        saveat = times,
    )
    ode_states = Array(ode_solution)
    mass_error = maximum(abs.(vec(sum(sde_states; dims = 1)) .- N))
    minimum_state = minimum(sde_states)
    @assert mass_error < 1e-8
    (; dt, mass_error, minimum_state)
end

# ╔═╡ 8c7c89a0-725d-49ae-b738-e8ab5c25d523
begin
    fig = Figure(size = (760, 420))
    ax = Axis(fig[1, 1], xlabel = "time", ylabel = "people",
        title = "Chemical-Langevin SIR")
    for (j, (name, color)) in enumerate(zip(["S", "I", "R"],
                                             [:steelblue, :firebrick, :seagreen]))
        lines!(ax, times, sde_states[j, :]; color, linewidth = 2,
            label = "$name (SDE)")
        lines!(ax, times, ode_states[j, :]; color, linestyle = :dash,
            linewidth = 1.2, label = "$name (ODE)")
    end
    axislegend(ax; position = :rt, nbanks = 2)
    fig
end

# ╔═╡ 30ea9b57-1526-44fc-ad5d-bf7bbd40fe8a
md"""
`mass_error` is a structural check.  `minimum_state` is deliberately reported:
an unconstrained Gaussian diffusion can cross a population boundary, where the
chemical-Langevin approximation is no longer biologically interpretable.
"""

# ╔═╡ 2b7578a8-0c2d-4fd9-a33c-cba54aefbd74
(
    course_unit = "ch05-stochastic-differential-equations",
    status = "complete",
    controls = (; dt, transmission, horizon),
)

# ╔═╡ Cell order:
# ╠═29190e7d-56ee-4647-b7b1-fbb2f376f62d
# ╠═9db38bed-cd2e-4947-b3b4-8cde0c84a3cb
# ╟─6d4136f6-94ce-4d85-8d95-e81e0b3e38d4
# ╠═9d41926f-e1e8-4c26-b1f0-3911f7c027d6
# ╠═c4cb108c-38c4-48c1-a5c2-3cfa26bdb7e0
# ╠═f96fe6ed-d5ee-471a-9cac-5b7d24fb4aec
# ╠═3e4cffab-56fc-4be7-bb54-5fad730f5de9
# ╠═8c7c89a0-725d-49ae-b738-e8ab5c25d523
# ╟─30ea9b57-1526-44fc-ad5d-bf7bbd40fe8a
# ╠═2b7578a8-0c2d-4fd9-a33c-cba54aefbd74
