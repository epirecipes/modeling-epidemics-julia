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

# ╔═╡ 9b88b3f1-dc31-4f18-a80a-56fd2e3122fd
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 4f2bf6c4-4a59-4df8-81f3-4bd680fbeb09
begin
    using CairoMakie
    using InfiniteOpt
    using Ipopt
    using PlutoUI
end

# ╔═╡ e27e88d6-8712-4f61-aaf3-9a1455d19eed
md"""
# Advanced optimal control

The control `u(t)` is a bounded function, not one scalar intervention time.
`InfiniteOpt` states the ODE-constrained problem and transcribes it with
orthogonal collocation before Ipopt solves the resulting nonlinear program.

**Classification:** deterministic, continuous-time, continuous-state,
population-level.
"""

# ╔═╡ 623f5f4a-8d98-4472-8fb7-95a539940bb9
@bind intervention_budget Slider(0.0:0.5:3.0; default = 2.0, show_value = true)

# ╔═╡ 5a6f01f3-6bf3-4292-bd74-7b6afac01442
function solve_control(budget; support_count = 15)
    tf, βc, γ, umax = 20.0, 0.5, 0.25, 0.5
    model = InfiniteModel(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)
    set_optimizer_attribute(model, "max_iter", 400)
    set_optimizer_attribute(model, "tol", 1e-7)

    @infinite_parameter(
        model,
        t in [0.0, tf],
        num_supports = support_count,
        derivative_method = OrthogonalCollocation(2),
    )
    @variable(model, S >= 0, Infinite(t), start = 0.99)
    @variable(model, I >= 0, Infinite(t), start = 0.01)
    @variable(model, C >= 0, Infinite(t), start = 0.0)
    @variable(model, 0 <= u <= umax, Infinite(t), start = 0.0)
    @constraint(model, S(0.0) == 0.99)
    @constraint(model, I(0.0) == 0.01)
    @constraint(model, C(0.0) == 0.0)
    @constraint(model, deriv(S, t) == -(1 - u) * βc * S * I)
    @constraint(model, deriv(I, t) == (1 - u) * βc * S * I - γ * I)
    @constraint(model, deriv(C, t) == (1 - u) * βc * S * I)
    @constraint(model, integral(u, t) <= budget)
    @objective(model, Min, C(tf) + 1e-3 * integral(u^2, t))
    constant_over_collocation(u, t)
    optimize!(model)

    return (
        status = termination_status(model),
        objective = objective_value(model),
        time = supports(t),
        S = value.(S),
        I = value.(I),
        C = value.(C),
        u = value.(u),
    )
end

# ╔═╡ f87f151a-2b18-4fdd-9837-ea9c78fe0ed3
control_solution = solve_control(intervention_budget)

# ╔═╡ 0c98c0a4-8a38-43f5-a31f-2dba06f78fb9
begin
    ts = control_solution.time
    fig = Figure(size = (900, 380))
    ax1 = Axis(fig[1, 1], xlabel = "time", ylabel = "fraction")
    lines!(ax1, ts, control_solution.S; label = "S")
    lines!(ax1, ts, control_solution.I; label = "I")
    lines!(ax1, ts, control_solution.C; label = "C")
    axislegend(ax1; position = :lt)
    ax2 = Axis(fig[1, 2], xlabel = "time", ylabel = "u(t)")
    stairs!(ax2, ts, control_solution.u; color = :darkorange)
    hlines!(ax2, [0.5]; linestyle = :dash)
    fig
end

# ╔═╡ b9c62f13-0b8e-4ef4-8abc-1fcd8a9f5c2a
begin
    identity_error = maximum(abs.(0.99 .- control_solution.S .- control_solution.C))
    bound_violation = maximum(max.(0.0, control_solution.u .- 0.5))
    diagnostic = (
        course_unit = "a08-advanced-optimal-control",
        status = "complete",
        solver_status = control_solution.status,
        objective = control_solution.objective,
        budget = intervention_budget,
        incidence_identity_error = identity_error,
        control_bound_violation = bound_violation,
        minimum_state = minimum(vcat(control_solution.S, control_solution.I, control_solution.C)),
    )
    diagnostic
end

# ╔═╡ 1fbb57d7-b34f-41c8-b6e2-cd24fa3ef8e0
md"""
Increasing the budget gives the optimizer more intervention pressure, subject
to the pointwise bound. Collocation and Ipopt provide a local numerical
solution; repeat with more supports to assess transcription sensitivity.

The nonlinear SIR dynamics are deliberately **not** passed to
`SDDP.LinearPolicyGraph`: current SDDP requires affine constraints/objectives.
An SDDP proxy must be a separately documented linear approximation, not a claim
to solve this nonlinear control problem.
"""

# ╔═╡ Cell order:
# ╠═9b88b3f1-dc31-4f18-a80a-56fd2e3122fd
# ╠═4f2bf6c4-4a59-4df8-81f3-4bd680fbeb09
# ╟─e27e88d6-8712-4f61-aaf3-9a1455d19eed
# ╠═623f5f4a-8d98-4472-8fb7-95a539940bb9
# ╠═5a6f01f3-6bf3-4292-bd74-7b6afac01442
# ╠═f87f151a-2b18-4fdd-9837-ea9c78fe0ed3
# ╠═0c98c0a4-8a38-43f5-a31f-2dba06f78fb9
# ╠═b9c62f13-0b8e-4ef4-8abc-1fcd8a9f5c2a
# ╟─1fbb57d7-b34f-41c8-b6e2-cd24fa3ef8e0
