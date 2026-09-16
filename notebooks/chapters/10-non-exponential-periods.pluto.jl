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

# ╔═╡ 888131b7-4f62-4e02-86a4-d184bb87452c
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 136fa74d-cf4f-43e6-bd44-6082fa5bae9c
begin
    using CairoMakie
    using DelayDiffEq
    using OrdinaryDiffEq
    using PlutoUI
end

# ╔═╡ a9d7368a-5a1e-4ec0-9c7a-864fbc64f28e
md"""
# Shaping the infectious period

The chapter showed that holding ``R_0`` and the mean infectious period fixed, and
changing only the *shape* of the period distribution, moves the epidemic peak by
more than a factor of two.

The slider sets the number of stages ``n`` in the method of stages. The infectious
period is Erlang with shape ``n`` and rate ``n\gamma``: the mean stays at
``1/\gamma`` for every ``n``, while the coefficient of variation is
``1/\sqrt{n}``.

``n = 1`` is the ordinary SIR model. Increase it and watch the curve climb towards
the fixed-delay solution, which is the ``n \to \infty`` limit.
"""

# ╔═╡ 9fd5330a-94c0-4925-970d-e97580365b4f
@bind stages PlutoUI.Slider(1:1:40; default = 1, show_value = true)

# ╔═╡ 1dc28053-af21-44d5-b266-d0c0ad3078f9
@bind mean_period PlutoUI.Slider(2.0:0.5:8.0; default = 4.0, show_value = true)

# ╔═╡ 591157e3-529a-48a0-92b4-8e15dd5f9747
begin
    βc = 0.5; N = 1000.0; tspan = (0.0, 60.0); times = collect(0.0:0.5:60.0)
    γ = 1/mean_period

    function staged!(du, u, p, t)
        n = p.n
        S = u[1]; I = @view u[2:n+1]
        infection = βc*S*sum(I)/N
        du[1] = -infection
        du[2] = infection - n*γ*I[1]
        @inbounds for j in 2:n
            du[j+1] = n*γ*(I[j-1] - I[j])
        end
        du[end] = n*γ*I[n]
    end

    function run_stages(n)
        u₀ = zeros(n + 2); u₀[1] = 990.0; u₀[2] = 10.0
        sol = solve(ODEProblem(staged!, u₀, tspan, (; n)), Tsit5();
                    saveat = times, abstol = 1e-10, reltol = 1e-10)
        [sum(x[2:n+1]) for x in sol.u]
    end

    function delayed!(du, u, h, p, t)
        now = βc*u[1]*u[2]/N
        past = h(p, t - mean_period)
        old = t < mean_period ? 0.0 : βc*past[1]*past[2]/N
        du[1] = -now; du[2] = now - old; du[3] = old
    end
    hist(p, t) = [990.0, 10.0, 0.0]

    current = run_stages(stages)
    exponential = run_stages(1)
    delayed = solve(DDEProblem(delayed!, [990.0, 10.0, 0.0], hist, tspan, nothing;
                               constant_lags = [mean_period]),
                    MethodOfSteps(Tsit5()); saveat = times, abstol = 1e-9, reltol = 1e-9)
    nothing
end

# ╔═╡ da1736c1-45cd-4320-ae90-11edb6589b05
let
    figure = Figure(size = (760, 380))
    axis = Axis(figure[1, 1], xlabel = "Time", ylabel = "Infectious",
                title = "n = $(stages),  CV = $(round(1/sqrt(stages); digits = 3)),  mean period = $(mean_period)")
    lines!(axis, times, exponential; color = :grey50, linewidth = 2, linestyle = :dot,
           label = "Exponential (n = 1)")
    lines!(axis, times, delayed[2, :]; color = :firebrick, linewidth = 2, linestyle = :dash,
           label = "Fixed delay (n to infinity)")
    lines!(axis, times, current; color = :steelblue, linewidth = 2.5,
           label = "Erlang n = $(stages)")
    axislegend(axis; position = :rt, framevisible = false)
    figure
end

# ╔═╡ 3af90166-ca42-4869-a048-96dce77ad384
(
    course_unit = "ch10-non-exponential-periods",
    status = "complete",
    controls = (; stages, mean_period),
    cv_of_period = round(1/sqrt(stages); digits = 3),
    peak_erlang = round(maximum(current); digits = 2),
    peak_exponential = round(maximum(exponential); digits = 2),
    peak_fixed_delay = round(maximum(delayed[2, :]); digits = 2),
    fraction_of_gap_closed =
        round((maximum(current) - maximum(exponential)) /
              (maximum(delayed[2, :]) - maximum(exponential)); digits = 3),
)

# ╔═╡ 25d8e2cc-37a7-4c6d-be47-e94bb2034377
md"""
The figure holds ``R_0`` fixed: raising ``n`` never changes how transmissible the
infection is, only how the infectious period is distributed.
`fraction_of_gap_closed` reports how far the Erlang model has travelled from the
exponential case towards the fixed-delay limit. Note how slowly it approaches 1,
which is the ``1/\sqrt{n}`` convergence in action.
"""

# ╔═╡ Cell order:
# ╠═888131b7-4f62-4e02-86a4-d184bb87452c
# ╠═136fa74d-cf4f-43e6-bd44-6082fa5bae9c
# ╟─a9d7368a-5a1e-4ec0-9c7a-864fbc64f28e
# ╠═9fd5330a-94c0-4925-970d-e97580365b4f
# ╠═1dc28053-af21-44d5-b266-d0c0ad3078f9
# ╠═591157e3-529a-48a0-92b4-8e15dd5f9747
# ╠═da1736c1-45cd-4320-ae90-11edb6589b05
# ╠═3af90166-ca42-4869-a048-96dce77ad384
# ╟─25d8e2cc-37a7-4c6d-be47-e94bb2034377
