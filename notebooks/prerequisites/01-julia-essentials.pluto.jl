### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 93ff1999-cabb-4389-b948-e224c2282e10
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 8157289c-06f3-4b3d-af0e-1f40e112cef9
begin
    using CairoMakie
    using PlutoUI
    using StableRNGs
end

# ╔═╡ 8a2b1d07-6d4b-4ce1-a20a-a38e7945e1b7
md"""
# Julia essentials for modelers

*Course unit `p01-julia-essentials`.* This companion notebook makes the course
building blocks interactive: the in-place SIR kernel, the rate-to-proportion
conversion, and reproducible randomness. Move the sliders to build intuition,
then read the [book chapter](../../prerequisites/01-julia-essentials.html) for
the full narrative.
"""

# ╔═╡ df29fbed-5108-40a1-8156-dcecc6912141
md"""
## The model in one line

With force of infection ``\lambda = \beta c I / N`` the SIR rates are

```math
\frac{dS}{dt} = -\lambda S, \qquad \frac{dI}{dt} = \lambda S - \gamma I, \qquad \frac{dR}{dt} = \gamma I .
```

The three right-hand sides sum to zero, so ``S+I+R=N`` is conserved. For a
discrete step we convert a rate to a probability with
``p(r,\delta t)=1-e^{-r\,\delta t}``.
"""

# ╔═╡ 1496a1de-cd77-428f-b283-2542157c5265
function sir_ode!(du, u, p, t)
    S, I, R = u                  # state: positional
    λ = p.β * p.c * I / p.N      # parameters: by field
    du[1] = -λ * S              # dS/dt
    du[2] = λ * S - p.γ * I     # dI/dt
    du[3] = p.γ * I             # dR/dt
    return nothing
end

# ╔═╡ d765bb92-a948-4680-b0a3-8463f0216c58
rate_to_proportion(r, δt) = -expm1(-r * δt)   # numerically stable 1 - exp(-r δt)

# ╔═╡ e6e241f0-cfbf-451b-a320-195464e7f2e0
md"""
## Explore ``R_0`` and the force of infection

Vary the transmission probability ``\beta``, the contact rate ``c`` and the
initial number infected. The reactive output recomputes the reproduction number
``R_0=\beta c/\gamma`` (with ``\gamma=0.25``) and the instantaneous derivatives.
"""

# ╔═╡ 07b520ab-e1b9-4a6f-baf5-5972fa4eb066
@bind β_ui PlutoUI.Slider(0.01:0.005:0.10; default = 0.05, show_value = true)

# ╔═╡ 69adda8a-fffa-43e8-98a7-01a1280dc334
@bind c_ui PlutoUI.Slider(1.0:1.0:20.0; default = 10.0, show_value = true)

# ╔═╡ 61a48731-15bd-4825-ab04-28cbde813542
@bind I0_ui PlutoUI.Slider(1:1:400; default = 10, show_value = true)

# ╔═╡ d99acdeb-7b39-47a9-8593-3ba7c32b62a6
begin
    p_ui = (β = β_ui, c = c_ui, γ = 0.25, N = 1000.0)
    u_ui = [1000.0 - I0_ui, float(I0_ui), 0.0]
    du_ui = zeros(3)
    sir_ode!(du_ui, u_ui, p_ui, 0.0)
    (R0 = p_ui.β * p_ui.c / p_ui.γ,
     force_of_infection = p_ui.β * p_ui.c * I0_ui / p_ui.N,
     dSdt = du_ui[1], dIdt = du_ui[2], dRdt = du_ui[3])
end

# ╔═╡ 7df8f556-06fc-4a9d-b19a-d962c27f2851
md"""
## Rate-to-proportion, reactively

Slide the rate ``r`` and watch how quickly the per-step event probability
saturates towards 1.
"""

# ╔═╡ c0199edb-1cca-42e4-a777-6cc9d1f1bd14
@bind r_ui PlutoUI.Slider(0.1:0.1:3.0; default = 1.0, show_value = true)

# ╔═╡ f5a07aad-0bcf-43c1-8d04-e93363f955bd
let
    δts = 0.0:0.02:2.0
    fig = Figure(size = (620, 360))
    ax = Axis(fig[1, 1]; xlabel = "time step δt", ylabel = "1 − exp(−r δt)",
              title = "Rate-to-proportion at r = $(round(r_ui; digits = 2))")
    lines!(ax, δts, rate_to_proportion.(r_ui, δts); linewidth = 3, color = :firebrick)
    ylims!(ax, 0, 1)
    fig
end

# ╔═╡ ed5977e3-e705-4f5d-807f-1d2c20997cf6
md"""
## Reproducible randomness

A *local* `StableRNG` gives bit-identical draws across machines and Julia
versions. Change the seed to get a different — but reproducible — stream.
"""

# ╔═╡ fd63fc57-50f4-4259-b7d7-43958eea9d07
@bind seed_ui PlutoUI.Slider(1:100; default = 24, show_value = true)

# ╔═╡ 62033df1-95ac-44ba-9071-671f33f9500e
rand(StableRNG(seed_ui), 5)

# ╔═╡ 11f39bc7-94cc-4fcb-8679-9a6f6c00c72a
md"""
## Validation

At the baseline state the derivatives must sum to zero (conservation), and the
rate conversion must match its Taylor limit ``p \approx r\,\delta t`` as
``\delta t \to 0``.
"""

# ╔═╡ db4704f3-5963-4861-b56a-94fb4acb5a68
let
    du = zeros(3)
    sir_ode!(du, [990.0, 10.0, 0.0], (β = 0.05, c = 10.0, γ = 0.25, N = 1000.0), 0.0)
    (conserved = isapprox(sum(du), 0.0; atol = 1e-12),
     taylor = isapprox(rate_to_proportion(0.25, 1e-6), 0.25e-6; rtol = 1e-3))
end

# ╔═╡ 6af3306e-5158-482a-b965-8d0fddd47aa1
md"""
## Takeaways

- In-place kernels `f!(du, u, p, t)` write derivatives into `du` and return
  `nothing`. SciML also accepts an out-of-place `f(u, p, t)`; this course uses
  the in-place form everywhere.
- State is read by position (`S, I, R = u`), parameters by field (`p.β`) — the
  convention every later chapter follows.
- Broadcasting (`.`) evaluates a scalar function on a whole grid without a loop.
- `expm1` avoids cancellation error in ``1 - e^{-x}``.
- Seed a local `StableRNG` to fix the random stream across machines and Julia
  versions; what you compute from it can still shift in the last digits.
"""

# ╔═╡ a5c3fead-9b7a-4e6f-aa8e-669439ced39f
(course_unit = "p01-julia-essentials", status = "complete")

# ╔═╡ Cell order:
# ╠═93ff1999-cabb-4389-b948-e224c2282e10
# ╠═8157289c-06f3-4b3d-af0e-1f40e112cef9
# ╟─8a2b1d07-6d4b-4ce1-a20a-a38e7945e1b7
# ╟─df29fbed-5108-40a1-8156-dcecc6912141
# ╠═1496a1de-cd77-428f-b283-2542157c5265
# ╠═d765bb92-a948-4680-b0a3-8463f0216c58
# ╟─e6e241f0-cfbf-451b-a320-195464e7f2e0
# ╠═07b520ab-e1b9-4a6f-baf5-5972fa4eb066
# ╠═69adda8a-fffa-43e8-98a7-01a1280dc334
# ╠═61a48731-15bd-4825-ab04-28cbde813542
# ╠═d99acdeb-7b39-47a9-8593-3ba7c32b62a6
# ╟─7df8f556-06fc-4a9d-b19a-d962c27f2851
# ╠═c0199edb-1cca-42e4-a777-6cc9d1f1bd14
# ╠═f5a07aad-0bcf-43c1-8d04-e93363f955bd
# ╟─ed5977e3-e705-4f5d-807f-1d2c20997cf6
# ╠═fd63fc57-50f4-4259-b7d7-43958eea9d07
# ╠═62033df1-95ac-44ba-9071-671f33f9500e
# ╟─11f39bc7-94cc-4fcb-8679-9a6f6c00c72a
# ╠═db4704f3-5963-4861-b56a-94fb4acb5a68
# ╟─6af3306e-5158-482a-b965-8d0fddd47aa1
# ╠═a5c3fead-9b7a-4e6f-aa8e-669439ced39f
