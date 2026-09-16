### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ b5346d49-231d-4996-9740-1bf819a85255
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 0129f82b-6b57-4128-acf6-5002dcd69f64
begin
    using CairoMakie
    using Distributions
    using PlutoUI
    using StableRNGs
end

# ╔═╡ 39a9278b-2347-489c-bef6-5c3132de6e80
md"""
# Discrete-time models

*Course unit `ch02-discrete-time-models`.* Compare the deterministic function
map with the stochastic Markov chain as you change the step width, ensemble
size, and contact rate. The map is the mean-field shadow of the chain. See the
[book chapter](../../chapters/02-discrete-time-models.html).
"""

# ╔═╡ 50b7348e-4495-4b6e-8bbf-80cf154c002c
md"""
## Two updates

Per-step probabilities ``p_I = 1 - e^{-\beta c I \delta t / N}`` and
``p_R = 1 - e^{-\gamma \delta t}``. The map moves expectations,

```math
S_{t+1} = S_t - p_I S_t, \quad I_{t+1} = I_t + p_I S_t - p_R I_t, \quad R_{t+1} = R_t + p_R I_t,
```

while the chain draws ``Y_t \sim \mathrm{Bin}(S_t, p_I)`` and
``Z_t \sim \mathrm{Bin}(I_t, p_R)``.
"""

# ╔═╡ ec1b9acc-6a7a-4666-b167-54e8e27263dc
rate_to_proportion(r, δt) = -expm1(-r * δt)

# ╔═╡ b34d43be-7024-40f0-9382-73fa40b5dbbc
function map_step(state, p, δt)
    S, I, R = state
    inf = rate_to_proportion(p.β * p.c * I / p.N, δt) * S
    rec = rate_to_proportion(p.γ, δt) * I
    return (S - inf, I + inf - rec, R + rec)
end

# ╔═╡ 31a272f9-1e6b-4d63-bd9c-7b9ff5b0117c
function chain_step(state, p, δt, rng)
    S, I, R = state
    Y = rand(rng, Binomial(round(Int, S), rate_to_proportion(p.β * p.c * I / p.N, δt)))
    Z = rand(rng, Binomial(round(Int, I), rate_to_proportion(p.γ, δt)))
    return (S - Y, I + Y - Z, R + Z)
end

# ╔═╡ 401ec509-3607-4ec7-ae3d-7b67a1bcaa3f
function iterate_model(step, u0, p, δt, nsteps; extra = ())
    traj = Matrix{Float64}(undef, 3, nsteps + 1)
    traj[:, 1] .= u0
    s = Tuple(float.(u0))
    for k in 1:nsteps
        s = step(s, p, δt, extra...)
        traj[:, k + 1] .= s
    end
    return traj
end

# ╔═╡ bde32edc-6bb3-4cc4-be67-92775bb8de83
md"""
## Compare map and chain

Horizon fixed at 40 days; seeds `1:nsims` keep the ensemble reproducible.
"""

# ╔═╡ 180bd361-54a1-4324-8edf-08fac7a85f99
@bind δt_ui PlutoUI.Slider([1.0, 0.5, 0.25]; default = 1.0, show_value = true)

# ╔═╡ bc04ca12-a6e3-427f-b3c0-a2bd73de4c68
@bind nsims_ui PlutoUI.Slider(50:50:400; default = 200, show_value = true)

# ╔═╡ c0de5689-1c3f-4daf-ac94-f88b188da985
@bind c_ui PlutoUI.Slider(4.0:1.0:14.0; default = 10.0, show_value = true)

# ╔═╡ 50496868-d2a8-4641-8216-417dd62fb286
begin
    u0_ch02 = [990, 10, 0]
    p_ch02 = (β = 0.05, c = c_ui, γ = 0.25, N = 1000.0)
    nsteps_ch02 = round(Int, 40 / δt_ui)
    map_traj = iterate_model((s, p, dt) -> map_step(s, p, dt), u0_ch02, p_ch02, δt_ui, nsteps_ch02)
    chain_ens = [iterate_model(chain_step, u0_ch02, p_ch02, δt_ui, nsteps_ch02;
                               extra = (StableRNG(s),)) for s in 1:nsims_ui]
    chain_I = reduce(hcat, (t[2, :] for t in chain_ens))
    chain_mean_I = vec(sum(chain_I; dims = 2)) ./ nsims_ui
    (R0 = p_ch02.β * p_ch02.c / p_ch02.γ, steps = nsteps_ch02,
     map_peak = round(maximum(map_traj[2, :]); digits = 1),
     chain_mean_peak = round(maximum(chain_mean_I); digits = 1))
end

# ╔═╡ 62086158-c21c-436d-a51b-432f3c950fef
let
    ts = range(0, 40; length = nsteps_ch02 + 1)
    fig = Figure(size = (720, 420))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "infectious I(t)",
              title = "Map vs stochastic chain (δt = $(δt_ui))")
    for j in 1:nsims_ui
        lines!(ax, ts, chain_I[:, j]; color = (:grey, 0.1), linewidth = 1)
    end
    lines!(ax, ts, chain_mean_I; color = :firebrick, linewidth = 3, label = "chain mean")
    lines!(ax, ts, map_traj[2, :]; color = :black, linestyle = :dash, linewidth = 3, label = "function map")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ 03e8961e-df3e-4d34-8cdb-513a681cc24e
md"""
## Validation and extinction

Both formulations conserve ``N``. The extinction probability — the fraction of
chains reaching ``I=0`` before day 40 — is invisible to the deterministic map.
"""

# ╔═╡ 08497c2a-28cb-40b5-9b23-7a15ae7964c3
let
    map_ok = all(k -> isapprox(sum(map_traj[:, k]), 1000.0; atol = 1e-8), 1:nsteps_ch02 + 1)
    chain_ok = all(t -> all(k -> sum(t[:, k]) == 1000, 1:nsteps_ch02 + 1), chain_ens)
    extinct = count(t -> any(==(0), t[2, 2:end]), chain_ens) / nsims_ui
    (map_conserves = map_ok, chain_conserves = chain_ok,
     extinction_probability = round(extinct; digits = 3))
end

# ╔═╡ b91dba49-ea00-4baa-bec8-22f29dbe6e9c
md"""
## Takeaways

- The map is the chain's one-step conditional expectation, so the red mean and
  dashed map run close together. They do not coincide: the update is nonlinear,
  so iterating the expectation is not the same as averaging the chain.
- The chain shows chance extinction and run-to-run variability whose size
  *relative to the population* shrinks as the population grows.
- As ``\delta t \to 0`` the map converges to the ODE of *Foundations*, while the
  chain converges to the continuous-time Markov chain. These are two different
  limits.
"""

# ╔═╡ 5480ae15-f262-45a1-a701-824ab7f2b365
(course_unit = "ch02-discrete-time-models", status = "complete")

# ╔═╡ Cell order:
# ╠═b5346d49-231d-4996-9740-1bf819a85255
# ╠═0129f82b-6b57-4128-acf6-5002dcd69f64
# ╟─39a9278b-2347-489c-bef6-5c3132de6e80
# ╟─50b7348e-4495-4b6e-8bbf-80cf154c002c
# ╠═ec1b9acc-6a7a-4666-b167-54e8e27263dc
# ╠═b34d43be-7024-40f0-9382-73fa40b5dbbc
# ╠═31a272f9-1e6b-4d63-bd9c-7b9ff5b0117c
# ╠═401ec509-3607-4ec7-ae3d-7b67a1bcaa3f
# ╟─bde32edc-6bb3-4cc4-be67-92775bb8de83
# ╠═180bd361-54a1-4324-8edf-08fac7a85f99
# ╠═bc04ca12-a6e3-427f-b3c0-a2bd73de4c68
# ╠═c0de5689-1c3f-4daf-ac94-f88b188da985
# ╠═50496868-d2a8-4641-8216-417dd62fb286
# ╠═62086158-c21c-436d-a51b-432f3c950fef
# ╟─03e8961e-df3e-4d34-8cdb-513a681cc24e
# ╠═08497c2a-28cb-40b5-9b23-7a15ae7964c3
# ╟─b91dba49-ea00-4baa-bec8-22f29dbe6e9c
# ╠═5480ae15-f262-45a1-a701-824ab7f2b365
