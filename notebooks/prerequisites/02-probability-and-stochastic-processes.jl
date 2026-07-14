### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 0d39fe3a-5ff5-4462-8f95-c4cabd4678b6
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ 6113098c-3003-402f-8ffd-17a711bddbb4
begin
    using CairoMakie
    using Distributions
    using PlutoUI
    using StableRNGs
end

# ╔═╡ 283ae481-21e4-45da-bf65-1eb72e6b4de0
md"""
# Probability and stochastic processes

*Course unit `p02-probability-and-stochastic-processes`.* This notebook lets you
drive the binomial SIR Markov chain: change the ensemble size, the contact rate
and the step width, and watch the realisations, their mean, and the extinction
probability respond. See the
[book chapter](../../prerequisites/02-probability-and-stochastic-processes.html)
for the derivations.
"""

# ╔═╡ 0aefeda9-ec9c-4b49-a928-231d5583fc12
md"""
## The chain

Per-step probabilities from the SIR rates,

```math
p_I = 1 - e^{-\beta c I_t \delta t / N}, \qquad p_R = 1 - e^{-\gamma \delta t},
```

then a Markov update by binomial draws,

```math
Y_t \sim \mathrm{Binomial}(S_t, p_I), \quad Z_t \sim \mathrm{Binomial}(I_t, p_R),
```

with ``S_{t+1}=S_t-Y_t``, ``I_{t+1}=I_t+Y_t-Z_t``, ``R_{t+1}=R_t+Z_t``. Only the
current counts enter — the Markov property — and ``N`` is conserved.
"""

# ╔═╡ 9fb0dd13-b4ab-4039-b567-9540d8231db7
rate_to_proportion(r, δt) = -expm1(-r * δt)

# ╔═╡ cedac9f1-e04c-4e1b-8628-85c899af2b99
function chain_step(state, p, δt, rng)
    S, I, R = state
    β, c, γ, N = p
    Y = rand(rng, Binomial(S, rate_to_proportion(β * c * I / N, δt)))
    Z = rand(rng, Binomial(I, rate_to_proportion(γ, δt)))
    return (S - Y, I + Y - Z, R + Z)
end

# ╔═╡ 5923a043-6ed2-4b64-b582-fe065dddeaf3
function simulate(u0, p, δt, nsteps, rng)
    traj = Matrix{Int}(undef, 3, nsteps + 1)
    traj[:, 1] .= u0
    s = Tuple(u0)
    for k in 1:nsteps
        s = chain_step(s, p, δt, rng)
        traj[:, k + 1] .= s
    end
    return traj
end

# ╔═╡ 0a3e0067-928b-42bc-bec1-f25234e91fea
md"""
## Run an ensemble

`nsims` realisations, contact rate `c`, step width `δt`; the horizon is fixed at
40 days. Seeds `1:nsims` keep every run reproducible.
"""

# ╔═╡ 47a8f9de-5695-42a3-b207-e37a466263b4
@bind nsims_ui Slider(50:50:400; default = 200, show_value = true)

# ╔═╡ 7e5d72b2-5c34-488c-ab6d-8f3426859919
@bind c_ui Slider(4.0:1.0:14.0; default = 10.0, show_value = true)

# ╔═╡ 8045f3d4-d63e-4eb0-a9fa-51d93f48f3a2
@bind δt_ui Slider([1.0, 0.5, 0.25]; default = 1.0, show_value = true)

# ╔═╡ 1444c183-bf9c-460d-abfb-db9e69df5bb6
begin
    p_ens = (0.05, c_ui, 0.25, 1000)
    nsteps_ens = round(Int, 40 / δt_ui)
    ensemble = [simulate([990, 10, 0], p_ens, δt_ui, nsteps_ens, StableRNG(s))
                for s in 1:nsims_ui]
    I_paths = reduce(hcat, (t[2, :] for t in ensemble))
    mean_I = vec(sum(I_paths; dims = 2)) ./ nsims_ui
    (R0 = 0.05 * c_ui / 0.25, ensemble_size = nsims_ui, steps = nsteps_ens)
end

# ╔═╡ 162a7518-370e-4f57-a92a-3600e32820ac
let
    ts = range(0, 40; length = nsteps_ens + 1)
    fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "infectious I(t)",
              title = "Binomial SIR chain: $(nsims_ui) realisations, R₀ = $(round(0.05 * c_ui / 0.25; digits = 2))")
    for j in 1:nsims_ui
        lines!(ax, ts, I_paths[:, j]; color = (:grey, 0.12), linewidth = 1)
    end
    lines!(ax, ts, mean_I; color = :firebrick, linewidth = 3, label = "ensemble mean")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ a8142103-5016-46a4-b762-c445c9683343
md"""
## Validation and extinction

A single binomial draw's sample mean approaches ``np``; every realisation
conserves ``N``; and we read off the fraction of epidemics that fade out
(``I=0``) before day 40 — a quantity the deterministic model cannot produce.
"""

# ╔═╡ 352abd99-d250-4145-ac72-da4d5cd1e8da
let
    rng = StableRNG(7)
    samples = [rand(rng, Binomial(990, 0.02)) for _ in 1:100_000]
    mean_ok = isapprox(sum(samples) / length(samples), 990 * 0.02; rtol = 0.02)
    cons_ok = all(t -> all(k -> sum(t[:, k]) == 1000, 1:nsteps_ens + 1), ensemble)
    extinct = count(t -> any(==(0), t[2, 2:end]), ensemble) / nsims_ui
    (binomial_mean_ok = mean_ok, conservation_ok = cons_ok,
     extinction_probability = round(extinct; digits = 3))
end

# ╔═╡ edc31565-4196-46f4-b72d-9a0303cc6eef
md"""
## Takeaways

- A Binomial is a sum of independent Bernoulli trials; its variance ``np(1-p)``
  peaks near ``p=0.5``.
- Demographic stochasticity lets small epidemics fade out by chance.
- The ensemble mean tracks the deterministic solution studied in *Foundations*.
- Finer steps reduce the within-step double-event artefact.
"""

# ╔═╡ ca537bf4-2518-4e57-9034-8bfdf48ecc47
(course_unit = "p02-probability-and-stochastic-processes", status = "complete")
# ╔═╡ Cell order:
# ╠═0d39fe3a-5ff5-4462-8f95-c4cabd4678b6
# ╠═6113098c-3003-402f-8ffd-17a711bddbb4
# ╟─283ae481-21e4-45da-bf65-1eb72e6b4de0
# ╟─0aefeda9-ec9c-4b49-a928-231d5583fc12
# ╠═9fb0dd13-b4ab-4039-b567-9540d8231db7
# ╠═cedac9f1-e04c-4e1b-8628-85c899af2b99
# ╠═5923a043-6ed2-4b64-b582-fe065dddeaf3
# ╟─0a3e0067-928b-42bc-bec1-f25234e91fea
# ╠═47a8f9de-5695-42a3-b207-e37a466263b4
# ╠═7e5d72b2-5c34-488c-ab6d-8f3426859919
# ╠═8045f3d4-d63e-4eb0-a9fa-51d93f48f3a2
# ╠═1444c183-bf9c-460d-abfb-db9e69df5bb6
# ╠═162a7518-370e-4f57-a92a-3600e32820ac
# ╟─a8142103-5016-46a4-b762-c445c9683343
# ╠═352abd99-d250-4145-ac72-da4d5cd1e8da
# ╟─edc31565-4196-46f4-b72d-9a0303cc6eef
# ╠═ca537bf4-2518-4e57-9034-8bfdf48ecc47
