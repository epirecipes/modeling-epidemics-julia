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

# ╔═╡ dd240ddc-772e-4d48-988d-9dadcb9e76d4
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ bc602d1c-5800-453c-9675-b3903196e69e
begin
    using EpiModelingCourse
    using Distributions
    using StableRNGs
    using CairoMakie
    using PlutoUI
end

# ╔═╡ e9344bf1-1d6a-45dd-b5dd-20e6bf7ab6c2
md"""
# Observation models and simulated data

Companion notebook for the `ch07-observation-models` unit. The latent SIR
process is fixed; the sliders below control the **observation model** that turns
latent daily incidence into reported case counts. Everything is reproducible: a
`StableRNG` seed makes each realisation deterministic.
"""

# ╔═╡ fe4b3563-7c27-4b3f-9868-7de49e37e4e1
md"""
**Reporting fraction ρ** (share of latent cases reported):
"""

# ╔═╡ ea3ff3f7-17ee-4ca0-aac6-0fa2e49cf6b0
@bind reporting_fraction PlutoUI.Slider(0.1:0.1:1.0; default = 1.0, show_value = true)

# ╔═╡ 5a136e57-895a-42d7-ba11-6fc20b0b7338
md"""**Observation model:**"""

# ╔═╡ 5a6e54c2-b394-4848-88b1-c120f496dfc3
@bind observation_model PlutoUI.Select(["Poisson", "Negative-binomial"]; default = "Poisson")

# ╔═╡ c2cdb3a3-8cc6-4ba9-8657-9082f7beebc9
md"""**Negative-binomial dispersion r** (smaller = more overdispersed):"""

# ╔═╡ 653b11de-f8d1-452d-9a86-44a3e42e2f9a
@bind dispersion PlutoUI.Select([1.0 => "r = 1", 5.0 => "r = 5", 10.0 => "r = 10",
    50.0 => "r = 50 (≈ Poisson)"]; default = 10.0)

# ╔═╡ ba54c943-4e2c-4e71-9f15-a772e1732b7e
md"""**Random seed:**"""

# ╔═╡ f1eef571-92f0-4307-9573-05a9c6fb9026
@bind seed PlutoUI.Slider(1:100; default = 42, show_value = true)

# ╔═╡ 55c0f9e4-bb85-4eb7-a4c5-25e31cf1c316
begin
    latent_df = synthetic_incidence(; tspan = (0.0, 40.0), saveat = 1.0,
        observation_noise = :none)
    obstimes = latent_df.time[2:end]
    latent = latent_df.incidence[2:end]
end;

# ╔═╡ bd16e29e-873e-40ef-8942-bc288b56d5cc
begin
    expected = reporting_fraction .* latent
    rng = StableRNG(seed)
    observed = if observation_model == "Poisson"
        rand.(rng, Poisson.(expected))
    else
        [rand(rng, NegativeBinomial(dispersion, dispersion / (dispersion + μ)))
         for μ in expected]
    end
end;

# ╔═╡ f62aaa37-7bab-4fb8-ac52-1039f915a335
let
    fig = Figure(size = (760, 400))
    ax = Axis(fig[1, 1]; xlabel = "time (days)", ylabel = "new cases per day",
        title = "$(observation_model) reports, ρ = $(reporting_fraction)")
    barplot!(ax, collect(obstimes), observed; color = (:darkorange, 0.55),
        label = "reported")
    lines!(ax, collect(obstimes), latent; color = :firebrick, linewidth = 3,
        label = "latent incidence")
    lines!(ax, collect(obstimes), expected; color = :steelblue, linewidth = 2,
        linestyle = :dash, label = "ρ × latent")
    axislegend(ax; position = :rt)
    fig
end

# ╔═╡ f0d97e4a-b563-4e72-8ad4-dbecadcd0ca4
md"""
The dashed line is the expected report count `ρ × latent`. Notice that the
Poisson scatter hugs its mean, while the negative-binomial scatter widens as `r`
shrinks. Reducing `ρ` shrinks the *apparent* epidemic even though the underlying
process is unchanged — the identifiability problem revisited in Chapter 10.
"""

# ╔═╡ 669eb73c-17ed-4f52-9657-3be8a5a1c8a5
(
    course_unit = "ch07-observation-models",
    status = "complete",
    reporting_fraction = reporting_fraction,
    observation_model = observation_model,
    dispersion = dispersion,
    seed = seed,
    total_latent = round(sum(latent); digits = 1),
    total_reported = sum(observed),
    var_over_mean = round(var(observed) / mean(observed); digits = 2),
)

# ╔═╡ Cell order:
# ╠═dd240ddc-772e-4d48-988d-9dadcb9e76d4
# ╠═bc602d1c-5800-453c-9675-b3903196e69e
# ╟─e9344bf1-1d6a-45dd-b5dd-20e6bf7ab6c2
# ╟─fe4b3563-7c27-4b3f-9868-7de49e37e4e1
# ╠═ea3ff3f7-17ee-4ca0-aac6-0fa2e49cf6b0
# ╟─5a136e57-895a-42d7-ba11-6fc20b0b7338
# ╠═5a6e54c2-b394-4848-88b1-c120f496dfc3
# ╟─c2cdb3a3-8cc6-4ba9-8657-9082f7beebc9
# ╠═653b11de-f8d1-452d-9a86-44a3e42e2f9a
# ╟─ba54c943-4e2c-4e71-9f15-a772e1732b7e
# ╠═f1eef571-92f0-4307-9573-05a9c6fb9026
# ╠═55c0f9e4-bb85-4eb7-a4c5-25e31cf1c316
# ╠═bd16e29e-873e-40ef-8942-bc288b56d5cc
# ╠═f62aaa37-7bab-4fb8-ac52-1039f915a335
# ╟─f0d97e4a-b563-4e72-8ad4-dbecadcd0ca4
# ╠═669eb73c-17ed-4f52-9657-3be8a5a1c8a5
