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

# ╔═╡ ee675fee-7940-410c-8375-e87a2f8fc33a
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ d836fa61-6d05-45ee-9e7b-adef29b704d6
begin
    using CairoMakie
    using OrdinaryDiffEq
    using PlutoUI
end

# ╔═╡ f4371c5d-e954-449b-bde3-1737c6d31780
md"""
# Age of infection: shaping the infectious period

The chapter showed that a constant recovery hazard makes the infectious period
exponential, and that the age-of-infection PDE reproduces the ODE exactly in that
case.

Here you can change the *shape* of the hazard. The sliders set a delay before
recovery becomes possible and the hazard applied afterwards. A long delay with a
sharp hazard approximates a **fixed** infectious period; zero delay recovers the
exponential case.

Watch what happens to the epidemic peak while ``R_0`` is held roughly constant:
concentrating recovery around a single time shortens the generation interval and
makes the epidemic faster and sharper.
"""

# ╔═╡ c887f3ac-d689-4ee1-b449-7de9b0540da8
@bind delay PlutoUI.Slider(0.0:0.5:8.0; default = 0.0, show_value = true)

# ╔═╡ 2a96a048-9a8f-4023-b442-afc50bae3ef1
@bind hazard PlutoUI.Slider(0.25:0.25:4.0; default = 0.25, show_value = true)

# ╔═╡ deba5f71-18f8-4e22-a4ab-a418bcf3ef6b
@bind bins PlutoUI.Slider(100:100:800; default = 400, show_value = true)

# ╔═╡ 189ec02c-892c-4cf9-97f0-085b51914927
function age_of_infection!(du, u, p, t)
    (; Δa, n, βvec, γvec, N) = p
    S = u[1]
    I = @view u[2:n+1]
    λ = sum(βvec .* I) * Δa / N
    du[1] = -λ * S
    du[2] = -(I[1] - λ*S)/Δa - γvec[1]*I[1]
    @inbounds for j in 2:n
        du[j+1] = -(I[j] - I[j-1])/Δa - γvec[j]*I[j]
    end
    du[end] = sum(γvec .* I)*Δa + I[n]
    return nothing
end

# ╔═╡ 3a2ebdf1-61f3-4817-a559-4f9438a664ba
begin
    N = 1000.0; amax = 40.0; times = collect(0.0:0.5:40.0)
    function run_model(γ_profile; n = bins)
        Δa = amax/n
        ages = range(Δa/2, amax - Δa/2, length = n)
        p = (Δa = Δa, n = n, βvec = fill(0.5, n), γvec = γ_profile.(ages), N = N)
        u₀ = zeros(n + 2); u₀[1] = 990.0; u₀[2] = 10.0/Δa
        sol = solve(ODEProblem(age_of_infection!, u₀, (0.0, 40.0), p), Tsit5();
                    saveat = times, abstol = 1e-9, reltol = 1e-9)
        prevalence = [sum(x[2:n+1])*Δa for x in sol.u]
        total = [x[1] + sum(x[2:n+1])*Δa + x[end] for x in sol.u]
        return (; prevalence, total)
    end
    shaped = run_model(a -> a < delay ? 0.0 : hazard)
    exponential = run_model(a -> 0.25)
    drift = maximum(shaped.total) - minimum(shaped.total)
end

# ╔═╡ 78f67b4e-a15b-4b12-abcb-74d8e0db279e
let
    figure = Figure(size = (760, 380))
    axis = Axis(figure[1, 1], xlabel = "Time", ylabel = "Infectious",
                title = "delay = $(delay), hazard = $(hazard), bins = $(bins)")
    lines!(axis, times, exponential.prevalence; color = :steelblue, linewidth = 2,
           linestyle = :dash, label = "Exponential (γ = 0.25)")
    lines!(axis, times, shaped.prevalence; color = :firebrick, linewidth = 2,
           label = "Shaped hazard")
    axislegend(axis; position = :rt, framevisible = false)
    figure
end

# ╔═╡ b18c828a-a125-4988-bfb9-e078051ae37c
(
    course_unit = "ch04-partial-differential-equations",
    status = "complete",
    controls = (; delay, hazard, bins),
    peak_shaped = round(maximum(shaped.prevalence); digits = 3),
    peak_exponential = round(maximum(exponential.prevalence); digits = 3),
    peak_time_shaped = times[argmax(shaped.prevalence)],
    peak_time_exponential = times[argmax(exponential.prevalence)],
    population_drift = round(drift; digits = 9),
)

# ╔═╡ Cell order:
# ╠═ee675fee-7940-410c-8375-e87a2f8fc33a
# ╠═d836fa61-6d05-45ee-9e7b-adef29b704d6
# ╟─f4371c5d-e954-449b-bde3-1737c6d31780
# ╠═c887f3ac-d689-4ee1-b449-7de9b0540da8
# ╠═2a96a048-9a8f-4023-b442-afc50bae3ef1
# ╠═deba5f71-18f8-4e22-a4ab-a418bcf3ef6b
# ╠═189ec02c-892c-4cf9-97f0-085b51914927
# ╠═3a2ebdf1-61f3-4817-a559-4f9438a664ba
# ╠═78f67b4e-a15b-4b12-abcb-74d8e0db279e
# ╠═b18c828a-a125-4988-bfb9-e078051ae37c
