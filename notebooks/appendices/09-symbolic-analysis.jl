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

# ╔═╡ e6e72ddd-0fc4-4798-b9f1-f811a2b8702f
begin
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
end

# ╔═╡ d0e9244d-7bb1-43fc-87a1-f79aee62739f
begin
    using CairoMakie
    using LinearAlgebra
    using OrdinaryDiffEq
    using PlutoUI
    using Symbolics
end

# ╔═╡ eb1d5d3d-b8e6-48cf-96d9-6dd7063dc652
md"""
# Symbolic equilibrium and stability analysis

The appendix derived, for the SIR model with demography, that the endemic
equilibrium exists and is stable exactly when ``R_0 = \beta c/(\gamma+\mu) > 1``,
and that the approach to it can be oscillatory.

Here you can watch the eigenvalues move. The sliders change the parameters; the
symbolic expressions are evaluated at those values, the Routh--Hurwitz verdict is
recomputed, and the eigenvalues are plotted in the complex plane.

Cross ``R_0 = 1`` and watch an eigenvalue cross into the right half-plane.
"""

# ╔═╡ 40d9ff7e-8f42-46e8-a3ab-41e7c07dacb6
@bind transmission PlutoUI.Slider(0.10:0.01:0.90; default = 0.50, show_value = true)

# ╔═╡ 9fdb247f-c197-497c-a46d-8d2151652a81
@bind death_rate PlutoUI.Slider(0.002:0.002:0.060; default = 0.010, show_value = true)

# ╔═╡ f103b25c-f003-44a2-a4ef-3aea3bf0b9aa
@bind recovery PlutoUI.Slider(0.10:0.05:0.60; default = 0.25, show_value = true)

# ╔═╡ fb781900-b2dc-455c-89dc-c5f3863c481b
begin
    @variables S I λ βc γ μ N
    f₂ = [μ*N - βc*S*I/N - μ*S, βc*S*I/N - γ*I - μ*I]
    J₂ = Symbolics.jacobian(f₂, [S, I])
    S_star = N*(γ + μ)/βc
    I_star = μ*N*(βc - γ - μ)/(βc*(γ + μ))
    J₂_endemic = Symbolics.simplify.(
        Symbolics.substitute.(J₂, (Dict(S => S_star, I => I_star),)); expand = true)
    nothing
end

# ╔═╡ 3b7a8970-4c7d-4411-b844-b317ad6c552e
begin
    values = Dict(βc => transmission, γ => recovery, μ => death_rate, N => 1000.0)
    ev(expr) = Float64(Symbolics.value(Symbolics.substitute(expr, values)))
    R₀ = transmission / (recovery + death_rate)
    Jn = [ev(J₂_endemic[i, j]) for i in 1:2, j in 1:2]
    eigenvalues = eigvals(Jn)
    trace_value = Jn[1, 1] + Jn[2, 2]
    det_value = Jn[1, 1]*Jn[2, 2] - Jn[1, 2]*Jn[2, 1]
    routh_hurwitz = trace_value < 0 && det_value > 0
    oscillatory = any(!iszero, imag.(eigenvalues))
    nothing
end

# ╔═╡ 508b3437-be60-4426-bd6f-d0fc179c75c3
let
    figure = Figure(size = (760, 360))
    axis = Axis(figure[1, 1], xlabel = "Re(λ)", ylabel = "Im(λ)",
                title = "R₀ = $(round(R₀; digits = 3)) — " *
                        (R₀ > 1 ? "endemic equilibrium exists" : "not positive; disease-free state is stable"))
    vlines!(axis, [0.0]; color = :grey, linestyle = :dash)
    scatter!(axis, real.(eigenvalues), imag.(eigenvalues);
             color = routh_hurwitz ? :seagreen : :firebrick, markersize = 14)
    figure
end

# ╔═╡ 8e62c4a4-7fb2-43a9-9b0d-bb9f003c3278
(
    course_unit = "a09-symbolic-analysis",
    status = "complete",
    controls = (; transmission, recovery, death_rate),
    R₀ = round(R₀; digits = 4),
    trace = round(trace_value; digits = 6),
    determinant = round(det_value; digits = 8),
    routh_hurwitz_stable = routh_hurwitz,
    approach_is_oscillatory = oscillatory,
)

# ╔═╡ Cell order:
# ╠═e6e72ddd-0fc4-4798-b9f1-f811a2b8702f
# ╠═d0e9244d-7bb1-43fc-87a1-f79aee62739f
# ╟─eb1d5d3d-b8e6-48cf-96d9-6dd7063dc652
# ╠═40d9ff7e-8f42-46e8-a3ab-41e7c07dacb6
# ╠═9fdb247f-c197-497c-a46d-8d2151652a81
# ╠═f103b25c-f003-44a2-a4ef-3aea3bf0b9aa
# ╠═fb781900-b2dc-455c-89dc-c5f3863c481b
# ╠═3b7a8970-4c7d-4411-b844-b317ad6c552e
# ╠═508b3437-be60-4426-bd6f-d0fc179c75c3
# ╠═8e62c4a4-7fb2-43a9-9b0d-bb9f003c3278
