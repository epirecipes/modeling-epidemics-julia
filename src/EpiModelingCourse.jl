module EpiModelingCourse

using CairoMakie
using CSV
using DataFrames
using Distributions
using OrdinaryDiffEq
using StableRNGs

export SIRParameters,
    SIRState,
    BASELINE_PARAMETERS,
    BASELINE_SIR_PARAMETERS,
    BASELINE_INITIAL_STATE,
    BASELINE_STATE,
    rate_to_proportion,
    sir_ode!,
    sir_rhs!,
    sir_cumulative_ode!,
    sir_cumulative_rhs!,
    synthetic_incidence,
    generate_synthetic_incidence,
    load_eyam_data,
    COURSE_PALETTE,
    course_palette,
    course_theme

"""
    SIRParameters(; beta=0.05, contact_rate=10.0, gamma=0.25, N=1000)

Parameters for the homogeneous SIR model. `beta` is the per-contact
transmission probability, `contact_rate` is the contact rate, `gamma` is
the recovery rate, and `N` is the population size.
"""
struct SIRParameters
    beta::Float64
    contact_rate::Float64
    gamma::Float64
    N::Float64
end

function Base.getproperty(parameters::SIRParameters, name::Symbol)
    if name === :β
        return getfield(parameters, :beta)
    elseif name === :c
        return getfield(parameters, :contact_rate)
    elseif name === :γ
        return getfield(parameters, :gamma)
    end
    return getfield(parameters, name)
end

function _validate_parameters(parameters::SIRParameters)
    values = (
        parameters.beta,
        parameters.contact_rate,
        parameters.gamma,
        parameters.N,
    )
    all(isfinite, values) ||
        throw(ArgumentError("SIR parameters must be finite"))
    0 <= parameters.beta <= 1 ||
        throw(ArgumentError("beta must be between 0 and 1"))
    parameters.contact_rate >= 0 ||
        throw(ArgumentError("contact_rate must be nonnegative"))
    parameters.gamma >= 0 ||
        throw(ArgumentError("gamma must be nonnegative"))
    parameters.N > 0 ||
        throw(ArgumentError("N must be positive"))
    return nothing
end

function SIRParameters(
    beta::Real,
    contact_rate::Real,
    gamma::Real,
    N::Real,
)
    parameters = SIRParameters(
        Float64(beta),
        Float64(contact_rate),
        Float64(gamma),
        Float64(N),
    )
    _validate_parameters(parameters)
    return parameters
end

function SIRParameters(;
    beta::Real = 0.05,
    contact_rate::Real = 10.0,
    gamma::Real = 0.25,
    N::Real = 1000,
)
    return SIRParameters(beta, contact_rate, gamma, N)
end

"""
    SIRState(; S=990, I=10, R=0)

Initial or current SIR state. The fields are the susceptible, infectious,
and recovered population counts. The model kernels operate on ordinary
three-element vectors in the order `[S, I, R]`; this type is a readable way
to construct such a state.
"""
struct SIRState
    S::Float64
    I::Float64
    R::Float64
end

function Base.getproperty(state::SIRState, name::Symbol)
    if name === :s
        return getfield(state, :S)
    elseif name === :i
        return getfield(state, :I)
    elseif name === :r
        return getfield(state, :R)
    end
    return getfield(state, name)
end

function _validate_state(state::SIRState)
    values = (state.S, state.I, state.R)
    all(isfinite, values) ||
        throw(ArgumentError("SIR state values must be finite"))
    all(>=(0), values) ||
        throw(ArgumentError("SIR state values must be nonnegative"))
    return nothing
end

function SIRState(S::Real, I::Real, R::Real)
    state = SIRState(Float64(S), Float64(I), Float64(R))
    _validate_state(state)
    return state
end

function SIRState(; S::Real = 990, I::Real = 10, R::Real = 0)
    return SIRState(S, I, R)
end

const BASELINE_PARAMETERS = SIRParameters()
const BASELINE_SIR_PARAMETERS = BASELINE_PARAMETERS
const BASELINE_INITIAL_STATE = SIRState()
const BASELINE_STATE = BASELINE_INITIAL_STATE

"""
    rate_to_proportion(rate, dt)

Convert a constant continuous-time rate into the probability of an event
over an interval of length `dt`, `1 - exp(-rate * dt)`.
"""
function rate_to_proportion(rate::Real, dt::Real)
    isfinite(rate) && isfinite(dt) ||
        throw(ArgumentError("rate and dt must be finite"))
    rate >= 0 || throw(ArgumentError("rate must be nonnegative"))
    dt >= 0 || throw(ArgumentError("dt must be nonnegative"))
    return -expm1(-Float64(rate) * Float64(dt))
end

function _parameters(p::SIRParameters)
    return p
end

function _parameters(p::NamedTuple)
    hasproperty(p, :beta) || hasproperty(p, :β) ||
        throw(ArgumentError("parameter tuple must contain beta or β"))
    hasproperty(p, :contact_rate) || hasproperty(p, :c) ||
        throw(ArgumentError("parameter tuple must contain contact_rate or c"))
    hasproperty(p, :gamma) || hasproperty(p, :γ) ||
        throw(ArgumentError("parameter tuple must contain gamma or γ"))
    hasproperty(p, :N) || throw(ArgumentError("parameter tuple must contain N"))
    beta = hasproperty(p, :beta) ? p.beta : p.β
    contact_rate = hasproperty(p, :contact_rate) ? p.contact_rate : p.c
    gamma = hasproperty(p, :gamma) ? p.gamma : p.γ
    return SIRParameters(beta, contact_rate, gamma, p.N)
end

function _parameters(p::Tuple)
    length(p) == 4 ||
        throw(ArgumentError("positional parameters must be (beta, contact_rate, gamma, N)"))
    return SIRParameters(p[1], p[2], p[3], p[4])
end

function _parameters(p)
    throw(ArgumentError("parameters must be SIRParameters or a four-value tuple"))
end

function _state(state::SIRState)
    return state
end

function _state(state::NamedTuple)
    hasproperty(state, :S) || hasproperty(state, :s) ||
        throw(ArgumentError("state tuple must contain S or s"))
    hasproperty(state, :I) || hasproperty(state, :i) ||
        throw(ArgumentError("state tuple must contain I or i"))
    hasproperty(state, :R) || hasproperty(state, :r) ||
        throw(ArgumentError("state tuple must contain R or r"))
    S = hasproperty(state, :S) ? state.S : state.s
    I = hasproperty(state, :I) ? state.I : state.i
    R = hasproperty(state, :R) ? state.R : state.r
    return SIRState(S, I, R)
end

function _state(state::Tuple)
    length(state) == 3 ||
        throw(ArgumentError("positional state must be (S, I, R)"))
    return SIRState(state[1], state[2], state[3])
end

function _state(state)
    state isa AbstractVector && length(state) == 3 ||
        throw(ArgumentError("initial_state must be an SIRState or three values"))
    return SIRState(state[1], state[2], state[3])
end

@inline function _rates(p)
    parameters = _parameters(p)
    transmission = parameters.beta * parameters.contact_rate
    return transmission, parameters.gamma, parameters.N
end

"""
    sir_ode!(du, u, p, t)

In-place right-hand side for the standard SIR ODE:
`S' = -beta * contact_rate * S * I / N`,
`I' = beta * contact_rate * S * I / N - gamma * I`, and
`R' = gamma * I`.

`p` may be a `SIRParameters`, a named tuple with fields
`(beta, contact_rate, gamma, N)`, or the corresponding four-value tuple.
"""
function sir_ode!(du, u, p, t)
    transmission, gamma, N = _rates(p)
    S, I = u[1], u[2]
    incidence = transmission * S * I / N
    du[1] = -incidence
    du[2] = incidence - gamma * I
    du[3] = gamma * I
    return nothing
end

const sir_rhs! = sir_ode!

"""
    sir_cumulative_ode!(du, u, p, t)

In-place SIR right-hand side with a fourth state `C` for cumulative
incidence. The first three derivatives are those of [`sir_ode!`](@ref), and
`C' = beta * contact_rate * S * I / N`.
"""
function sir_cumulative_ode!(du, u, p, t)
    transmission, gamma, N = _rates(p)
    S, I = u[1], u[2]
    incidence = transmission * S * I / N
    du[1] = -incidence
    du[2] = incidence - gamma * I
    du[3] = gamma * I
    du[4] = incidence
    return nothing
end

const sir_cumulative_rhs! = sir_cumulative_ode!

const COURSE_PALETTE = (
    susceptible = :steelblue,
    infected = :firebrick,
    recovered = :seagreen,
    incidence = :darkorange,
)

"""Return the immutable default palette used by course figures."""
course_palette() = COURSE_PALETTE

"""
    course_theme(; palette=COURSE_PALETTE)

Construct a CairoMakie `Theme` for course figures. Constructing the theme
does not change Makie's global theme; pass the result to `with_theme` or
`set_theme` explicitly at the call site.
"""
function course_theme(; palette = COURSE_PALETTE)
    colors = [
        palette.susceptible,
        palette.infected,
        palette.recovered,
        palette.incidence,
    ]
    return Theme(
        palette = (color = colors,),
        Axis = (
            xlabelsize = 18,
            ylabelsize = 18,
            titlesize = 20,
        ),
    )
end

function _time_grid(tspan, saveat)
    tspan isa Tuple && length(tspan) == 2 ||
        throw(ArgumentError("tspan must be a two-element tuple"))
    t0, t1 = tspan
    t0 isa Real && t1 isa Real ||
        throw(ArgumentError("tspan values must be real"))
    isfinite(t0) && isfinite(t1) ||
        throw(ArgumentError("tspan values must be finite"))
    t1 > t0 || throw(ArgumentError("tspan must have t1 > t0"))

    if saveat isa Real
        isfinite(saveat) && saveat > 0 ||
            throw(ArgumentError("saveat must be a positive finite number"))
        times = collect(t0:saveat:t1)
        isempty(times) && (times = [t0])
        times[end] < t1 && push!(times, t1)
        return Float64.(times), (Float64(t0), Float64(t1))
    elseif saveat isa AbstractVector
        isempty(saveat) && throw(ArgumentError("saveat cannot be empty"))
        all(x -> x isa Real && isfinite(x), saveat) ||
            throw(ArgumentError("saveat values must be finite real numbers"))
        times = Float64.(saveat)
        times[1] == t0 ||
            throw(ArgumentError("saveat must start at tspan[1]"))
        times[end] <= t1 ||
            throw(ArgumentError("saveat cannot exceed tspan[2]"))
        all(diff(times) .> 0) ||
            throw(ArgumentError("saveat values must be strictly increasing"))
        return times, (Float64(t0), Float64(t1))
    end
    throw(ArgumentError("saveat must be a positive number or vector"))
end

function _observed_counts(rng, expected::Vector{Float64}, observation_noise::Symbol)
    if observation_noise === :poisson
        return [rand(rng, Poisson(mean)) for mean in expected]
    elseif observation_noise === :none
        return round.(Int, expected)
    end
    throw(ArgumentError("observation_noise must be :poisson or :none"))
end

"""
    synthetic_incidence(; kwargs...)

Solve the cumulative-incidence SIR ODE and sample reproducible observations.
The returned `DataFrame` has the stable schema
`time`, `S`, `I`, `R`, `cumulative_incidence`, `incidence`,
`expected_incidence`, and `observed_incidence`. `incidence` is the latent
increment between consecutive observation times (zero in the first row);
`observed_incidence` is a Poisson count by default. Set
`observation_noise=:none` to inspect rounded noiseless observations.

The random stream is local to this function (`StableRNG(seed)`), so global
random state is not changed.
"""
function synthetic_incidence(;
    seed::Integer = 1,
    parameters = BASELINE_PARAMETERS,
    initial_state = BASELINE_INITIAL_STATE,
    tspan = (0.0, 30.0),
    saveat = 1.0,
    observation_noise::Symbol = :poisson,
    reporting_rate::Real = 1.0,
    abstol::Real = 1e-8,
    reltol::Real = 1e-8,
)
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    parameters = _parameters(parameters)
    initial_state = _state(initial_state)
    isapprox(
        initial_state.S + initial_state.I + initial_state.R,
        parameters.N;
        atol = 1e-8,
        rtol = 0,
    ) || throw(ArgumentError("initial state must sum to N"))
    isfinite(reporting_rate) && reporting_rate >= 0 ||
        throw(ArgumentError("reporting_rate must be finite and nonnegative"))
    isfinite(abstol) && abstol > 0 || throw(ArgumentError("abstol must be positive"))
    isfinite(reltol) && reltol > 0 || throw(ArgumentError("reltol must be positive"))
    times, span = _time_grid(tspan, saveat)

    u0 = [initial_state.S, initial_state.I, initial_state.R, 0.0]
    problem = ODEProblem(sir_cumulative_ode!, u0, span, parameters)
    solution = solve(problem, Tsit5(); saveat = times, abstol = abstol, reltol = reltol)

    states = Array(solution)
    cumulative = Float64.(states[4, :])
    incidence = [0.0; diff(cumulative)]
    expected = reporting_rate .* incidence
    observed = _observed_counts(StableRNG(seed), expected, observation_noise)

    return DataFrame(
        time = Float64.(solution.t),
        S = Float64.(states[1, :]),
        I = Float64.(states[2, :]),
        R = Float64.(states[3, :]),
        cumulative_incidence = cumulative,
        incidence = incidence,
        expected_incidence = expected,
        observed_incidence = observed,
    )
end

const generate_synthetic_incidence = synthetic_incidence

"""
    load_eyam_data([path])

Load and validate the committed Eyam plague capstone data. The returned
`DataFrame` has columns `time_months`, `S`, `I`, and `R`; every row must
conserve the documented population of 261.
"""
function load_eyam_data(
    path::AbstractString = normpath(joinpath(@__DIR__, "..", "data", "eyam_plague_1666.csv")),
)
    isfile(path) || throw(ArgumentError("Eyam data file does not exist: $path"))
    data = CSV.read(path, DataFrame)
    names(data) == ["time_months", "S", "I", "R"] ||
        throw(ArgumentError("Eyam data must have columns time_months, S, I, and R"))
    nrow(data) == 8 || throw(ArgumentError("Eyam data must contain 8 observations"))
    all(diff(data.time_months) .> 0) ||
        throw(ArgumentError("Eyam observation times must be strictly increasing"))
    all(data.S .+ data.I .+ data.R .== 261) ||
        throw(ArgumentError("Eyam observations must conserve N = 261"))
    return data
end

end
