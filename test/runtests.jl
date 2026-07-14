using EpiModelingCourse
using Test
using CairoMakie
using DataFrames
using OrdinaryDiffEq

@testset "baseline SIR definitions" begin
    @test nameof(EpiModelingCourse) == :EpiModelingCourse
    @test BASELINE_PARAMETERS == SIRParameters(0.05, 10.0, 0.25, 1000)
    @test BASELINE_INITIAL_STATE == SIRState(990, 10, 0)
    @test BASELINE_INITIAL_STATE.S + BASELINE_INITIAL_STATE.I +
          BASELINE_INITIAL_STATE.R == BASELINE_PARAMETERS.N
    @test BASELINE_SIR_PARAMETERS === BASELINE_PARAMETERS
    @test BASELINE_STATE === BASELINE_INITIAL_STATE
    @test BASELINE_PARAMETERS.β == BASELINE_PARAMETERS.beta
    @test BASELINE_PARAMETERS.c == BASELINE_PARAMETERS.contact_rate
    @test BASELINE_PARAMETERS.γ == BASELINE_PARAMETERS.gamma
    @test BASELINE_INITIAL_STATE.s == BASELINE_INITIAL_STATE.S
end

@testset "rate conversion and ODE kernels" begin
    @test rate_to_proportion(0.25, 1.0) ≈ 1 - exp(-0.25)
    @test rate_to_proportion(0, 0) == 0
    @test_throws ArgumentError rate_to_proportion(-1, 1)
    @test_throws ArgumentError rate_to_proportion(1, -1)
    @test_throws ArgumentError rate_to_proportion(NaN, 1)
    @test_throws ArgumentError SIRParameters(beta = -0.1)
    @test_throws ArgumentError SIRParameters(gamma = Inf)
    @test_throws ArgumentError SIRState(S = -1)

    du = zeros(3)
    sir_ode!(du, [990.0, 10.0, 0.0], BASELINE_PARAMETERS, 0.0)
    @test sum(du) ≈ 0
    @test du[1] < 0
    @test du[3] > 0

    cumulative_du = zeros(4)
    sir_cumulative_ode!(
        cumulative_du,
        [990.0, 10.0, 0.0, 0.0],
        (beta = 0.05, contact_rate = 10.0, gamma = 0.25, N = 1000),
        0.0,
    )
    @test cumulative_du[4] ≈ -cumulative_du[1]
    @test sir_rhs! === sir_ode!
    @test sir_cumulative_rhs! === sir_cumulative_ode!

    problem = ODEProblem(
        sir_ode!,
        [990.0, 10.0, 0.0],
        (0.0, 20.0),
        BASELINE_PARAMETERS,
    )
    solution = solve(problem, Tsit5(); saveat = 1.0)
    @test all(sum(solution.u[i]) ≈ 1000 for i in eachindex(solution.u))
    @test all(x -> all(x .>= 0), solution.u)
end

@testset "synthetic incidence data" begin
    data1 = synthetic_incidence(
        seed = 42,
        tspan = (0.0, 10.0),
        saveat = 1.0,
    )
    data2 = generate_synthetic_incidence(
        seed = 42,
        tspan = (0.0, 10.0),
        saveat = 1.0,
    )
    @test data1 isa DataFrame
    @test names(data1) == [
        "time",
        "S",
        "I",
        "R",
        "cumulative_incidence",
        "incidence",
        "expected_incidence",
        "observed_incidence",
    ]
    @test data1 == data2
    custom = synthetic_incidence(
        seed = 42,
        parameters = (β = 0.05, c = 10.0, γ = 0.25, N = 1000),
        initial_state = (S = 990, I = 10, R = 0),
        tspan = (0.0, 2.0),
        saveat = 1.0,
    )
    @test all(
        custom[!, name] ≈ first(data1, 3)[!, name]
        for name in names(data1)
    )
    @test all(data1.observed_incidence .>= 0)
    @test all(data1.incidence .>= 0)
    @test all(diff(data1.cumulative_incidence) .>= -1e-8)
    @test data1.observed_incidence != synthetic_incidence(
        seed = 43,
        tspan = (0.0, 10.0),
        saveat = 1.0,
    ).observed_incidence

    noiseless = synthetic_incidence(
        seed = 1,
        tspan = (0.0, 2.0),
        saveat = 1.0,
        observation_noise = :none,
    )
    @test noiseless.observed_incidence == round.(Int, noiseless.expected_incidence)

    @test_throws ArgumentError synthetic_incidence(seed = -1)
    @test_throws ArgumentError synthetic_incidence(tspan = (0.0, 0.0))
    @test_throws ArgumentError synthetic_incidence(saveat = 0.0)
    @test_throws ArgumentError synthetic_incidence(reporting_rate = -1)
    @test_throws ArgumentError synthetic_incidence(observation_noise = :normal)
    @test_throws ArgumentError synthetic_incidence(
        parameters = SIRParameters(N = 1001),
    )
end

@testset "course theme" begin
    @test course_palette() === COURSE_PALETTE
    theme = course_theme()
    @test theme isa CairoMakie.Theme
end

@testset "Eyam capstone data" begin
    eyam = load_eyam_data()
    @test names(eyam) == ["time_months", "S", "I", "R"]
    @test size(eyam) == (8, 4)
    @test all(eyam.S .+ eyam.I .+ eyam.R .== 261)
    @test all(diff(eyam.time_months) .> 0)
    @test_throws ArgumentError load_eyam_data("missing-eyam-data.csv")
end
