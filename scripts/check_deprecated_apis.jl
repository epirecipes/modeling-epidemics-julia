const ROOT = normpath(joinpath(@__DIR__, ".."))
const SOURCE_DIRS = ["prerequisites", "chapters", "appendices", "slides", "notebooks", "src", "test"]
const FORBIDDEN_IMPORTS = [
    "DifferentialEquations",
    "DiffEqJump",
    "DiffEqSensitivity",
    "StochasticDelayDiffEq",
    "SimpleDiffEq",
    "SimJulia",
    "StatsPlots",
    "Plots",
]
const FORBIDDEN_PATTERNS = [
    r"^\s*(?:\w+\s*=\s*)?DiffEqFlux\.(?:FastChain|FastDense|sciml_train)\s*\(",
    r"^\s*(?:\w+\s*=\s*)?(?:FastChain|FastDense)\s*\(",
    r"@eval\s+FiniteStateProjection",
]

errors = String[]

for source_dir in SOURCE_DIRS
    directory = joinpath(ROOT, source_dir)
    isdir(directory) || continue
    for (parent, _, files) in walkdir(directory)
        for file in files
            (endswith(file, ".qmd") || endswith(file, ".jl")) || continue
            path = joinpath(parent, file)
            for (line_number, line) in enumerate(eachline(path))
                stripped = strip(line)
                if startswith(stripped, "using ") || startswith(stripped, "import ")
                    for package in FORBIDDEN_IMPORTS
                        occursin(Regex("\\b$package\\b"), stripped) &&
                            push!(errors, "$(relpath(path, ROOT)):$line_number imports $package")
                    end
                end
                for pattern in FORBIDDEN_PATTERNS
                    occursin(pattern, line) &&
                        push!(errors, "$(relpath(path, ROOT)):$line_number uses forbidden API pattern $pattern")
                end
            end
        end
    end
end

if !isempty(errors)
    foreach(error -> println(stderr, "ERROR: ", error), errors)
    exit(1)
end

println("No deprecated course API imports or compatibility monkey-patches found.")
