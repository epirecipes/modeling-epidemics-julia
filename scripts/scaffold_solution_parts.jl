using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
units = TOML.parsefile(joinpath(ROOT, "course-units.toml"))["units"]
units_by_id = Dict(unit["id"] => unit for unit in units)

groups = [
    ("foundations", "Prerequisite and foundational solutions", ["p01", "p02", "p03", "ch01", "ch02", "ch03", "ch04"]),
    ("stochastic", "Stochastic modeling solutions", ["ch05", "ch06", "ch07", "ch08", "ch09", "ch10"]),
    ("inference", "Observation and inference solutions", ["ch11", "ch12", "ch13"]),
    ("analysis-capstone", "Analysis, decisions, and capstone solutions", ["ch14", "ch15", "ch16", "ch17"]),
    ("advanced-sciml", "Advanced scientific machine learning solutions", ["a01", "a02"]),
    ("advanced-decisions", "Advanced decision and composition solutions", ["a03", "a04", "a05"]),
]

output_dir = joinpath(ROOT, "appendices", "solutions")
mkpath(output_dir)

for (slug, title, ids) in groups
    output = IOBuffer()
    println(output, "## $title")
    for id in ids
        unit = units_by_id[id]
        println(output)
        println(output, "### $(unit["title"])")
        for number in 1:2
            println(output)
            println(output, "#### Solution to Exercise $number {#solution-$id-ex$number}")
            println(output)
            println(output, "Solution pending.")
        end
    end
    write(joinpath(output_dir, "$slug.qmd"), String(take!(output)))
end

println("Created $(length(groups)) solution part files.")
