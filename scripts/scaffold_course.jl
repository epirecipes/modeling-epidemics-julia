using TOML
using UUIDs

const ROOT = normpath(joinpath(@__DIR__, ".."))
const UNIT_CONFIG = joinpath(ROOT, "course-units.toml")
const SECTION_BY_KIND = Dict(
    "prerequisite" => "prerequisites",
    "core" => "chapters",
    "advanced" => "appendices",
)

units = TOML.parsefile(UNIT_CONFIG)["units"]

function unit_stem(unit)
    number = replace(unit["id"], r"^[a-z]+" => "")
    "$(number)-$(unit["slug"])"
end

function unit_paths(unit)
    section = SECTION_BY_KIND[unit["kind"]]
    stem = unit_stem(unit)
    (
        book = joinpath(ROOT, section, "$stem.qmd"),
        slide = joinpath(ROOT, "slides", section, "$stem.qmd"),
        notebook = joinpath(ROOT, "notebooks", section, "$stem.jl"),
        section = section,
        stem = stem,
    )
end

function write_if_missing(path, content)
    isfile(path) && return false
    mkpath(dirname(path))
    write(path, content)
    true
end

function book_stub(unit, paths)
    id = unit["id"]
    course_unit = "$(id)-$(unit["slug"])"
    title = unit["title"]
    status = unit["status"]
    """
    ---
    course-unit: $course_unit
    status: $status
    ---

    # $title

    ::: {.callout-note}
    This teaching unit is scaffolded and awaiting its executable content pass.
    :::

    ## Learning goals

    Content pending.

    ## Model assumptions and classification

    Content pending.

    ## Equations

    Content pending.

    ## Julia implementation

    Content pending.

    ## Validation

    Content pending.

    ## Interpretation and limitations

    Content pending.

    ## Exercises

    ### Exercise 1 {#$id-ex1}

    Exercise pending.

    ### Exercise 2 {#$id-ex2}

    Exercise pending.

    ## Companion materials

    - [Slides](../slides/$(paths.section)/$(paths.stem).html)
    - [Pluto preview](../notebook-exports/$(paths.section)/$(paths.stem).html)
    - [Pluto source](../notebooks/$(paths.section)/$(paths.stem).jl)
    """
end

function slide_stub(unit)
    course_unit = "$(unit["id"])-$(unit["slug"])"
    """
    ---
    title: "$(unit["title"])"
    course-unit: $course_unit
    status: $(unit["status"])
    ---

    ## Why this unit matters

    Content pending.

    ## Model structure

    Content pending.

    ## Julia API

    Content pending.

    ## Key takeaways

    Content pending.

    ::: {.notes}
    Speaker notes pending.
    :::
    """
end

function notebook_stub(unit)
    activation_id, imports_id, title_id, control_id, output_id = [string(uuid4()) for _ in 1:5]
    course_unit = "$(unit["id"])-$(unit["slug"])"
    title = unit["title"]
    status = unit["status"]
    """
    ### A Pluto.jl notebook ###
    # v1.0.3

    using Markdown
    using InteractiveUtils

    # ╔═╡ $activation_id
    begin
        import Pkg
        Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io = devnull)
    end

    # ╔═╡ $imports_id
    using PlutoUI

    # ╔═╡ $title_id
    md\"\"\"
    # $title

    This companion notebook is scaffolded for the `$course_unit` teaching unit.
    \"\"\"

    # ╔═╡ $control_id
    @bind exploration Slider(1:10; default = 5, show_value = true)

    # ╔═╡ $output_id
    (
        course_unit = "$course_unit",
        status = "$status",
        exploration = exploration,
    )

    # ╔═╡ Cell order:
    # ╠═$activation_id
    # ╠═$imports_id
    # ╟─$title_id
    # ╠═$control_id
    # ╠═$output_id
    """
end

created = String[]
for unit in units
    paths = unit_paths(unit)
    write_if_missing(paths.book, book_stub(unit, paths)) && push!(created, relpath(paths.book, ROOT))
    write_if_missing(paths.slide, slide_stub(unit)) && push!(created, relpath(paths.slide, ROOT))
    write_if_missing(paths.notebook, notebook_stub(unit)) && push!(created, relpath(paths.notebook, ROOT))
end

solutions = IOBuffer()
println(solutions, "# Exercise solutions")
println(solutions)
println(solutions, "Solutions are completed during the exercise authoring pass.")
for unit in units
    id = unit["id"]
    println(solutions)
    println(solutions, "## $(unit["title"])")
    println(solutions)
    println(solutions, "### Solution to Exercise 1 {#solution-$id-ex1}")
    println(solutions)
    println(solutions, "Solution pending.")
    println(solutions)
    println(solutions, "### Solution to Exercise 2 {#solution-$id-ex2}")
    println(solutions)
    println(solutions, "Solution pending.")
end
solutions_path = joinpath(ROOT, "appendices", "solutions.qmd")
write_if_missing(solutions_path, String(take!(solutions))) && push!(created, relpath(solutions_path, ROOT))

resources = IOBuffer()
println(resources, "# Course resources {.unnumbered}")
println(resources)
println(resources, "Each teaching unit has a book page, slide deck, Pluto source, and static notebook preview.")
for unit in units
    paths = unit_paths(unit)
    println(resources)
    println(resources, "## $(unit["title"])")
    println(resources)
    println(resources, "- [Book]($(paths.section)/$(paths.stem).html)")
    println(resources, "- [Slides](slides/$(paths.section)/$(paths.stem).html)")
    println(resources, "- [Pluto preview](notebook-exports/$(paths.section)/$(paths.stem).html)")
    println(resources, "- [Pluto source](notebooks/$(paths.section)/$(paths.stem).jl)")
end
resources_path = joinpath(ROOT, "resources.qmd")
write_if_missing(resources_path, String(take!(resources))) && push!(created, relpath(resources_path, ROOT))

println("Created $(length(created)) files.")
foreach(path -> println("  ", path), created)
