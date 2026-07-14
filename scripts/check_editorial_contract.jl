using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SECTION_BY_KIND = Dict(
    "prerequisite" => "prerequisites",
    "core" => "chapters",
    "advanced" => "appendices",
)

units = TOML.parsefile(joinpath(ROOT, "course-units.toml"))["units"]
errors = String[]

function unit_stem(unit)
    number = replace(unit["id"], r"^[a-z]+" => "")
    "$(number)-$(unit["slug"])"
end

for unit in units
    section = SECTION_BY_KIND[unit["kind"]]
    stem = unit_stem(unit)
    book_path = joinpath(ROOT, section, "$stem.qmd")
    slide_path = joinpath(ROOT, "slides", section, "$stem.qmd")
    notebook_path = joinpath(ROOT, "notebooks", section, "$stem.jl")
    book = read(book_path, String)
    slide = read(slide_path, String)
    notebook = read(notebook_path, String)

    occursin("#| fig-cap:", book) ||
        push!(errors, "$(unit["id"]) book has no captioned generated figure.")
    occursin("::: {.notes}", slide) ||
        push!(errors, "$(unit["id"]) slide deck has no speaker notes.")
    if occursin("Figure(", slide) && !occursin("#| fig-cap:", slide)
        push!(errors, "$(unit["id"]) slide figure has no caption.")
    end
    occursin("@bind", notebook) ||
        push!(errors, "$(unit["id"]) Pluto notebook has no interactive control.")
end

solution_dir = joinpath(ROOT, "appendices", "solutions")
for path in sort(readdir(solution_dir; join = true))
    endswith(path, ".qmd") || continue
    occursin(r"\bSolution pending\b", read(path, String)) &&
        push!(errors, "$(relpath(path, ROOT)) still contains pending solutions.")
end

bibliography = read(joinpath(ROOT, "references.bib"), String)
bib_keys = Set(
    match_result.captures[1]
    for match_result in eachmatch(r"(?m)^@\w+\{([^,]+),", bibliography)
)
citation_keys = Set{String}()
for directory in ("prerequisites", "chapters", "appendices", "slides")
    for (parent, _, files) in walkdir(joinpath(ROOT, directory))
        for file in files
            endswith(file, ".qmd") || continue
            text = read(joinpath(parent, file), String)
            for match_result in eachmatch(r"@([a-z][a-z0-9_-]*[0-9]{4}[a-z]?)\b", text)
                push!(citation_keys, match_result.captures[1])
            end
        end
    end
end
for key in sort!(collect(setdiff(citation_keys, bib_keys)))
    push!(errors, "Citation @$key has no references.bib entry.")
end

if !isempty(errors)
    foreach(error -> println(stderr, "ERROR: ", error), errors)
    exit(1)
end

println("Editorial contract is valid: figures, notes, interactivity, solutions, and citations.")
