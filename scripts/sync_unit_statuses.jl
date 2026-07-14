using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIG_PATH = joinpath(ROOT, "course-units.toml")
const SECTION_BY_KIND = Dict(
    "prerequisite" => "prerequisites",
    "core" => "chapters",
    "advanced" => "appendices",
)

units = TOML.parsefile(CONFIG_PATH)["units"]

function unit_stem(unit)
    number = replace(unit["id"], r"^[a-z]+" => "")
    "$(number)-$(unit["slug"])"
end

function capture_status(path, pattern)
    text = read(path, String)
    result = match(pattern, text)
    isnothing(result) && error("Missing status metadata in $(relpath(path, ROOT))")
    result.captures[1]
end

statuses = Dict{String, String}()
for unit in units
    section = SECTION_BY_KIND[unit["kind"]]
    stem = unit_stem(unit)
    book_status = capture_status(
        joinpath(ROOT, section, "$stem.qmd"),
        r"(?m)^status:\s*(stub|draft|complete)\s*$",
    )
    slide_status = capture_status(
        joinpath(ROOT, "slides", section, "$stem.qmd"),
        r"(?m)^status:\s*(stub|draft|complete)\s*$",
    )
    notebook_status = capture_status(
        joinpath(ROOT, "notebooks", section, "$stem.jl"),
        r"status\s*=\s*\"(stub|draft|complete)\"",
    )
    artifact_statuses = [book_status, slide_status, notebook_status]
    statuses[unit["id"]] = all(==("complete"), artifact_statuses) ? "complete" :
        all(==("stub"), artifact_statuses) ? "stub" : "draft"
end

function write_statuses!(statuses)
    lines = readlines(CONFIG_PATH; keep = true)
    current_id = nothing
    for index in eachindex(lines)
        id_match = match(r"^id = \"([^\"]+)\"", lines[index])
        if !isnothing(id_match)
            current_id = id_match.captures[1]
        elseif !isnothing(current_id) && occursin(r"^status = \"", lines[index])
            newline = endswith(lines[index], "\n") ? "\n" : ""
            lines[index] = "status = \"$(statuses[current_id])\"$newline"
        end
    end
    write(CONFIG_PATH, join(lines))
end

write_statuses!(statuses)

counts = Dict(status => count(==(status), values(statuses)) for status in ("stub", "draft", "complete"))
println("Synchronized unit statuses: ", counts)
