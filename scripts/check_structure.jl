using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const VALID_STATUSES = Set(["stub", "draft", "complete"])
const SECTION_BY_KIND = Dict(
    "prerequisite" => "prerequisites",
    "core" => "chapters",
    "advanced" => "appendices",
)
const REQUIRED_BOOK_HEADINGS = [
    "## Learning goals",
    "## Model assumptions and classification",
    "## Equations",
    "## Julia implementation",
    "## Validation",
    "## Interpretation and limitations",
    "## Exercises",
    "## Companion materials",
]
const VALID_TRACKS = Set(["workshop", "self-study"])
const WORKSHOP_TIME_LIMIT_MINUTES = parse(Int, get(ENV, "WORKSHOP_TIME_LIMIT_MINUTES", "1000"))

units = TOML.parsefile(joinpath(ROOT, "course-units.toml"))["units"]
errors = String[]

length(units) == 24 || push!(errors, "Expected 24 teaching units, found $(length(units)).")

function unit_stem(unit)
    number = replace(unit["id"], r"^[a-z]+" => "")
    "$(number)-$(unit["slug"])"
end

function read_required(path)
    if !isfile(path)
        push!(errors, "Missing $(relpath(path, ROOT)).")
        return ""
    end
    read(path, String)
end

function capture_status(text, pattern, artifact)
    match_result = match(pattern, text)
    if isnothing(match_result)
        push!(errors, "$artifact is missing status metadata.")
        return ""
    end
    status = match_result.captures[1]
    status in VALID_STATUSES || push!(errors, "$artifact has invalid status '$status'.")
    status
end

solution_dir = joinpath(ROOT, "appendices", "solutions")
solution_paths = isdir(solution_dir) ? sort(filter(path -> endswith(path, ".qmd"), readdir(solution_dir; join = true))) : String[]
length(solution_paths) == 7 || push!(errors, "Expected 7 solution part files, found $(length(solution_paths)).")
solutions = join(read_required.(solution_paths), "\n")
resources = read_required(joinpath(ROOT, "resources.qmd"))
workshop_page = read_required(joinpath(ROOT, "workshop.qmd"))
seen_ids = Set{String}()
workshop_slots = Dict{Tuple{Int, Int}, String}()
workshop_minutes = Int[]

for unit in units
    id = unit["id"]
    slug = unit["slug"]
    full_id = "$id-$slug"
    overall_status = unit["status"]
    section = SECTION_BY_KIND[unit["kind"]]
    stem = unit_stem(unit)

    id in seen_ids && push!(errors, "Duplicate teaching unit ID: $id.")
    push!(seen_ids, id)
    overall_status in VALID_STATUSES || push!(errors, "$id has invalid overall status '$overall_status'.")

    track = get(unit, "track", nothing)
    if isnothing(track)
        push!(errors, "$id is missing a track field.")
    elseif !(track in VALID_TRACKS)
        push!(errors, "$id has invalid track '$track'.")
    elseif track == "workshop"
        for key in ("day", "session", "time_minutes")
            haskey(unit, key) || push!(errors, "$id is a workshop unit but has no '$key'.")
        end
        if haskey(unit, "day") && haskey(unit, "session")
            slot = (unit["day"], unit["session"])
            if haskey(workshop_slots, slot)
                push!(errors, "$id and $(workshop_slots[slot]) both claim day $(slot[1]) session $(slot[2]).")
            else
                workshop_slots[slot] = id
            end
        end
        haskey(unit, "time_minutes") && push!(workshop_minutes, unit["time_minutes"])
        occursin("$(section)/$(stem).html", workshop_page) ||
            push!(errors, "$id is a workshop unit but has no link on workshop.qmd.")
    else
        for key in ("day", "session", "time_minutes")
            haskey(unit, key) && push!(errors, "$id is self-study but has a workshop-only field '$key'.")
        end
    end

    book_path = joinpath(ROOT, section, "$stem.qmd")
    slide_path = joinpath(ROOT, "slides", section, "$stem.qmd")
    notebook_path = joinpath(ROOT, "notebooks", section, "$stem.jl")

    book = read_required(book_path)
    slide = read_required(slide_path)
    notebook = read_required(notebook_path)
    book_status = capture_status(book, r"(?m)^status:\s*(stub|draft|complete)\s*$", "$id book")
    slide_status = capture_status(slide, r"(?m)^status:\s*(stub|draft|complete)\s*$", "$id slide")
    notebook_status = capture_status(notebook, r"status\s*=\s*\"(stub|draft|complete)\"", "$id notebook")

    occursin("course-unit: $full_id", book) || push!(errors, "$id book metadata has the wrong course-unit.")
    occursin("course-unit: $full_id", slide) || push!(errors, "$id slide metadata has the wrong course-unit.")
    occursin("course_unit = \"$full_id\"", notebook) || push!(errors, "$id notebook metadata has the wrong course-unit.")
    occursin("Pkg.activate", notebook) || push!(errors, "$id notebook does not activate the shared course environment.")

    # The opening callout on every book page is a further copy of the schedule, so it
    # drifts silently on a reorder unless it is checked against course-units.toml.
    expected_callout = if track == "workshop" && haskey(unit, "day") && haskey(unit, "session")
        "**Workshop programme**: taught live on Day $(unit["day"]), Session $(unit["session"])."
    elseif track == "self-study"
        "**Self-study**: not taught at the workshop"
    end
    isnothing(expected_callout) || occursin(expected_callout, book) ||
        push!(errors, "$id book callout contradicts course-units.toml; expected \"$expected_callout\".")

    if overall_status == "complete" && (book_status != "complete" || slide_status != "complete" || notebook_status != "complete")
        push!(errors, "$id is overall complete but one or more companion artifacts are incomplete.")
    end

    for heading in REQUIRED_BOOK_HEADINGS
        occursin(heading, book) || push!(errors, "$id book is missing '$heading'.")
    end

    exercise_matches = collect(eachmatch(Regex("\\{#$id-ex[0-9]+\\}"), book))
    length(exercise_matches) == 2 || push!(errors, "$id must define exactly two exercise anchors.")
    for number in 1:2
        occursin("{#$id-ex$number}", book) || push!(errors, "$id is missing exercise $number.")
        occursin("{#solution-$id-ex$number}", solutions) || push!(errors, "$id is missing solution $number.")
    end

    occursin("$(section)/$(stem).html", resources) || push!(errors, "$id book link is missing from resources.qmd.")
    occursin("slides/$(section)/$(stem).html", resources) || push!(errors, "$id slide link is missing from resources.qmd.")
    occursin("notebook-exports/$(section)/$(stem).html", resources) || push!(errors, "$id notebook preview link is missing from resources.qmd.")

    book_status == "complete" && occursin(r"\bpending\b", lowercase(book)) &&
        push!(errors, "$id book is complete but still contains pending content.")
    book_status == "complete" && occursin(r"\bscaffold(?:ed)?\b", lowercase(book)) &&
        push!(errors, "$id book is complete but still contains scaffold text.")
    slide_status == "complete" && occursin(r"\bpending\b", lowercase(slide)) &&
        push!(errors, "$id slide is complete but still contains pending content.")
    slide_status == "complete" && occursin(r"\bscaffold(?:ed)?\b", lowercase(slide)) &&
        push!(errors, "$id slide is complete but still contains scaffold text.")
    notebook_status == "complete" && occursin(r"\bscaffold(?:ed)?\b", lowercase(notebook)) &&
        push!(errors, "$id notebook is complete but still contains scaffold text.")
end

solution_matches = collect(eachmatch(r"\{#solution-[a-z0-9]+-ex[0-9]+\}", solutions))
length(solution_matches) == 48 || push!(errors, "Expected 48 solution anchors, found $(length(solution_matches)).")

total_workshop_minutes = sum(workshop_minutes; init = 0)
total_workshop_minutes > WORKSHOP_TIME_LIMIT_MINUTES && push!(
    errors,
    "Workshop track totals $(total_workshop_minutes) minutes, over the $(WORKSHOP_TIME_LIMIT_MINUTES)-minute budget.",
)

for source_dir in ("prerequisites", "chapters", "appendices", "slides", "notebooks")
    for (parent, subdirectories, files) in walkdir(joinpath(ROOT, source_dir))
        filter!(directory -> !startswith(directory, ".") && directory != "_freeze", subdirectories)
        for directory in subdirectories
            endswith(directory, "_files") &&
                push!(errors, "Generated directory found in source tree: $(relpath(joinpath(parent, directory), ROOT)).")
        end
        for file in files
            (endswith(file, ".html") || endswith(file, ".log")) &&
                push!(errors, "Generated file found in source tree: $(relpath(joinpath(parent, file), ROOT)).")
        end
    end
end

if !isempty(errors)
    foreach(error -> println(stderr, "ERROR: ", error), errors)
    exit(1)
end

println("Course structure is valid: 24 teaching units and 48 exercise/solution pairs.")
