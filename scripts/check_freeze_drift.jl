# Compares the freeze in the working tree against the one committed at HEAD.
# Base64 image payloads are replaced by a placeholder and numbers are compared
# with a relative tolerance, so only content drift is reported.
# Set FREEZE_DRIFT_RTOL / FREEZE_DRIFT_ATOL to change the numeric tolerance.

const ROOT = normpath(joinpath(@__DIR__, ".."))
const TARGETS = ["_freeze", "slides/_freeze"]

# A macOS and a Linux run of the course agree to about 1e-4 on the noisiest
# unit.
const RTOL = parse(Float64, get(ENV, "FREEZE_DRIFT_RTOL", "1e-3"))
const ATOL = parse(Float64, get(ENV, "FREEZE_DRIFT_ATOL", "1e-9"))
const REPORTED_LINES = 4
const LINE_WIDTH = 160

# Training a neural ODE amplifies platform differences: macOS and Linux runs
# disagree by tens of percent on the fitted error metrics and print different
# solver warnings. solutions.qmd includes the SciML answers, so it inherits
# this. Only the presence of these units is checked.
const NON_REPRODUCIBLE = [
    "appendices/04-universal-differential-equations",
    "appendices/solutions",
]

const BASE64_PAYLOAD = r"base64,\s*[A-Za-z0-9+/=]+"
const NUMBER = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"

# Quarto versions differ in how deep a path they record for a unit's figure
# directory, and `freeze: auto` only rewrites the entry when the unit is
# re-executed, so a freeze holds a mix of depths. Compare the unit name only.
const SUPPORTING_ENTRY = r"^(\s*)\"([^\"]*_files)(?:/[^\"]*)?\"(,?)$"

cd(ROOT)

committed_paths() = split(
    read(`git ls-tree -r --name-only HEAD -- $TARGETS`, String), '\n'; keepempty = false)

is_execute_result(path) = endswith(path, "/execute-results/html.json")

# Quarto rewrites the freeze on every render, so anything else committed under
# it is left over from an earlier render and will be deleted by the next one.
is_regenerated_elsewhere(path) = contains(path, "/site_libs/")

function committed_text(path)
    read(`git show HEAD:$path`, String)
end

function normalized_lines(text)
    payloads_stripped = replace(text, BASE64_PAYLOAD => "base64,<image>")
    lines = split(replace(payloads_stripped, "\\n" => '\n'), '\n')
    replace.(lines, SUPPORTING_ENTRY => s"\1\"\2\"\3")
end

function skeleton_and_numbers(line)
    numbers = Float64[]
    skeleton = replace(line, NUMBER => function (match)
        value = tryparse(Float64, match)
        value === nothing && return match
        push!(numbers, value)
        "\0"
    end)
    skeleton, numbers
end

function lines_match(committed, current)
    committed == current && return true
    committed_skeleton, committed_numbers = skeleton_and_numbers(committed)
    current_skeleton, current_numbers = skeleton_and_numbers(current)
    committed_skeleton == current_skeleton || return false
    length(committed_numbers) == length(current_numbers) || return false
    all(isapprox(a, b; rtol = RTOL, atol = ATOL)
        for (a, b) in zip(committed_numbers, current_numbers))
end

function truncated(line)
    length(line) > LINE_WIDTH ? string(first(line, LINE_WIDTH), " […]") : line
end

function drifted_lines(committed, current)
    drifted = Tuple{Int,String,String}[]
    for index in 1:min(length(committed), length(current))
        lines_match(committed[index], current[index]) && continue
        push!(drifted, (index, committed[index], current[index]))
    end
    drifted
end

errors = String[]
stray = String[]
checked = 0
image_only = 0
skipped = 0

for path in committed_paths()
    if !is_execute_result(path)
        is_regenerated_elsewhere(path) && continue
        push!(stray, path)
        continue
    end

    if !isfile(joinpath(ROOT, path))
        push!(errors, "$path was not regenerated; the unit no longer renders.")
        continue
    end

    if any(unit -> contains(path, unit), NON_REPRODUCIBLE)
        global skipped += 1
        continue
    end

    global checked += 1
    committed = normalized_lines(committed_text(path))
    current = normalized_lines(read(joinpath(ROOT, path), String))
    drifted = drifted_lines(committed, current)

    if length(committed) != length(current)
        push!(errors, "$path: $(length(committed)) lines committed, " *
                      "$(length(current)) regenerated.")
    end

    if isempty(drifted)
        length(committed) == length(current) && (global image_only += 1)
        continue
    end

    report = ["$path: $(length(drifted)) line(s) drifted beyond " *
              "rtol=$(RTOL), atol=$(ATOL)."]
    for (index, was, now) in first(drifted, REPORTED_LINES)
        push!(report, "    line $index committed:   $(truncated(was))")
        push!(report, "    line $index regenerated: $(truncated(now))")
    end
    length(drifted) > REPORTED_LINES &&
        push!(report, "    … $(length(drifted) - REPORTED_LINES) more.")
    push!(errors, join(report, '\n'))
end

if !isempty(stray)
    directories = sort(unique(dirname.(stray)))
    listed = join(first(directories, 5), "\n    ")
    push!(errors, "$(length(stray)) file(s) committed under the freeze are not " *
                  "execution results; a full render deletes them:\n    " * listed *
                  (length(directories) > 5 ? "\n    …" : "") *
                  "\n    Remove them with: git rm -r --cached <directory>")
end

println("Freeze results compared: $(checked)")
println("Unchanged apart from image bytes: $(image_only)")
println("Present but not compared: $(skipped) ($(join(NON_REPRODUCIBLE, ", ")))")

if !isempty(errors)
    foreach(error -> println(stderr, "ERROR: ", error), errors)
    exit(1)
end

println("No freeze content drift.")
