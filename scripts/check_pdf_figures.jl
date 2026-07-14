const ROOT = normpath(joinpath(@__DIR__, ".."))
const SOURCE_DIRS = ["prerequisites", "chapters", "appendices", "slides"]

errors = String[]

for source_dir in SOURCE_DIRS
    for (parent, _, files) in walkdir(joinpath(ROOT, source_dir))
        for file in files
            endswith(file, ".qmd") || continue
            path = joinpath(parent, file)
            lines = readlines(path)
            for (line_number, line) in enumerate(lines)
                if occursin("\\(", line) || occursin("\\)", line) ||
                   strip(line) in ("\\[", "\\]")
                    push!(
                        errors,
                        "$(relpath(path, ROOT)):$line_number uses LaTeX math delimiters unsupported by Pandoc PDF output.",
                    )
                end
                stripped = strip(line)
                if startswith(stripped, "\$\$") &&
                   stripped != "\$\$" &&
                   !startswith(stripped, "\$\$ {#")
                    push!(
                        errors,
                        "$(relpath(path, ROOT)):$line_number must place display-math content after a standalone \$\$ delimiter.",
                    )
                end
            end
            in_julia = false
            block_start = 0
            block_lines = String[]
            for (line_number, line) in enumerate(lines)
                if !in_julia && strip(line) == "```{julia}"
                    in_julia = true
                    block_start = line_number
                    empty!(block_lines)
                elseif in_julia && strip(line) == "```"
                    text = join(block_lines, "\n")
                    if occursin("Figure(", text)
                        label_match = match(r"(?m)^#\|\s*label:\s*(\S+)\s*$", text)
                        if !isnothing(label_match) && !startswith(label_match.captures[1], "fig-")
                            push!(
                                errors,
                                "$(relpath(path, ROOT)):$block_start creates a figure with non-figure label '$(label_match.captures[1])'.",
                            )
                        end
                    end
                    in_julia = false
                elseif in_julia
                    push!(block_lines, line)
                end
            end
        end
    end
end

if !isempty(errors)
    foreach(error -> println(stderr, "ERROR: ", error), errors)
    exit(1)
end

println("PDF figure labels are valid.")
