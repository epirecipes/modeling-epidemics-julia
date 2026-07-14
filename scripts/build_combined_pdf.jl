using TOML
using SHA

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PDF_ROOT = joinpath(ROOT, "_pdf")
const DIST_DIR = joinpath(ROOT, "dist")
const BOOK_PDF = joinpath(PDF_ROOT, "book", "modeling-epidemics-with-julia.pdf")
const FINAL_PDF = joinpath(DIST_DIR, "modeling-epidemics-with-julia-complete.pdf")
const SECTION_BY_KIND = Dict(
    "prerequisite" => "prerequisites",
    "core" => "chapters",
    "advanced" => "appendices",
)

function required_tool(name)
    path = Sys.which(name)
    isnothing(path) && error("Required PDF tool is not available: $name")
    path
end

function unit_stem(unit)
    number = replace(unit["id"], r"^[a-z]+" => "")
    "$(number)-$(unit["slug"])"
end

function cleanup_source_artifacts()
    generated_extensions = Set([
        ".html",
        ".log",
        ".tex",
        ".aux",
        ".nav",
        ".snm",
        ".toc",
        ".out",
        ".vrb",
    ])
    for source_dir in ("prerequisites", "chapters", "appendices", "slides", "notebooks")
        for (parent, subdirectories, files) in walkdir(joinpath(ROOT, source_dir))
            filter!(subdirectories) do directory
                if startswith(directory, ".") || directory == "_freeze"
                    return false
                elseif endswith(directory, "_files")
                    rm(joinpath(parent, directory); recursive = true)
                    return false
                end
                true
            end
            for file in files
                splitext(file)[2] in generated_extensions && rm(joinpath(parent, file))
            end
        end
        for path in (
            joinpath(ROOT, "index_files"),
            joinpath(ROOT, "modeling-epidemics-with-julia.tex"),
        )
            ispath(path) && rm(path; recursive = isdir(path))
        end
    end
end

skip_render = "--skip-render" in ARGS
unknown_args = setdiff(ARGS, ["--skip-render"])
isempty(unknown_args) || error("Unknown arguments: $(join(unknown_args, ", "))")

if !skip_render
    cd(ROOT) do
        run(`quarto render --profile pdf --to pdf`)
    end
end

isfile(BOOK_PDF) || error("Book PDF was not generated: $BOOK_PDF")
units = TOML.parsefile(joinpath(ROOT, "course-units.toml"))["units"]
slide_pdfs = String[]
for unit in units
    section = SECTION_BY_KIND[unit["kind"]]
    path = joinpath(PDF_ROOT, "slides", section, "$(unit_stem(unit)).pdf")
    isfile(path) || error("Slide PDF was not generated: $path")
    push!(slide_pdfs, path)
end
length(slide_pdfs) == 24 || error("Expected 24 slide PDFs, found $(length(slide_pdfs))")

pdfunite = required_tool("pdfunite")
ghostscript = required_tool("gs")
mkpath(DIST_DIR)
mkpath(PDF_ROOT)

uncompressed = joinpath(PDF_ROOT, "modeling-epidemics-with-julia-complete-uncompressed.pdf")
isfile(uncompressed) && rm(uncompressed)
run(Cmd(vcat([pdfunite], [BOOK_PDF], slide_pdfs, [uncompressed])))

temporary_output = joinpath(PDF_ROOT, "modeling-epidemics-with-julia-complete.pdf")
isfile(temporary_output) && rm(temporary_output)
run(
    `$ghostscript -sDEVICE=pdfwrite -dCompatibilityLevel=1.7 -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages=true -dCompressFonts=true -sOutputFile=$temporary_output $uncompressed`,
)

isfile(temporary_output) || error("Ghostscript did not produce the combined PDF.")
filesize(temporary_output) > 0 || error("Combined PDF is empty.")
mv(temporary_output, FINAL_PDF; force = true)

checksum = bytes2hex(sha256(read(FINAL_PDF)))
checksum_path = "$FINAL_PDF.sha256"
write(checksum_path, "$checksum  $(basename(FINAL_PDF))\n")
cleanup_source_artifacts()

println("Created $(relpath(FINAL_PDF, ROOT)) ($(round(filesize(FINAL_PDF) / 1024^2; digits = 1)) MiB)")
println("SHA-256: $checksum")
