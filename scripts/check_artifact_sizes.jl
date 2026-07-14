const ROOT = normpath(joinpath(@__DIR__, ".."))

function directory_size(path; excluded = Set{String}())
    isdir(path) || return 0
    total = 0
    for (directory, subdirectories, files) in walkdir(path)
        filter!(name -> !(name in excluded), subdirectories)
        for file in files
            total += filesize(joinpath(directory, file))
        end
    end
    total
end

function megabytes(bytes)
    round(bytes / 1024^2; digits = 2)
end

freeze_bytes =
    directory_size(joinpath(ROOT, "_freeze")) +
    directory_size(joinpath(ROOT, "slides", "_freeze"))
notebook_bytes = directory_size(
    joinpath(ROOT, "notebook-exports");
    excluded = Set([".cache"]),
)

freeze_limit = parse(Int, get(ENV, "FREEZE_SIZE_LIMIT_MB", "100"))
notebook_limit = parse(Int, get(ENV, "NOTEBOOK_EXPORT_SIZE_LIMIT_MB", "250"))

println("Quarto freeze artifacts: $(megabytes(freeze_bytes)) MiB")
println("Pluto HTML exports: $(megabytes(notebook_bytes)) MiB")

errors = String[]
freeze_bytes > freeze_limit * 1024^2 &&
    push!(errors, "Quarto freeze artifacts exceed $(freeze_limit) MiB.")
notebook_bytes > notebook_limit * 1024^2 &&
    push!(errors, "Pluto HTML exports exceed $(notebook_limit) MiB.")

if !isempty(errors)
    foreach(error -> println(stderr, "ERROR: ", error), errors)
    exit(1)
end
