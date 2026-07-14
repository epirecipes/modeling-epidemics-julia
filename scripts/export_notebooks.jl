using PlutoSliderServer

const ROOT = normpath(joinpath(@__DIR__, ".."))
const NOTEBOOK_DIR = joinpath(ROOT, "notebooks")
const OUTPUT_DIR = joinpath(ROOT, "notebook-exports")
const CACHE_DIR = joinpath(OUTPUT_DIR, ".cache")

parallel_tasks = parse(Int, get(ENV, "PLUTO_EXPORT_TASKS", "1"))

PlutoSliderServer.export_directory(
    NOTEBOOK_DIR;
    Export_output_dir = OUTPUT_DIR,
    Export_cache_dir = CACHE_DIR,
    Export_create_index = false,
    Export_offer_binder = false,
    Export_number_of_parallel_tasks = parallel_tasks,
)
