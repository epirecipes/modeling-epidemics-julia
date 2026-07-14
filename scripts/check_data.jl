using CSV
using DataFrames
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CSV_PATH = joinpath(ROOT, "data", "eyam_plague_1666.csv")
const METADATA_PATH = joinpath(ROOT, "data", "eyam_plague_1666.toml")

metadata = TOML.parsefile(METADATA_PATH)
data = CSV.read(CSV_PATH, DataFrame)

expected_columns = ["time_months", "S", "I", "R"]
String.(names(data)) == expected_columns ||
    error("Unexpected Eyam data columns: $(names(data)).")
nrow(data) == 8 || error("Expected 8 Eyam observations, found $(nrow(data)).")
all(diff(data.time_months) .> 0) || error("Eyam observation times must increase.")
all(data.S .+ data.I .+ data.R .== metadata["population"]) ||
    error("Every Eyam row must conserve the modeled population.")

actual_checksum = bytes2hex(sha256(read(CSV_PATH)))
expected_checksum = metadata["csv_sha256"]
actual_checksum == expected_checksum ||
    error("Eyam CSV checksum differs from its provenance metadata.")

println("Capstone data is valid: 8 observations, N=$(metadata["population"]), checksum verified.")
