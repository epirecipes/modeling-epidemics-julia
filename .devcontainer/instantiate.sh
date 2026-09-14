#!/usr/bin/env bash
# Installs and precompiles the course packages.
#
# This runs as the devcontainer's onCreateCommand, which Codespaces also runs
# while building a prebuild, so with prebuilds enabled participants never wait
# for it. Without prebuilds it runs on first start and takes a long time.
set -euo pipefail

cd "$(dirname "$0")/.."

# Precompilation of this stack is memory-hungry; on a 2-core machine the
# default one task per core can exhaust the 8 GB image.
export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-2}"

julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=environments/ch10 -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Load what the early worksheets load, so the first cell a participant runs
# does not pay for loading these into a fresh session's cache.
julia --project=. -e 'using EpiModelingCourse, CairoMakie, OrdinaryDiffEq, DataFrames, StableRNGs; println("Course environment ready on Julia ", VERSION)'
