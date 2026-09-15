#!/usr/bin/env bash
# Runtime repair, run as the devcontainer's postCreateCommand.
#
# The published image ships a depot built from Project.toml and Manifest.toml as
# they stood when the image was last built. If the manifest has moved on since,
# this brings the codespace's packages back into line; when it has not, it costs
# a few seconds. The expensive build itself lives in build-depot.sh.
set -euo pipefail

cd "$(dirname "$0")/.."

# Precompilation of this stack is memory-hungry; the default of one task per
# core can exhaust a small machine.
export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-2}"

julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=environments/ch10 -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
