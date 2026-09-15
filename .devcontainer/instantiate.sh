#!/usr/bin/env bash
# Installs and precompiles the course packages.
#
# This runs as the devcontainer's onCreateCommand, which Codespaces also runs
# while building a prebuild, so with prebuilds enabled participants never wait
# for it. Without prebuilds it runs on first start and takes a long time.
set -euo pipefail

cd "$(dirname "$0")/.."

# --no-warmup skips the final load, for the postCreate run where the depot is
# usually already correct and the extra half minute buys nothing.
warmup=1
[ "${1:-}" = "--no-warmup" ] && warmup=0

# Precompilation of this stack is memory-hungry; on a 2-core machine the
# default one task per core can exhaust the 8 GB image.
export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-2}"

julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=environments/ch10 -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Trim anything orphaned, and reset the "we haven't cleaned this depot up for a
# bit" timer so that message does not greet the first render from a prebuild.
julia --project=. -e 'using Pkg; Pkg.gc()' >/dev/null

if [ "$warmup" = 1 ]; then
  # Load what the early worksheets load, so the first cell a participant runs
  # does not pay for loading these into a fresh session's cache.
  julia --project=. -e 'using EpiModelingCourse, CairoMakie, OrdinaryDiffEq, DataFrames, StableRNGs; println("Course environment ready on Julia ", VERSION)'

  # Quarto's Julia engine precompiles its worker the first time any document is
  # executed, which otherwise costs a participant roughly a minute and prints
  # into the output of whichever cell they happened to run first. The CI
  # workflow absorbs this the same way, with a throwaway render.
  cat > warmup.qmd <<'QMD'
---
title: Warmup
---

```{julia}
using CairoMakie, DataFrames
DataFrame(warm = 1:2)
lines(1:10, collect(1.0:10.0))
```
QMD
  quarto render warmup.qmd >/dev/null 2>&1 || echo "warm-up render failed; first worksheet render will be slower"
  rm -rf warmup.qmd warmup.html warmup_files _freeze/warmup
fi
