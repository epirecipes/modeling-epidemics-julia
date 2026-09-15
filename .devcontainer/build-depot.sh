#!/usr/bin/env bash
# Builds the Julia depot into the container image, so that a codespace starts by
# pulling an image that already has every package compiled rather than spending
# an hour compiling them. Run only from the Dockerfile, as the vscode user, so
# the depot lands in that user's home directory.
set -euo pipefail

export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-2}"
cd /opt/course

julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=environments/ch10 -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. -e 'using EpiModelingCourse, CairoMakie, OrdinaryDiffEq, DataFrames, StableRNGs; println("depot built on Julia ", VERSION)'

# Quarto's Julia engine compiles its own worker the first time any document is
# executed. Rendering a throwaway document here moves that cost into the image,
# the same trick the CI workflow uses to keep it out of a chapter's first cell.
cat > warmup.qmd <<'QMD'
---
title: Warmup
engine: julia
julia:
  exeflags: ["--project=/opt/course"]
---

```{julia}
using CairoMakie, DataFrames
DataFrame(warm = 1:2)
lines(1:10, collect(1.0:10.0))
```
QMD
quarto render warmup.qmd --to html
rm -rf warmup.qmd warmup.html warmup_files _freeze

julia --project=. -e 'using Pkg; Pkg.gc()' >/dev/null
