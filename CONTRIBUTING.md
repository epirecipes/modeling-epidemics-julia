# Contributing

This file is for people editing the course. If you are following the course,
see [Setting up](setup.qmd) instead.

## Requirements

- Julia 1.12
- Quarto 1.9 or newer

Two environments must be instantiated. The second serves chapter 8 alone,
which needs a package that conflicts with the rest of the course:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=environments/ch08 -e 'using Pkg; Pkg.instantiate()'
```

## The rule that catches people first

`_freeze/` and `slides/_freeze/` are **committed to the repository**, and
execution is set to `freeze: auto`. Any change that touches executed code, or
even the prose around it, must be rendered locally and the refreshed freeze
committed in the same commit as the source file:

```bash
quarto render chapters/06-agent-based-models.qmd
```

`scripts/check_freeze_drift.jl` compares the working tree against `HEAD`, so
after editing a `.qmd` it reports drift until the refreshed freeze is committed
alongside it. That is the check working, not a failure. Run it before
committing, not mid-edit.

It compares content rather than bytes: figures are embedded as base64 PNGs that
a Linux runner cannot reproduce byte for byte from a macOS render, so image
payloads are ignored and numbers are compared at `rtol=1e-3`.

## Every unit has four artifacts

Each teaching unit is a book page, a slide deck, a Pluto notebook and a pair of
solutions, and they must stay in lockstep. They are keyed by a shared
`<number>-<slug>` stem:

| Artifact | Path |
|---|---|
| Book page | `chapters/06-agent-based-models.qmd` |
| Slide deck | `slides/chapters/06-agent-based-models.qmd` |
| Pluto notebook | `notebooks/chapters/06-agent-based-models.jl` |
| Solutions | `appendices/solutions/<group>.qmd`, via `{#solution-ch06-ex1}` anchors |

[course-units.toml](course-units.toml) is the machine-readable map of the
course and is what every script reads first, so a change to the curriculum
structure starts there rather than in the files. Adding or renaming a unit
means touching all four artifacts plus `course-units.toml` and the chapter list
in `_quarto.yml`.

## Build

```bash
# Render the book and all slide decks.
quarto render

# Export every Pluto notebook to static HTML.
julia --project=. scripts/export_notebooks.jl

# Run shared model and data tests.
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Validation

CI runs all seven of these, and they are the fastest way to check your work.
Each exits nonzero with specific messages.

```bash
# The 24-unit artifact contract and 48 exercise/solution pairs.
julia --project=. scripts/check_structure.jl

# Figure captions, slide notes, notebook interactivity, solutions, citations.
julia --project=. scripts/check_editorial_contract.jl

# Obsolete package imports and third-party monkey-patches.
julia --project=. scripts/check_deprecated_apis.jl

# PDF-safe math delimiters and figure labels.
julia --project=. scripts/check_pdf_figures.jl

# The capstone dataset and its SHA-256.
julia --project=. scripts/check_data.jl

# Freeze and export size ceilings.
julia --project=. scripts/check_artifact_sizes.jl

# The committed freeze still matches a re-render.
julia --project=. scripts/check_freeze_drift.jl
```

## Self-contained PDF

The PDF build renders the complete book, renders all 24 decks with Beamer,
appends the decks after the book, embeds fonts, and compresses the result:

```bash
julia --project=. scripts/build_combined_pdf.jl
```

The output is `dist/modeling-epidemics-with-julia-complete.pdf`, with a
matching SHA-256 file. The build requires Quarto/TinyTeX, `pdfunite`, and
Ghostscript. To recombine already-rendered intermediate PDFs without
re-executing the course:

```bash
julia --project=. scripts/build_combined_pdf.jl --skip-render
```

## Publication

`.github/workflows/publish.yml` renders the committed freeze results with the
`publish` Quarto profile and deploys `_site/` through GitHub Pages. Configure
the repository's Pages source as **GitHub Actions** and push `main`.

## Reference material

The sibling `../sir-julia/` repository is outside this course repository and is
used only as local reference material. Its examples predate current APIs, so
verify against current package documentation before adapting anything from it.
