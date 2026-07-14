# Modeling Epidemics with Julia

An executable course that develops epidemic models from assumptions and equations through simulation, inference, and decision-making in Julia.

The course contains 3 optional prerequisites, 12 core chapters, and 9 advanced appendices. Each teaching unit has a Quarto book page, Reveal.js slides, a Pluto notebook, and two exercises with solutions.

## Local setup

Requirements:

- Julia 1.12
- Quarto 1.9 or newer

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=environments/ch08 -e 'using Pkg; Pkg.instantiate()'
```

## Build and validation

```bash
# Render the book and all slide decks.
quarto render

# Export every Pluto notebook to static HTML.
julia --project=. scripts/export_notebooks.jl

# Check the 24-unit artifact contract and 48 exercise/solution pairs.
julia --project=. scripts/check_structure.jl

# Verify the capstone dataset and checksum.
julia --project=. scripts/check_data.jl

# Reject obsolete package imports and third-party monkey-patches.
julia --project=. scripts/check_deprecated_apis.jl

# Check figure captions, slide notes, notebook interactivity, solutions, and citations.
julia --project=. scripts/check_editorial_contract.jl

# Run shared model and data tests.
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Self-contained PDF

The PDF build renders the complete book, renders all 24 decks with Beamer, appends the decks after the book, embeds fonts, and compresses the result:

```bash
julia --project=. scripts/build_combined_pdf.jl
```

The output is `dist/modeling-epidemics-with-julia-complete.pdf`, with a matching SHA-256 file. The build requires Quarto/TinyTeX, `pdfunite`, and Ghostscript. To recombine already-rendered intermediate PDFs without re-executing the course:

```bash
julia --project=. scripts/build_combined_pdf.jl --skip-render
```

Open the interactive notebooks with:

```bash
julia --project=. -e 'using Pluto; Pluto.run()'
```

## Publication

`.github/workflows/publish.yml` renders the committed freeze results with the `publish` Quarto profile and deploys `_site/` through GitHub Pages. After adding a GitHub remote, configure the repository's Pages source as **GitHub Actions** and push `main`.

The sibling `../sir-julia/` repository is outside this course repository and is used only as local reference material.
