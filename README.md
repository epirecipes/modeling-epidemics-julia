# Modeling Epidemics with Julia

An executable course that develops epidemic models from their assumptions and
equations through simulation, inference, and decision-making in Julia. Every
model is built up from what it assumes about a population, written as
equations, implemented in Julia, and then checked against something known.

## Start here

**Attending the workshop?** Work through [Setting up](https://epirecipes.github.io/modeling-epidemics-julia/setup.html) before Day 1.
It takes about half an hour, most of which is waiting for downloads, and doing
it in advance means the first session starts with you writing model code rather
than installing software. Then see [the schedule](https://epirecipes.github.io/modeling-epidemics-julia/workshop.html) for what is
taught when, and [the session resources](https://epirecipes.github.io/modeling-epidemics-julia/resources.html) for the chapter,
slides and worksheet of each taught session.

**Working through the course on your own?** Start with the same
[Setting up](https://epirecipes.github.io/modeling-epidemics-julia/setup.html) page, then read from [the introduction](https://epirecipes.github.io/modeling-epidemics-julia/)
onward. The workshop pages are safe to ignore: every chapter links its own
slides and notebook at the end, and each has two exercises with worked
solutions.

## What's inside

Twenty-four teaching units: 3 optional prerequisites covering Julia,
probability and calculus, 12 core chapters running from the SIR model to a
complete data-to-decision case study, and 9 advanced appendices extending the
same ideas to delays, scientific machine learning, probabilistic numerics,
control, and model composition.

Each unit comes as four things:

- a **chapter** that derives the model and implements it;
- a **slide deck** summarising it;
- a **Pluto notebook** with sliders, for changing a parameter and watching the
  figures move;
- **two exercises** with full solutions.

Ten of the twenty-four units are taught live over three days, and each of those
also has a **worksheet**: a partly written Quarto file you fill in during the
session.

## Getting the materials

Full instructions, including what to do when a step fails, are in
[Setting up](https://epirecipes.github.io/modeling-epidemics-julia/setup.html). In short:

```bash
git clone https://github.com/epirecipes/modeling-epidemics-julia.git
cd modeling-epidemics-julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then open that folder in VSCode with *File* → *Open Folder*, opening the
folder itself rather than an individual file. If you do not use git, download
the repository as a ZIP and unzip it somewhere you can find again.

The install step downloads every package the course uses, at the versions it
was tested with. Expect ten to twenty minutes, and do it before you travel
rather than on workshop wifi.

## Where things are

| Folder | What is in it |
|---|---|
| `prerequisites/`, `chapters/`, `appendices/` | the book pages, in reading order |
| `worksheets/` | the fill-in files for the taught sessions |
| `slides/` | the slide deck for every unit |
| `notebooks/` | the Pluto notebooks |
| `notebook-exports/` | static previews of those notebooks, no install needed |
| `data/` | the capstone dataset, with its checksum and provenance |

## Contributing

Corrections and improvements are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md)
for how to build, render and validate the course.

## License and citation

This book is developed by the epirecipes organization and builds upon material
from the [`epirecipes/sir-julia`](https://github.com/epirecipes/sir-julia)
repository.

The book, its source files, and accompanying code are distributed under the
MIT License. Copyright (c) 2020–2026 epirecipes.

### How to cite

```bibtex
@misc{epirecipes_modeling_epidemics_julia,
  author       = {Montes-Olivas, Sandra and
                  Frost, Simon D. W.},
  title        = {Modeling Epidemics with Julia: From Methods to Decisions},
  year         = {2026},
  url          = {https://github.com/epirecipes/modeling-epidemics-julia},
  note         = {MIT-licensed course materials}
}
```
