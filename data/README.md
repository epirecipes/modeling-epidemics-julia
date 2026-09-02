# Capstone data

`eyam_plague_1666.csv` contains eight SIR state observations for the second wave of the Eyam plague, beginning on 18 June 1666. Every row sums to the modeled population of 261.

The values are historical epidemiological facts packaged as the `Eyam` dataset in MultiBD 1.0.2, distributed under the Apache License 2.0. Cite both the package source and Raggett's original analysis:

- MultiBD: <https://github.com/msuchard/MultiBD>
- Raggett, G. F. (1982), "A stochastic model of the Eyam plague," *Journal of Applied Statistics* 9, 212–226. <https://doi.org/10.1080/02664768200000021>

The data are intentionally committed to the repository. Course rendering must not fetch live epidemiological data from the network.

See `eyam_plague_1666.toml` for field definitions, provenance, and modeling caveats.

# Language benchmark data

`language_benchmarks_2022.csv` and `language_benchmarks_2016.csv` hold the JuliaLang micro-benchmark timings used by the two motivating figures in the `p01-julia-essentials` slide deck. Columns are `language`, `benchmark`, `time` (seconds), with `lines` and `kb` for the source length of each language's implementation in the 2016 file. Neither file has a `time` entry for every language and benchmark pair, so the slide code keeps only the benchmarks C itself ran and normalises against those.

- 2022 timings: <https://github.com/JuliaLang/Microbenchmarks>, `bin/benchmarks.csv` at commit `a963d284b09d04b3e0374f6dd46ec4b039ed5569`. A header row was added; values are unchanged.
- 2016 timings and source lengths: retrieved through the Internet Archive from <https://web.archive.org/web/20230601010916/https://www.oceanographerschoice.com/>. Column names were lowercased; values are unchanged.

Source lengths stopped being published alongside the timings after 2016, which is why the two figures use different years.
