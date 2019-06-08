
## This file defines a benchmark suite with the tools provided by `PkgBenchmark` and `BenchmarkTools`.

To run the benchmarks, execute:
```julia
using PkgBenchmark
results = benchmarkpkg("LazySets")

#To compare current version to another tagged version, commit or branch:

julia> results = judge("LazySets", <tagged-version-or-branch>)

#To export the benchmark results to a Markdown file:

julia> export_markdown("results.md", results)


#To export the benchmark results to a JSON file:
julia> writeresults("results.json", results)
```
