## How to run the benchmarks

The files in this folder define a benchmark suite with the tools provided by [PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl) and [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl).

To run the benchmarks, execute:

```julia
julia> using PkgBenchmark

julia> results = benchmarkpkg("TaylorModels")
```

## How to compare benchmarks

To compare current version to another tagged version, commit or branch:

```julia
julia> results = judge("TaylorModels", <tagged-version-or-branch>)
```

## Exporting results

To export the benchmark results to a Markdown file:

```julia
julia> export_markdown("results.md", results)
```

To export the benchmark results to a JSON file:

```julia
julia> writeresults("results.json", results)
```

## Daisy benchmarks

The dictionary `RELPREC` contains the relative precision for each test case.
You can extract the results and store them in a `DataFrame` as follows:

```julia
using DataFrames

include(".../TaylorModels/benchmark/benchmarks.jl") # fix to your case

relprec_uni = DataFrame(Test = String[], Interval_Arithm = Interval[],
                   TM_order_2 = Interval[], TM_order_5 = Interval[], TM_order_10 = Interval[],
                   norm_TM_order_2 = Interval[], norm_TM_order_5 = Interval[], norm_TM_order_10 = Interval[])

relprec_multi = DataFrame(Test = String[], Interval_Arithm = Interval[],
                     TM_order_2 = Interval[], TM_order_5 = Interval[], TM_order_10 = Interval[],
                     norm_TM_order_2 = Interval[], norm_TM_order_5 = Interval[], norm_TM_order_10 = Interval[]);

for name in keys(RELPREC)
    data = (name, RELPREC[name]["order 2"]["Interval Arithmetic subs."],
                             RELPREC[name]["order 2"]["Taylor Model subs."],
                             RELPREC[name]["order 5"]["Taylor Model subs."],
                             RELPREC[name]["order 10"]["Taylor Model subs."],
                             RELPREC[name]["order 2"]["Normalized Taylor Model subs."],
                             RELPREC[name]["order 5"]["Normalized Taylor Model subs."],
                             RELPREC[name]["order 10"]["Normalized Taylor Model subs."])
    DAISY[name]["numvars"] == 1 ? push!(relprec_uni, data) : push!(relprec_multi, data)
end
```
The results are now stored in the `relprec_uni` (resp. `relprec_multi` data frames). They can be displayed in a Jupyter notebook or in the REPL doing

```julia
using PrettyTables

pretty_table(relprec_uni)

pretty_table(relprec_multi)
```

The tables with runtimes can be generated similarly.

```julia
using TaylorModels, DataFrames

RES = results.benchmarkgroup["Daisy"]

runtime_uni = DataFrame(Test = String[], Interval_Arithm = Float64[],
                   TM_order_2 = Float64[], TM_order_5 = Float64[], TM_order_10 = Float64[],
                   norm_TM_order_2 = Float64[], norm_TM_order_5 = Float64[], norm_TM_order_10 = Float64[])

runtime_multi = DataFrame(Test = String[], Interval_Arithm = Float64[],
                     TM_order_2 = Float64[], TM_order_5 = Float64[], TM_order_10 = Float64[],
                     norm_TM_order_2 = Float64[], norm_TM_order_5 = Float64[], norm_TM_order_10 = Float64[]);

for name in keys(RELPREC)
    data = (name, time(RES[name]["order 2"]["Interval Arithmetic subs."])/1e3,
                  time(RES[name]["order 2"]["Taylor Model subs."])/1e3,
                  time(RES[name]["order 5"]["Taylor Model subs."])/1e3,
                  time(RES[name]["order 10"]["Taylor Model subs."])/1e3,
                  time(RES[name]["order 2"]["Normalized Taylor Model subs."])/1e3,
                  time(RES[name]["order 5"]["Normalized Taylor Model subs."])/1e3,
                  time(RES[name]["order 10"]["Normalized Taylor Model subs."])/1e3)
    DAISY[name]["numvars"] == 1 ? push!(runtime_uni, data) : push!(runtime_multi, data)
end
```
