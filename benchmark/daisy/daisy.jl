# ======================
# Initialization
# ======================

# The benchmarks in this file are taken from
# [Project Daisy](https://github.com/malyzajko/daisy/blob/master/testcases/).
SUITE["Daisy"] = BenchmarkGroup()
RESULTS = SUITE["Daisy"]

# Dictionaries to store data and reference values
DAISY = Dict{String, Any}()

include("setup.jl")
include("univariate.jl")
include("multivariate.jl")

# ===============================
# Run benchmarks
# ===============================
for name in keys(DAISY)
    RESULTS[name] = BenchmarkGroup()
    RELPREC[name] = Dict()
    func = DAISY[name]["f"]
    dom = DAISY[name]["dom"]
    ref = DAISY[name]["ref"]
    for ord in [2, 5, 10]
        RESULTS[name]["order $ord"] = BenchmarkGroup()
        RELPREC[name]["order $ord"] = Dict()

        RESULTS[name]["order $ord"]["Taylor Model subs."] = @benchmarkable bounds_TM($func, $dom, $ord)
        RELPREC[name]["order $ord"]["Taylor Model subs."] = relative_precision(bounds_TM(func, dom, ord), ref)

        RESULTS[name]["order $ord"]["Normalized Taylor Model subs."] = @benchmarkable bounds_TM_NORM($func, $dom, $ord)
        RELPREC[name]["order $ord"]["Normalized Taylor Model subs."] = relative_precision(bounds_TM_NORM(func, dom, ord), ref)

        RESULTS[name]["order $ord"]["Interval Arithmetic subs."] = @benchmarkable bounds_IA($func, $dom)
        RELPREC[name]["order $ord"]["Interval Arithmetic subs."] = relative_precision(bounds_IA(func, dom), ref)
    end
end
