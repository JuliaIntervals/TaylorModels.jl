# ======================
# Initialization
# ======================

# The benchmarks in this file are taken from
# [Project Daisy](https://github.com/malyzajko/daisy/blob/master/testcases/).
SUITE["Daisy"] = BenchmarkGroup()
RESULTS = SUITE["Daisy"]

# Dictionaries to store data and reference values
DAISY_1D = Dict{String, Any}()
DAISY_ND = Dict{String, Any}()

include("setup.jl")
include("univariate.jl")
include("multivariate.jl")

# ===============================
# Run one-dimensional benchmarks
# ===============================
for name in keys(DAISY_1D)
    RESULTS[name] = BenchmarkGroup()
    RELPREC[name] = Dict()
    dom = DAISY_1D[name]["dom"]
    ref = DAISY_1D[name]["ref"]
    for m in [2, 5, 10]
        x = Taylor1(m)
        p = DAISY_1D[name]["f"](x)

        RESULTS[name]["order $m"] = BenchmarkGroup()
        RELPREC[name]["order $m"] = Dict()

        RESULTS[name]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
        approx = evaluate(p, dom)
        RELPREC[name]["order $m"]["evaluate"] = relative_precision(approx, ref)

        RESULTS[name]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
        approx = normalize_and_evaluate(p, dom)
        RELPREC[name]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
    end
end

# ===============================
# Run multivariate benchmarks
# ===============================
for name in keys(DAISY_ND)
    RESULTS[name] = BenchmarkGroup()
    RELPREC[name] = Dict()
    dom = DAISY_ND[name]["dom"]
    ref = DAISY_ND[name]["ref"]
    for m in [2, 5, 10]
        varnames, numvars = DAISY_ND[name]["vars"], DAISY_ND[name]["numvars"]
        vars = set_variables(Float64, varnames, order=m, numvars=numvars)
        f = DAISY_ND[name]["f"]
        p = f(vars...)

        RESULTS[name]["order $m"] = BenchmarkGroup()
        RELPREC[name]["order $m"] = Dict()

        RESULTS[name]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom) setup=(
            vars = set_variables(Float64, $varnames, order=$m, numvars=$numvars); p = $f(vars...))
        approx = evaluate(p, dom)
        RELPREC[name]["order $m"]["evaluate"] = relative_precision(approx, ref)

        RESULTS[name]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom) setup=(
            vars = set_variables(Float64, $varnames, order=$m, numvars=$numvars); p =$f(vars...))
        approx = normalize_and_evaluate(p, dom)
        RELPREC[name]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
    end
end
