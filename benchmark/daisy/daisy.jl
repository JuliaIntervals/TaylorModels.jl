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
    for ord in [2, 5, 10]
        x = TaylorModel1(Taylor1(2*ord), 0..0, Interval(mid(dom)), dom)
        p = DAISY_1D[name]["f"](x)

        RESULTS[name]["order $ord"] = BenchmarkGroup()
        RELPREC[name]["order $ord"] = Dict()

        RESULTS[name]["order $ord"]["evaluate"] = @benchmarkable evaluate($p, $dom)
        approx = evaluate(p, dom)
        RELPREC[name]["order $ord"]["evaluate"] = relative_precision(approx, ref)

        RESULTS[name]["order $ord"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($(p.pol), $dom)
        approx = normalize_and_evaluate(p.pol, dom)
        RELPREC[name]["order $ord"]["normalize and evaluate"] = relative_precision(approx, ref)
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
    for ord in [2, 5, 10]
        varnames, numvars = DAISY_ND[name]["vars"], DAISY_ND[name]["numvars"]

        vars = set_variables(Float64, varnames, order=ord, numvars=numvars)
        f = DAISY_ND[name]["f"]
        p = f(vars...)

        RESULTS[name]["order $ord"] = BenchmarkGroup()
        RELPREC[name]["order $ord"] = Dict()

        RESULTS[name]["order $ord"]["evaluate"] = @benchmarkable evaluate($p, $dom) setup=(
            vars = set_variables(Float64, $varnames, order=$ord, numvars=$numvars); p = $f(vars...))
        approx = evaluate(p, dom)
        RELPREC[name]["order $ord"]["evaluate"] = relative_precision(approx, ref)

        RESULTS[name]["order $ord"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom) setup=(
            vars = set_variables(Float64, $varnames, order=$ord, numvars=$numvars); p =$f(vars...))
        approx = normalize_and_evaluate(p, dom)
        RELPREC[name]["order $ord"]["normalize and evaluate"] = relative_precision(approx, ref)
    end
end
