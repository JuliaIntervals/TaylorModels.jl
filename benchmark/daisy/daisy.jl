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
        x0 = Interval(mid(dom))
        x = TaylorModel1(ord, x0, dom)
        p = DAISY_1D[name]["f"](x)
        xnorm = normalize_taylor(x.pol, dom - x0, true)
        xnormTM = TaylorModel1(xnorm, 0..0, 0..0, -1..1)
        q = DAISY_1D[name]["f"](xnormTM)
        RESULTS[name]["order $ord"] = BenchmarkGroup()
        RELPREC[name]["order $ord"] = Dict()

        RESULTS[name]["order $ord"]["evaluate"] = @benchmarkable evaluate($p, $(dom-x0))
        approx = evaluate(p, dom-x0)
        RELPREC[name]["order $ord"]["evaluate"] = relative_precision(approx, ref)

        RESULTS[name]["order $ord"]["normalize and evaluate"] = @benchmarkable evaluate($q, $(-1..1))
        approx = evaluate(q, -1..1)
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

        set_variables(Float64, "vars", order=2ord, numvars=numvars)
        vars = [TaylorModelN(i, ord, IntervalBox(mid(dom)), dom) for i in 1:numvars]
        x0 = mid(dom)
        zeroBox = IntervalBox(0..0, numvars)
        sym_dom = IntervalBox(-1..1, numvars)
        xnorm = [normalize_taylor(x.pol, dom - x0, true) for x in vars]
        xnormTM = [TaylorModelN(x_norm, 0..0, zeroBox, sym_dom) for x_norm in xnorm]
        f = DAISY_ND[name]["f"]
        p = f(vars...)
        q = f(xnormTM...)
        RESULTS[name]["order $ord"] = BenchmarkGroup()
        RELPREC[name]["order $ord"] = Dict()

        RESULTS[name]["order $ord"]["evaluate"] = @benchmarkable evaluate($p, $(dom-x0)) setup=(
            vars = set_variables(Float64, $varnames, order=$ord, numvars=$numvars); p = $f(vars...))
        approx = evaluate(p, dom-x0)
        RELPREC[name]["order $ord"]["evaluate"] = relative_precision(approx, ref)

        RESULTS[name]["order $ord"]["normalize and evaluate"] = @benchmarkable evaluate($q, $(sym_dom)) setup=(
            xnormTM = set_variables(Float64, $varnames, order=$ord, numvars=$numvars); q =$f(xnormTM...))
        approx = evaluate(q, sym_dom)
        RELPREC[name]["order $ord"]["normalize and evaluate"] = relative_precision(approx, ref)
    end
end
