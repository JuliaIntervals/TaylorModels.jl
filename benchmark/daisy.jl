# ==========
# Daisy
# ==========

# The benchmarks in this file are taken from
# [Project Daisy](https://github.com/malyzajko/daisy/blob/master/testcases/).
SUITE["Daisy"] = BenchmarkGroup()
DAISY = SUITE["Daisy"]

# The following domains are used throughout the tests in Tables 3-5 in
# [1] Althoff, M., Grebenyuk, D., & Kochdumper, N. (2018). Implementation of Taylor models in CORA 2018.
#     In Proc. of the 5th International Workshop on Applied Verification for Continuous and Hybrid Systems.
a = Interval(-4.5, -0.3)
b = Interval(0.4, 0.9)
c = Interval(3.8, 7.8)
d = Interval(8.0, 10.0)
e = Interval(-10.0, 8.0)
f = Interval(1.0, 2.0)

# This function measures the relative precision of the result in a more informative way than
# taking the scalar overestimation because it evaluates the precision of the lower and the
# upper range bounds separately, see Eq. (20) in [1].
function relative_precision(x, x_ref)
    x_low, x_high = inf(x), sup(x)
    x_ref_low, x_ref_high = inf(x_ref), sup(x_ref)
    rel_low = -(x_low - x_ref_low) / (x_ref_high - x_ref_low)
    rel_high = (x_high - x_ref_high) / (x_ref_high - x_ref_low)
    return 100 * Interval(rel_low, rel_high)
end

# Helper function to evaluate after normalization
symmetric_domain(dom::Interval) = Interval(-1, 1)
symmetric_domain(dom::IntervalBox) = IntervalBox(-1..1, length(dom))

function normalize_and_evaluate(p, dom)
    dom_sym = symmetric_domain(dom)
    evaluate(normalize_taylor(p, dom, true), dom_sym)
end

# Dictionary to hold the vector of relative precision intervals for each benchmark
RELPREC = Dict{String, Any}()

# ======
# sine
# ======

DAISY["sin"] = BenchmarkGroup()
RELPREC["sin"] = Dict()
ref = Interval(-0.9998434996853951, 2.724232859130124) # IntervalOptimisation
dom = a

for m in [2, 5, 10]
    x = Taylor1(m)
    p = sin(x)

    DAISY["sin"]["order $m"] = BenchmarkGroup()
    RELPREC["sin"]["order $m"] = Dict()

    DAISY["sin"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["sin"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["sin"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["sin"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# bspline0
# ==========

DAISY["bspline0"] = BenchmarkGroup()
RELPREC["bspline0"] = Dict()
ref = Interval(0.36616666317087393, 27.729165285894855) # MOSEK deg 5
dom = a

for m in [2, 5, 10]
    x = Taylor1(m)
    p = (1 - x) * (1 - x) * (1 - x) / 6.0

    DAISY["bspline0"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline0"]["order $m"] = Dict()

    DAISY["bspline0"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline0"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline0"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline0"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# bspline1
# ==========

DAISY["bspline1"] = BenchmarkGroup()
RELPREC["bspline1"] = Dict()
ref = Interval(-65.14582571514087, 0.5631666653631329) # MOSEK deg 5
dom = a

for m in [2, 5, 10]
    x = Taylor1(m)
    p = (3*x*x*x - 6*x*x + 4) / 6.0

    DAISY["bspline1"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline1"]["order $m"] = Dict()

    DAISY["bspline1"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline1"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline1"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline1"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# bspline2
# ==========

DAISY["bspline2"] = BenchmarkGroup()
RELPREC["bspline2"] = Dict()
ref = Interval(0.07407407744904644, 53.60416654726812) # MOSEK deg 5
dom = a

for m in [2, 5, 10]
    x = Taylor1(m)
    p = (-3*x*x*x  + 3*x*x + 3*x + 1) / 6.0

    DAISY["bspline2"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline2"]["order $m"] = Dict()

    DAISY["bspline2"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline2"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# bspline3
# ==========

DAISY["bspline3"] = BenchmarkGroup()
RELPREC["bspline3"] = Dict()
ref = Interval(0.004500048981347225, 15.18749986989248)  # MOSEK deg 5
dom = a

for m in [2, 5, 10]
    x = Taylor1(m)
    p = -x*x*x / 6.0

    DAISY["bspline3"]["order $m"] = BenchmarkGroup()
    RELPREC["bspline3"]["order $m"] = Dict()

    DAISY["bspline3"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["bspline3"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["bspline3"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["bspline3"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# himmilbeau
# ==========

DAISY["himmilbeau"] = BenchmarkGroup()
RELPREC["himmilbeau"] = Dict()
ref = Interval(85.46830677734748, 221.7338939301446) # MOSEK deg 6
#ref = Interval(85.46812223067747, 221.73401457671844) # SDPA deg 4
dom = a × b

for m in [2, 5, 10]
    x1, x2 = set_variables(Float64, "x1 x2", order=m, numvars=2)
    p = (x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7)

    DAISY["himmilbeau"]["order $m"] = BenchmarkGroup()
    RELPREC["himmilbeau"]["order $m"] = Dict()

    DAISY["himmilbeau"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["himmilbeau"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["himmilbeau"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["himmilbeau"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# kepler1
# ==========

DAISY["kepler1"] = BenchmarkGroup()
RELPREC["kepler1"] = Dict()
ref = Interval(-5.255935494810441, 7.321362422825775)  # MOSEK deg 6
dom = a×b×c×d×e×f

for m in [2, 5, 10]
    x1, x2, x3, x4, x5, x6 = set_variables(Float64, "x1 x2 x3 x4 x5 x6", order=m, numvars=6)
    p = x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6 + x1 * (-x1 + x2 + x3 - x4 + x5 + x6)

    DAISY["kepler1"]["order $m"] = BenchmarkGroup()
    RELPREC["kepler1"]["order $m"] = Dict()

    DAISY["kepler1"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["kepler1"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["kepler1"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["kepler1"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# kepler2
# ==========

DAISY["kepler2"] = BenchmarkGroup()
RELPREC["kepler2"] = Dict()
dom = a×b×c×d
ref = Interval(-195.36974909125482, 78.3669520644375)  # MOSEK deg 6

for m in [2, 5, 10]
    x1, x2, x3, x4 = set_variables(Float64,"x1 x2 x3 x4", order=m, numvars=4)
    p = x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4) +
        x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4

    DAISY["kepler2"]["order $m"] = BenchmarkGroup()
    RELPREC["kepler2"]["order $m"] = Dict()

    DAISY["kepler2"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["kepler2"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["kepler2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["kepler2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ==========
# kepler3
# ==========

DAISY["kepler3"] = BenchmarkGroup()
RELPREC["kepler3"] = Dict()
dom = a×b×c×d×e×f
ref = Interval(-309.8484155131222, 17.982082401462407) # MOSEK deg 6

for m in [2, 5, 10]
    x1, x2, x3, x4, x5, x6 = set_variables(Float64, "x1 x2 x3 x4 x5 x6", order=m, numvars=6)
    p =  x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6) +
         x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) +x3* x6 * (x1 + x2 - x3 + x4 + x5 - x6) -
         x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6

    DAISY["kepler3"]["order $m"] = BenchmarkGroup()
    RELPREC["kepler3"]["order $m"] = Dict()

    DAISY["kepler3"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["kepler3"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["kepler2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["kepler2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ============
# Rigidbody1
# ============

DAISY["Rigidbody1"] = BenchmarkGroup()
RELPREC["Rigidbody1"] = Dict()
dom = a×b×c
ref = Interval(-20.786552979420335, -0.540012836551535) # MOSEK deg 6

for m in [2, 5, 10]
    x1, x2, x3 = set_variables(Float64, "x1 x2 x3", order=m, numvars=3)
    p = -x1*x2 - 2*x2*x3 - x1 - x3

    DAISY["Rigidbody1"]["order $m"] = BenchmarkGroup()
    RELPREC["Rigidbody1"]["order $m"] = Dict()

    DAISY["Rigidbody1"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["Rigidbody1"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["Rigidbody1"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["Rigidbody1"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end

# ============
# Rigidbody2
# ============

DAISY["Rigidbody2"] = BenchmarkGroup()
RELPREC["Rigidbody2"] = Dict()
dom = a×b×c
ref = Interval(68.81138021006673, 359.98566570476504)  # MOSEK deg 6

for m in [2, 5, 10]
    x1, x2, x3 = set_variables(Float64, "x1 x2 x3", order=m, numvars=3)
    p = 2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2

    DAISY["Rigidbody2"]["order $m"] = BenchmarkGroup()
    RELPREC["Rigidbody2"]["order $m"] = Dict()

    DAISY["Rigidbody2"]["order $m"]["evaluate"] = @benchmarkable evaluate($p, $dom)
    approx = evaluate(p, dom)
    RELPREC["Rigidbody2"]["order $m"]["evaluate"] = relative_precision(approx, ref)

    DAISY["Rigidbody2"]["order $m"]["normalize and evaluate"] = @benchmarkable normalize_and_evaluate($p, $dom)
    approx = normalize_and_evaluate(p, dom)
    RELPREC["Rigidbody2"]["order $m"]["normalize and evaluate"] = relative_precision(approx, ref)
end
