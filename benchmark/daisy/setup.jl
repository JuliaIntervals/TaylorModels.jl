# ======================
# Domains
# ======================
# The following domains are used throughout the tests in Tables 3-5 in
# [1] Althoff, M., Grebenyuk, D., & Kochdumper, N. (2018). Implementation of Taylor models in CORA 2018.
#     In Proc. of the 5th International Workshop on Applied Verification for Continuous and Hybrid Systems.
a = Interval(-4.5, -0.3)
b = Interval(0.4, 0.9)
c = Interval(3.8, 7.8)
d = Interval(8.0, 10.0)
e = Interval(-10.0, 8.0)
f = Interval(1.0, 2.0)

# ======================
# Normalization
# ======================
symmetric_domain(dom::Interval) = Interval(-1, 1)
symmetric_domain(dom::IntervalBox) = IntervalBox(-1..1, length(dom))

# Helper function to evaluate after normalization
function normalize_and_evaluate(p, dom)
    dom_sym = symmetric_domain(dom)
    evaluate(normalize_taylor(p, dom, true), dom_sym)
end

# ======================
# Relative precision
# ======================
# Dictionary to hold the vector of relative precision intervals for each benchmark
RELPREC = Dict{String, Any}()

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
