# ======================
# Domains
# ======================

# The following domains are used throughout the tests in Tables 3-5 in
# [1] Althoff, M., Grebenyuk, D., & Kochdumper, N. (2018). Implementation of Taylor models in CORA 2018.
#     In Proc. of the 5th International Workshop on Applied Verification for Continuous and Hybrid Systems.
a = interval(-4.5, -0.3)
b = interval(0.4, 0.9)
c = interval(3.8, 7.8)
d = interval(8.0, 10.0)
e = interval(-10.0, 8.0)
f = interval(1.0, 2.0)

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
    return 100 * interval(rel_low, rel_high)
end

# ==========================
# Methods to compute bounds
# ==========================

# taylor model in one variable
function bounds_TM(func::Function, dom::Interval, ord::Int)
    x0 = interval(mid(dom))
    x = TaylorModel1(ord, x0, dom)
    return evaluate(func(x), dom - x0)
end

# taylor model in N variables
function bounds_TM(func::Function, dom::AbstractVector{<:Interval}, ord)
    x0 = mid.(dom)
    set_variables(Float64, "x", order=2ord, numvars=N)
    x = [TaylorModelN(i, ord, x0, dom) for i in eachindex(dom)]
    return evaluate(func(x...), dom - x0)
end

# normalized taylor model in one variable
function bounds_TM_NORM(func::Function, dom::Interval, ord::Int)
    x0 = interval(mid(dom))
    x = TaylorModel1(ord, x0, dom)
    xnorm = normalize_taylor(x.pol, dom - x0, true)
    xnormTM = TaylorModel1(xnorm, 0..0, 0..0, -1..1)
    return evaluate(func(xnormTM), -1..1)
end

# normalized taylor model in N variables
function bounds_TM_NORM(func::Function, dom::AbstractVector{<:Interval}, ord::Int)
    x0 = mid.(dom)
    N = length(dom)
    set_variables(Float64, "x", order=2ord, numvars=N)

    zeroBox = fill(0..0, N)
    symBox = fill(-1..1, N)

    x = [TaylorModelN(i, ord, x0, dom) for i in eachindex(dom)]
    xnorm = [normalize_taylor(xi.pol, dom - x0, true) for xi in x]
    xnormTM = [TaylorModelN(xi_norm, 0..0, zeroBox, symBox) for xi_norm in xnorm]
    return evaluate(func(xnormTM...), symBox)
end

# interval arithmetic substitution
function bounds_IA(f, dom)
    return f(dom...)
end
