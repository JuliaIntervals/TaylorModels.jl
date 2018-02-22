"""
integrate a TaylorNModel with respect to the variable number `which`.
Optionally adds `x0` to the result.
"""
function integrate(f::TaylorNModel, which=1, x0=0)

    p = integrate(f.p, which)  # not necessary if an already complete Taylor series, in which case p2 == f.p

    n = degree(f.p)
    high_order_term = f.p[n]  # a HomogeneousPolynomial

    Δ = ( bound(high_order_term, f.x0, f.I) + f.Δ ) * diam(f.I[which])

    g = TaylorNModel(f.n, f.x0, f.I, p, Δ)
    g.p[0] += x0  # constant term

    return g

end


bound(f::HomogeneousPolynomial, x0, I) = f( [(I - x0)...] )
# applies the hom poly to the bounds
