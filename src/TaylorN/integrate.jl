"""
integrate a TaylorNModel with respect to the variable number `which`.
Optionally adds `x0` to the result.
"""
function integrate(f::TaylorNModel, which=1, x0=0)

    p = integrate(f.p, which)  # not necessary if an already complete Taylor series, in which case p2 == f.p

    Δ = integral_bound(f, which)

    g = TaylorNModel(f.n, f.x0, f.I, p, Δ)

    g.p[0] += x0  # constant term

    # for k in 0:(x0.order)
    #     TaylorSeries.add!(g.p, g.p, x0, k)
    # end

    return g

end

doc"""
Bound the integral of a `TaylorNModel` `f` with respect to the variable `which`.
"""
function integral_bound(f::TaylorNModel, which)

    high_order_term = f.p[end]  # a HomogeneousPolynomial

    Δ = ( bound(high_order_term, f.x0, f.I) + f.Δ ) * diam(f.I[which])

    return Δ
end


bound(f::HomogeneousPolynomial, x0, I) = f( [(I - x0)...] )
# applies the hom poly to the bounds
