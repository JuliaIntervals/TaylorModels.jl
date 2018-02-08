# bounds.jl

"""
   boundarem(f::Function, polf::Taylor1, fTIend::Interval, x0::Interval, ii::Interval)

Bound the absoulte remainder of the polynomial approximation of `f` given
by the Taylor polynomial `polf` around `x0` on the interval `ii`.
It requires an interval bound `fTIend` for the coefficient of Laplace bound. If
`fTIend` has a definite sign, then it is monotonic in the intervals [ii.lo, x0]
and [x0.hi, ii.hi], which is exploited; otherwise, it is used as the returned
value.

"""
function boundarem(f::Function, polf::Taylor1, polfI::Taylor1,
        x0::Interval, ii::Interval)

    fTIend = polfI[end]
    if (sup(fTIend) ≤ 0 || inf(fTIend) ≥ 0)
        # Absolute remainder is monotonic
        a = interval(ii.lo)
        b = interval(ii.hi)
        Δlo = f(a) - polf(a-x0)
        Δhi = f(b) - polf(b-x0)
        Δx0 = f(x0) - polf[0]
        Δ = interval( min(inf(Δlo), inf(Δx0), inf(Δhi)),
                      max(sup(Δlo), sup(Δx0), sup(Δhi)) )
    else
        # Laplace bound
        Δ = fTIend * (ii-x0)^(polf.order+1)
    end
    return Δ
end


"""
   boundrrem(f::Function, polf::Taylor1, fTIend::Interval, x0::Interval, ii::Interval)

Bound the relative remainder of the polynomial approximation of `f` given
by the Taylor polynomial `polf` around `x0` on the interval `ii`.
It requires an interval bound `fTIend` for the coefficient of Laplace bound. If
`fTIend` has a definite sign, then it is monotonic in the interval `ii`,
which is exploited; otherwise, it is used as the returned value.

"""
function boundrrem(f::Function, polf::Taylor1, polfI::Taylor1,
        x0::Interval, ii::Interval)

    _order = polf.order + 1
    fTIend = polfI[end]
    a = interval(ii.lo)
    b = interval(ii.hi)
    if (sup(fTIend) ≤ 0 || inf(fTIend) ≥ 0) && isempty(a ∩ x0) && isempty(b ∩ x0)
        # Error is monotonic
        lodenom = (a-x0)^_order
        Δlo = f(a) - polf(a-x0)
        Δlo = Δlo / lodenom
        hidenom = (b-x0)^_order
        Δhi = f(b) - polf(b-x0)
        Δhi = Δhi / hidenom
        Δ = interval( min(inf(Δlo), inf(Δhi)), max(sup(Δlo), sup(Δhi)) )
    else
        # Laplace coefficient
        Δ = polfI[_order]
    end
    return Δ
end


"""
    bound_taylor1(fT::Taylor1, ii::Interval)

Compute a *tight* polynomial bound for the Taylor polynomial `fT`
for the interval `ii`.

Note: Algorithm 2.1.1 corresponds to `evaluate(fT, ii)` or simply `fT(ii).
This function uses the roots of the derivative of `ft` to obtain a
tighter bound.

"""
function bound_taylor1(fT::Taylor1, ii::Interval)

    # Check if the derivative is monotonous
    fTd  = TaylorSeries.derivative(fT)
    (sup(fTd(ii)) ≤ 0 || inf(fTd(ii)) ≥ 0) && return bound_taylor1(fT, fTd, ii)

    # Compute roots of the derivative using the second derivative
    fTd2 = TaylorSeries.derivative(fTd)
    rootsder = find_roots(x->fTd(x), x->fTd2(x), ii)

    num_roots = length(rootsder)
    # num_roots == 0 && return bound_taylor1(fT, fTd, ii)

    # Obtain the range exploiting monotonicity in the intervals between the roots
    vi = typeof(rootsder[end].interval)[]
    for ind in 1:num_roots+1
        # Define the interval of interest
        if ind == 1
            iit = interval(ii.lo, inf(rootsder[1].interval))
        elseif ind == num_roots+1
            iit = interval(sup(rootsder[num_roots].interval), ii.hi)
        else
            iit = interval(sup(rootsder[ind-1].interval), inf(rootsder[ind].interval))
        end

        # Exploit monotonicity
        if fTd(iit) ≥ 0.0           # fT is increasing
            push!(vi, @interval(fT(iit.lo), fT(iit.hi)))
        elseif fTd(iit) ≤ 0.0       # fT is decreasing
            push!(vi, @interval(fT(iit.hi), fT(iit.lo)))
        end
    end

    # Compute the hull of all elements in `vi`
    rangepoly = vi[end]
    for ind in eachindex(rootsder)
        rangepoly = hull(rangepoly, vi[ind])
    end

    return rangepoly
end


"""
    bound_taylor1(fT::Taylor1, fTd::Taylor1, ii::Interval)

Compute a *tight* polynomial bound for the Taylor polynomial `fT`
in the interval `ii`, considering whether its derivative `ftd` has
a definite sign.

"""
function bound_taylor1(fT::Taylor1{T}, fTd::Taylor1{T}, ii::Interval{T}) where {T}
    #
    if inf(fTd(ii)) ≥ 0
        return @interval(fT(ii.lo), fT(ii.hi))
    elseif sup(fTd(ii)) ≤ 0
        return @interval(fT(ii.hi), fT(ii.lo))
    end
    return fT(ii)
end
function bound_taylor1(fT::Taylor1{Interval{T}}, fTd::Taylor1{Interval{T}},
        ii::Interval{S}) where {T,S}
    #
    if inf(fTd(ii)) ≥ 0
        return hull(fT(ii.lo), fT(ii.hi))
    elseif sup(fTd(ii)) ≤ 0
        return hull(fT(ii.hi), fT(ii.lo))
    end
    return fT(ii)
end
