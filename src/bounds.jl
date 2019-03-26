# bounds.jl

"""
   bound_absrem(f::Function, polf::Taylor1, polfI::Taylor1, x0::Interval, I::Interval)

Bound the absolute remainder of the polynomial approximation of `f` given
by the Taylor polynomial `polf` around `x0` on the interval `I`. It requires
the interval extension `polfI` of the polynomial that approximates `f` for
the whole interval `I`, in order to compute the Lagrange remainder.

If `polfI[end]` has a definite sign, then it is monotonic in the intervals
[I.lo, x0] and [x0.hi, I.hi], which is exploited; otherwise, it is used
to compute the Lagrange remainder.

"""
function bound_absrem(f::Function, polf::Taylor1, polfI::Taylor1,
        x0, I::Interval)

    _order = get_order(polf) + 1
    fTIend = polfI[_order]
    if (sup(fTIend) < 0 || inf(fTIend) > 0)
        # Absolute remainder is monotonic
        a = interval(I.lo)
        b = interval(I.hi)
        Δlo = f(a) - polf(a-x0)
        # Δlo = f(a) - bound_taylor1(polf, a-x0)
        Δhi = f(b) - polf(b-x0)
        # Δhi = f(b) - bound_taylor1(polf, b-x0)
        Δx0 = f(x0) - polf[0]
        Δ = hull(hull(Δlo, Δx0), Δhi)
    else
        # Lagrange bound
        Δ = fTIend * (I-x0)^_order
    end
    return Δ
end


"""
   bound_relrem(f::Function, polf::Taylor1, polfI::Taylor1, x0::Interval, I::Interval)

Bound the relative remainder of the polynomial approximation of `f` given
by the Taylor polynomial `polf` around `x0` on the interval `I`. It requires
an the interval extension `polfI` of a polynomial that approximates `f` for
the whole interval `I`, in order to compute the Lagrange remainder.

If `polfI[end]` has a definite sign, then it is monotonic in the interval `I`,
which is exploited; otherwise, the last coefficients bounds the relative
remainder.

"""
function bound_relrem(f::Function, polf::Taylor1, polfI::Taylor1,
        x0::Interval, I::Interval)

    _order = get_order(polf) + 1
    fTIend = polfI[_order+1]
    a = interval(I.lo)
    b = interval(I.hi)
    if (sup(fTIend) < 0 || inf(fTIend) > 0) && isempty(a ∩ x0) && isempty(b ∩ x0)
        # Error is monotonic
        denom_lo = (a-x0)^_order
        Δlo = f(a) - polf(a-x0)
        # Δlo = f(a) - bound_taylor1(polf, a-x0)
        Δlo = Δlo / denom_lo
        denom_hi = (b-x0)^_order
        Δhi = f(b) - polf(b-x0)
        # Δhi = f(b) - bound_taylor1(polf, b-x0)
        Δhi = Δhi / denom_hi
        Δ = hull(Δlo, Δhi)
    else
        # Lagrange coefficient
        Δ = polfI[_order]
    end
    return Δ
end


"""
    bound_taylor1(fT::Taylor1, I::Interval)

Compute a *tight* polynomial bound for the Taylor polynomial `fT`
in the interval `I`.

Note: Algorithm 2.1.1 corresponds to `evaluate(fT, I)` or simply `fT(I).
This function uses the roots of the derivative of `ft` to obtain a
tighter bound.

"""
function bound_taylor1(fT::Taylor1, I::Interval)

    # Check if the fT is monotonous (the derivative has a given sign)
    fTd  = TaylorSeries.derivative(fT)
    range_deriv = fTd(I)
    (sup(range_deriv) ≤ 0 || inf(range_deriv) ≥ 0) && return bound_taylor1(fT, fTd, I)

    # Compute roots of the derivative using the second derivative
    # Fix some sort of relative tolerance for Newton root search
    fTd2 = TaylorSeries.derivative(fTd)
    rootsder = roots(fTd, fTd2, I, Newton, 1.0e-5*mag(I))

    # Bound the range of fT using the roots and end points
    num_roots = length(rootsder)
    num_roots == 0 && return fT(I)
    rangepoly = hull( fT(interval(I.lo)), fT(interval(I.hi)) )
    for ind in 1:num_roots
        rangepoly = hull(rangepoly, fT(rootsder[ind].interval))
    end

    return rangepoly
end


"""
    bound_taylor1(fT::Taylor1, fTd::Taylor1, I::Interval)

Compute a *tight* polynomial bound for the Taylor polynomial `fT`
in the interval `I`, considering whether its derivative `ftd` has
a definite sign.

"""
function bound_taylor1(fT::Taylor1{T}, fTd::Taylor1{T}, I::Interval{T}) where {T}
    #
    if inf(fTd(I)) ≥ 0
        return @interval(fT(I.lo), fT(I.hi))
    elseif sup(fTd(I)) ≤ 0
        return @interval(fT(I.hi), fT(I.lo))
    end
    return fT(I)
end
function bound_taylor1(fT::Taylor1{Interval{T}}, fTd::Taylor1{Interval{T}},
        I::Interval{S}) where {T,S}
    #
    if inf(fTd(I)) ≥ 0
        return hull(fT(I.lo), fT(I.hi))
    elseif sup(fTd(I)) ≤ 0
        return hull(fT(I.hi), fT(I.lo))
    end
    return fT(I)
end
