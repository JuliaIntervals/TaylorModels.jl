# bounds.jl

"""
   bound_remainder(::Type{TaylorModel1}f::Function, polf::Taylor1, polfI::Taylor1, x0::Interval, I::Interval)

Bound the absolute remainder of the polynomial approximation of `f` given
by the Taylor polynomial `polf` around `x0` on the interval `I`. It requires
the interval extension `polfI` of the polynomial that approximates `f` for
the whole interval `I`, in order to compute the Lagrange remainder.

If `polfI[end]` has a definite sign, then it is monotonic in the intervals
[I.lo, x0] and [x0.hi, I.hi], which is exploited; otherwise, it is used
to compute the Lagrange remainder. This corresponds to Prop 2.2.1 in Mioara
Joldes PhD thesis (pp 52).

"""
function bound_remainder(::Type{TaylorModel1}, f::Function, polf::Taylor1, polfI::Taylor1, x0, I::Interval)

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
   bound_remainder(::Type{RTaylorModel1}, f::Function, polf::Taylor1, polfI::Taylor1, x0::Interval, I::Interval)

Bound the relative remainder of the polynomial approximation of `f` given
by the Taylor polynomial `polf` around `x0` on the interval `I`. It requires
an the interval extension `polfI` of a polynomial that approximates `f` for
the whole interval `I`, in order to compute the Lagrange remainder.

If `polfI[end]` has a definite sign, then it is monotonic in the interval `I`,
which is exploited; otherwise, the last coefficients bounds the relative
remainder. This corresponds to Prop 2.3.7 in Mioara Joldes' PhD thesis (pp 67).

"""
function bound_remainder(::Type{RTaylorModel1}, f::Function, polf::Taylor1, polfI::Taylor1, x0, I::Interval)

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
    rootsder = roots(x->fTd(x), x->fTd2(x), I, Newton, 1.0e-5*mag(I))

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
        return Interval(fT(I.lo), fT(I.hi))
    elseif sup(fTd(I)) ≤ 0
        return Interval(fT(I.hi), fT(I.lo))
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

"""
    bound_taylor1(fT::TaylorModel1, I=domain(fT)::Interval)

Compute a *tight* polynomial bound for the Taylor model `fT`
in the interval `I`, considering whether its derivative `ftd` has
a definite sign.

"""
bound_taylor1(fT::TaylorModel1, I=domain(fT)::Interval) = bound_taylor1(polynomial(fT), I)

"""
    linear_dominated_bounder(fT::TaylorModel1, ϵ=1e-3::Float64, max_iter=5::Int)

Compute a tighter polynomial bound for the Taylor model `fT` by the linear
dominated bounder algorithm. The linear dominated algorithm is applied until
the bound of `fT` gets tighter than `ϵ` or the number of steps reachs `max_iter`.
The returned bound corresponds to the improved polynomial bound with the remainder
of the `TaylorModel1` included.
"""
function linear_dominated_bounder(fT::TaylorModel1{S, T}; ϵ=1e-3, max_iter=5) where {S, T}
    d = 1.
    Dn = fT.dom
    Dm = Dn
    Pm = Taylor1(copy(fT.pol.coeffs))
    bound = interval(0.)
    n_iter = 0
    while d > ϵ && n_iter < max_iter
        if n_iter == 0
            x0 = mid(Dn) - mid(fT.x0)
        else
            x0 = mid(Dn) - mid(Dm)
        end
        c = mid(Dn)
        update!(Pm, x0)
        linear = Taylor1(Pm.coeffs[1:2], Pm.order)
        non_linear = Pm - linear
        if S <: Interval
            Li = mid(linear[1])
        else
            Li = linear[1]
        end
        I1 = linear(Dn - c)
        Ih = non_linear(Dn - c)
        bound = I1.lo + Ih
        d = diam(bound)
        n_iter += 1
        if Li == 0
            break
        elseif Li > 0
            new_hi = min(Dn.lo + (d / abs(Li)), Dn.hi)
            Dm = Dn
            Dn = interval(Dn.lo, new_hi)
        else
            new_lo = max(Dn.hi - (d / abs(Li)), Dn.lo)
            Dm = Dn
            Dn = interval(new_lo, Dn.hi)
        end
    end
    hi = fT.pol(fT.dom - fT.x0).hi
    return interval(bound.lo, hi) + fT.rem
end

"""
    linear_dominated_bounder(fT::TaylorModelN, ϵ=1e-3::Float64, max_iter=5::Int)

Compute a tighter polynomial bound for the Taylor model `fT` by the linear
dominated bounder algorithm. The linear dominated algorithm is applied until
the bound of `fT` gets tighter than `ϵ` or the number of steps reachs `max_iter`.
The returned bound corresponds to the improved polynomial bound with the remainder
of the `TaylorModelN` included.
"""
function linear_dominated_bounder(fT::TaylorModelN{R, S, T}; ϵ=1e-5, max_iter=5) where {R, S, T}
    d = 1.
    Dn = fT.dom
    Dm = Dn
    Pm = TaylorN(deepcopy(fT.pol.coeffs))
    bound = interval(0.)
    n_iter = 0
    while d > ϵ && n_iter < max_iter
        if n_iter == 0
            x0 = Array(mid(Dn) - mid(fT.x0))
        else
            x0 = Array(mid(Dn) - mid(Dm))
        end
        c = mid(Dn)
        update!(Pm, x0)
        linear_part = Pm[1]
        if S <: Interval
            linear_coeffs = Float64[mid(coeff) for coeff in linear_part.coeffs]
        else
            linear_coeffs = Float64[coeff for coeff in linear_part.coeffs]
        end
        linear = TaylorN(deepcopy(Pm.coeffs[1:2]), Pm.order)
        non_linear = Pm - linear
        centered_domain = Dn .- c
        I1 = linear(centered_domain)
        Ih = non_linear(centered_domain)
        bound = I1.lo + Ih
        d = diam(bound)
        n_iter += 1
        new_boxes = Interval[]
        for (idx, box) in enumerate(Dn)
            Li = linear_coeffs[idx]
            if Li == 0
                Dni = box
            elseif Li < 0
                lo = max(box.hi - (d / abs(Li)), box.lo)
                Dni = interval(lo, box.hi)
            else
                hi = min(box.lo + (d / abs(Li)), box.hi)
                Dni = interval(box.lo, hi)
            end
            push!(new_boxes, Dni)
        end
        Dm = Dn
        Dn = IntervalBox(new_boxes...)
    end
    hi = fT.pol(fT.dom - fT.x0).hi
    return interval(bound.lo, hi) + fT.rem
end

"""
    quadratic_fast_bounder(fT::TaylorModel)

Compute a *tighter* polynomial bound by the quadratic fast bounder.
The returned bound corresponds to the "improved" polynomial bound
with the remainder of the `TaylorModel` included. This "improved" bound
can be one of the following two:
    1) An improved bound: if the domain of `fT` has a local minimizer,
       then an improved bound is returned.
    2) Original TaylorModel bound: if the local minimizer criteria is not
       satisfied, then the original bound of `fT` is returned.

This algorithm is a slight modification to the Makino & Berz algorithm.
For this algorithm the linear part is bounded by solving a simple
set of linear equations (compared to the iterative procedure by Makino & Berz).
"""
function quadratic_fast_bounder(fT::TaylorModel1)
    P = fT.pol
    bound_tm = fT(fT.dom - fT.x0)
    if signbit(P[2])
        return bound_tm
    else
        x0 = -P[1] / (2 * P[2])
        x = Taylor1(fT.pol.order)
        Qx0 = (x - x0) * P[2] * (x - x0)
        bound_qfb = (P - Qx0)(fT.dom - fT.x0)
        hi = P(fT.dom - fT.x0).hi
        bound_qfb = interval(bound_qfb.lo, hi) + fT.rem
        bound = bound_qfb ∩ bound_tm
        return bound
    end
end

function quadratic_fast_bounder(fT::TaylorModelN)
    P = fT.pol
    H = Matrix(TaylorSeries.hessian(P))
    bound_tm = fT(fT.dom - fT.x0)
    if isposdef(H)
        P1 = -P[1].coeffs
        xn = H \ P1
        x = set_variables("x", numvars=length(xn))
        Qxn = 0.5 * (x - xn)' * H * (x - xn)
        bound_qfb = (P - Qxn)(fT.dom - fT.x0)
        hi = P(fT.dom - fT.x0).hi
        bound_qfb = interval(bound_qfb.lo, hi) + fT.rem
        bound = bound_qfb ∩ bound_tm
        return bound
    else
        return bound_tm
    end
end
