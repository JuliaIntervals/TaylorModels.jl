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
    bb = sup(fTIend) < 0 || inf(fTIend) > 0
    return _monot_bound_remainder(TaylorModel1, Val(bb), f, polf, polfI, x0, I)
end
function bound_remainder(::Type{TaylorModel1}, f::Function, polf::Taylor1{TaylorN{T}}, polfI::Taylor1, x0, I::Interval) where {T}
    # The domain of the TaylorN part is assumed to be
    # the normalized symmetric box
    _order = get_order(polf) + 1
    fTIend = polfI[_order]
    bb = sup(fTIend) < 0 || inf(fTIend) > 0
    return _monot_bound_remainder(TaylorModel1, Val(bb), f, polf, polfI, x0, I)
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
    a = Interval(inf(I))
    b = Interval(sup(I))
    bb = (sup(fTIend) < 0 || inf(fTIend) > 0) && isempty(a ∩ x0) && isempty(b ∩ x0)
    return _monot_bound_remainder(RTaylorModel1, Val(bb), f, polf, polfI, x0, I)
end


"""
    _monot_bound_remainder(::Type{TaylorModel1}, ::Val{true}, f::Function, polf::Taylor1, polfI::Taylor1, x0, I::Interval)

Computes the remainder exploiting monotonicity; see Prop 2.2.1 in Mioara Joldes' PhD thesis (pp 52).
"""
@inline function _monot_bound_remainder(::Type{TaylorModel1}, ::Val{true}, f::Function,
        polf::Taylor1, polfI::Taylor1, x0, I::Interval)
    # Absolute remainder is monotonic
    a = Interval(inf(I))
    b = Interval(sup(I))
    Δlo = f(a) - polf(a-x0)
    # Δlo = f(a) - bound_taylor1(polf, a-x0)
    Δhi = f(b) - polf(b-x0)
    # Δhi = f(b) - bound_taylor1(polf, b-x0)
    Δx0 = f(x0) - polf[0]
    return hull(Δlo, Δx0, Δhi)
end
@inline function _monot_bound_remainder(::Type{TaylorModel1}, ::Val{true}, f::Function,
        polf::Taylor1{TaylorN{T}}, polfI::Taylor1, x0, I::Interval) where {T}
    # Absolute remainder is monotonic
    a = Interval(inf(I))
    b = Interval(sup(I))
    symIbox = IntervalBox(-1 .. 1, get_numvars())
    Δlo = (f(a) - polf(a-x0))(symIbox)
    # Δlo = f(a) - bound_taylor1(polf, a-x0)
    Δhi = (f(b) - polf(b-x0))(symIbox)
    # Δhi = f(b) - bound_taylor1(polf, b-x0)
    Δx0 = (f(x0) - polf[0])(symIbox)
    return hull(Δlo, Δx0, Δhi)
end
"""
    _monot_bound_remainder(::Type{TaylorModel1}, ::Val{false}, f::Function, polf::Taylor1, polfI::Taylor1, x0, I::Interval)

Computes the remainder using Lagrange bound
"""
@inline function _monot_bound_remainder(::Type{TaylorModel1}, ::Val{false}, f::Function,
        polf, polfI, x0, I)
    _order = get_order(polf) + 1
    fTIend = polfI[_order]
    # Lagrange bound
    return fTIend * (I-x0)^_order
end

"""
    _monot_bound_remainder(::Type{RTaylorModel1}, ::Val{true}, f::Function, polf::Taylor1, polfI::Taylor1, x0, I::Interval)

Computes the remainder exploiting monotonicity; see Prop 2.3.7 in Mioara Joldes' PhD thesis (pp 67).
"""
@inline function _monot_bound_remainder(::Type{RTaylorModel1}, ::Val{true}, f::Function,
        polf::Taylor1, polfI::Taylor1, x0, I::Interval)

    _order = get_order(polf) + 1
    a = Interval(inf(I))
    b = Interval(sup(I))
    # Error is monotonic
    denom_lo = (a-x0)^_order
    Δlo = f(a) - polf(a-x0)
    # Δlo = f(a) - bound_taylor1(polf, a-x0)
    Δlo = Δlo / denom_lo
    denom_hi = (b-x0)^_order
    Δhi = f(b) - polf(b-x0)
    # Δhi = f(b) - bound_taylor1(polf, b-x0)
    Δhi = Δhi / denom_hi
    return hull(Δlo, Δhi)
end
@inline function _monot_bound_remainder(::Type{RTaylorModel1}, ::Val{true}, f::Function,
        polf::Taylor1{TaylorN{T}}, polfI::Taylor1, x0, I::Interval) where {T}

    _order = get_order(polf) + 1
    a = Interval(inf(I))
    b = Interval(sup(I))
    symIbox = IntervalBox(-1 .. 1, get_numvars())
    # Error is monotonic
    denom_lo = (a-x0)^_order
    Δlo = (f(a) - polf(a-x0))(symIbox)
    # Δlo = f(a) - bound_taylor1(polf, a-x0)
    Δlo = Δlo / denom_lo
    denom_hi = (b-x0)^_order
    Δhi = (f(b) - polf(b-x0))(symIbox)
    # Δhi = f(b) - bound_taylor1(polf, b-x0)
    Δhi = Δhi / denom_hi
    return hull(Δlo, Δhi)
end
@inline function _monot_bound_remainder(::Type{RTaylorModel1}, ::Val{false}, f::Function,
        polf::Taylor1, polfI::Taylor1, x0, I::Interval)
    _order = get_order(polf) + 1
    return polfI[_order]
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
    rangepoly = hull( fT(Interval(inf(I))), fT(Interval(sup(I))) )
    @inbounds for ind in 1:num_roots
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
    I_lo = inf(I)
    I_hi = sup(I)
    if inf(fTd(I)) ≥ 0
        return Interval(fT(I_lo), fT(I_hi))
    elseif sup(fTd(I)) ≤ 0
        return Interval(fT(I_hi), fT(I_lo))
    end
    return fT(I)
end
function bound_taylor1(fT::Taylor1{Interval{T}}, fTd::Taylor1{Interval{T}},
        I::Interval{S}) where {T,S}
    #
    I_lo = inf(I)
    I_hi = sup(I)
    if inf(fTd(I)) ≥ 0
        return hull(fT(I_lo), fT(I_hi))
    elseif sup(fTd(I)) ≤ 0
        return hull(fT(I_hi), fT(I_lo))
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
function linear_dominated_bounder(fT::TaylorModel1{T, S}; ϵ=1e-3, max_iter=5) where {T, S}
    d = one(T)
    dom = domain(fT)
    x0 = expansion_point(fT)
    pol = polynomial(fT)
    Pm = deepcopy(pol)
    bound = zero(x0)
    hi = sup(pol(dom - x0))

    n_iter = 0
    while d > ϵ && n_iter < max_iter
        x0 = mid(dom - x0)
        c = mid(dom)
        update!(Pm, x0)
        non_linear = nonlinear_polynomial(Pm)
        linear = Pm - non_linear
        if T <: Interval
            Li = mid(linear[1])
        else
            Li = linear[1]
        end
        I1 = linear(dom - c)
        Ih = non_linear(dom - c)
        bound = inf(I1) + Ih
        d = diam(bound)
        n_iter += 1
        dom_lo = inf(dom)
        dom_hi = sup(dom)
        if Li == 0
            break
        elseif Li > 0
            new_hi = min(dom_lo + (d / abs(Li)), dom_hi)
            x0 = dom
            dom = Interval(dom_lo, new_hi)
        else
            new_lo = max(dom_hi - (d / abs(Li)), dom_lo)
            x0 = dom
            dom = Interval(new_lo, dom_hi)
        end
    end

    return Interval(inf(bound), hi) + remainder(fT)
end

"""
    linear_dominated_bounder(fT::TaylorModelN, ϵ=1e-3::Float64, max_iter=5::Int)

Compute a tighter polynomial bound for the Taylor model `fT` by the linear
dominated bounder algorithm. The linear dominated algorithm is applied until
the bound of `fT` gets tighter than `ϵ` or the number of steps reachs `max_iter`.
The returned bound corresponds to the improved polynomial bound with the remainder
of the `TaylorModelN` included.
"""
function linear_dominated_bounder(fT::TaylorModelN{N,T,S}; ϵ=1e-5, max_iter=5) where {N, T, S}
    d = one(T)
    dom = domain(fT)
    x0 = expansion_point(fT)
    pol = polynomial(fT)
    Pm = deepcopy(pol)
    bound = zero(Interval{S})
    pol_hi = sup(pol(dom - x0))

    n_iter = 0
    x00 = Array{S}(undef, N)
    new_boxes = fill(bound, N)
    linear_coeffs = Array{Float64}(undef, N)
    while d > ϵ && n_iter < max_iter
        x00 .= mid(dom - x0)
        c = mid(dom)
        update!(Pm, x00)
        linear_part = Pm[1]
        if T <: Interval
            @. linear_coeffs = mid(linear_part.coeffs)
        else
            @. linear_coeffs = linear_part.coeffs
        end
        non_linear = nonlinear_polynomial(Pm)
        linear = Pm - non_linear
        centered_domain = dom .- c
        I1 = linear(centered_domain)
        Ih = non_linear(centered_domain)
        bound = inf(I1) + Ih
        d = diam(bound)
        n_iter += 1
        for (idx, box) in enumerate(dom)
            Li = linear_coeffs[idx]
            box_lo = inf(box)
            box_hi = sup(box)
            if Li == 0
                domi = box
            elseif Li < 0
                lo = max(box_hi - (d / abs(Li)), box_lo)
                domi = Interval(lo, box_hi)
            else
                hi = min(box_lo + (d / abs(Li)), box_hi)
                domi = Interval(box_lo, hi)
            end
            new_boxes[idx] = domi
        end
        x0 = dom
        dom = IntervalBox(new_boxes...)
    end

    return Interval(inf(bound), pol_hi) + remainder(fT)
end

"""
    quadratic_fast_bounder(fT::TaylorModel1)

Compute a *tighter* polynomial bound by the quadratic fast bounder.
The returned bound corresponds to the "improved" polynomial bound
with the remainder of the `TaylorModel1` included. This "improved" bound
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
    P = polynomial(fT)
    dom = domain(fT)
    x00 = expansion_point(fT)
    bound_tm = fT(dom - x00)
    signbit(P[2]) && return bound_tm

    cent_dom = dom - x00
    x0 = -P[1] / (2 * P[2])
    x = Taylor1(get_order(P))
    Qx0 = (x - x0) * P[2] * (x - x0)
    bound_qfb = (P - Qx0)(cent_dom)
    hi = sup(P(cent_dom))
    bound_qfb = Interval(inf(bound_qfb), hi) + remainder(fT)

    return bound_qfb ∩ bound_tm
end

"""
    quadratic_fast_bounder(fT::TaylorModelN)

Compute a *tighter* polynomial bound by the quadratic fast bounder.
The returned bound corresponds to the "improved" polynomial bound
with the remainder of the `TaylorModelN` included. This "improved" bound
can be one of the following two:
    1) An improved bound: if the domain of `fT` has a local minimizer,
       then an improved bound is returned.
    2) Original TaylorModel bound: if the local minimizer criteria is not
       satisfied, then the original bound of `fT` is returned.

This algorithm is a slight modification to the Makino & Berz algorithm.
For this algorithm the linear part is bounded by solving a simple
set of linear equations (compared to the iterative procedure by Makino & Berz).
"""
function quadratic_fast_bounder(fT::TaylorModelN)
    @assert get_numvars() == get_numvars(fT)
    P = polynomial(fT)
    dom = domain(fT)
    x0 = expansion_point(fT)
    H = Matrix(TaylorSeries.hessian(P))
    bound_tm = fT(dom - x0)
    if isposdef(H)
        P1 = -P[1].coeffs
        xn = H \ P1
        x = get_variables()#set_variables("x", numvars=length(xn))
        Qxn = 0.5 * (x - xn)' * H * (x - xn)
        bound_qfb = (P - Qxn)(dom - x0)
        hi = sup(P(dom - x0))
        bound_qfb = interval(inf(bound_qfb), hi) + remainder(fT)
        bound = bound_qfb ∩ bound_tm
        return bound
    end

    return bound_tm
end
