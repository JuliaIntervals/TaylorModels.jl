# rpa_functions.jl

# α_mid=0.484375 is used to preferentially round-down the mid point
# when the mid point is not exactly representable.
const α_mid = 0.484375 # == 31/64


"""
   _rpaar(f::Function, x0::Interval, I::Interval, _order::Integer)

Rigurous polynomial approximation (RPA) with absolute remainder
for the function `f` on the interval `I`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
by `bound_absrem`(@ref)
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function _rpaar(f::Function, x0::Interval{T}, I::Interval{T}, _order::Integer) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( I+Taylor1(Interval{T}, _order+1) )
    Δ = bound_absrem(f, polf, polfI, x0, I)
    return TaylorModel1(polf, Δ, x0, I)
end
function _rpaar(f::Function, x0::T, I::Interval{T}, _order::Integer) where {T}
    polf  = f( x0+Taylor1(T, _order) )
    polfI = f( I+Taylor1(Interval{T}, _order+1) )
    Δ = bound_absrem(f, polf, polfI, Interval(x0), I)
    return TaylorModel1(polf, Δ, Interval(x0), I)
end


"""
   _rparr(f::Function, x0::Interval, I::Interval, _order::Integer)

Rigurous polynomial approximation (RPA) with relative remainder
for the function `f` on the interval `I`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
by `bound_relrem`(@ref)
exploiting monotonicity if possible, otherwise, it uses the Lagrange
coefficient.

"""
function _rparr(f::Function, x0, I::Interval{T}, _order::Integer) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( I+Taylor1(Interval{T}, _order+2) )
    Δ = bound_relrem(f, polf, polfI, x0, I)
    return RTaylorModel1(polf, Δ, x0, I)
end


"""
   rpa(g::Function, tmf::TaylorModel1)
   rpa(g::Function, tmf::TaylorModelN)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function rpa(g::Function, tmf::TaylorModel1)
    _order = get_order(tmf)

    # # Avoid overestimations:
    # if tmf == TaylorModel1(_order, tmf.x0, tmf.dom)
    #     # ... if `tmf` is the independent variable
    #     return _rpaar(g, tmf.x0, tmf.dom, _order)
    # elseif tmf == TaylorModel1(constant_term(tmf.pol), _order, tmf.x0, tmf.dom)
    #     # ... in case `tmf` is a simple constant polynomial
    #     range_g = bound_taylor1(g(tmf.pol), tmf.dom-tmf.x0) + remainder(tmf)
    #     return TaylorModel1(range_g, _order, tmf.x0, tmf.dom)
    # end

    f_pol = tmf.pol
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    I = tmf.dom

    # Range of tmf including remainder (Δf)
    # range_tmf = bound_taylor1(f_pol, I-x0) + Δf
    range_tmf = f_pol(I-x0) + Δf

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    tmg = _rpaar(g, f_pol0, range_tmf, _order)

    # Use original independent variable
    # tm1 = copy(tmf)
    # tm1.pol[0] = zero(f_pol0)
    tm1 = tmf - f_pol0   # OVER-ESTIMATION; IMPROVE
    tmres = tmg( tm1 )

    # Final remainder
    Δ = remainder(tmres) + remainder(tmg)
    return TaylorModel1(tmres.pol, Δ, x0, I)
end

function rpa(g::Function, tmf::TaylorModelN{N,T,S}) where {N,T,S}
    _order = get_order(tmf)

    # # Avoid overestimations
    # if tmf == TaylorModelN(constant_term(tmf.pol), _order, tmf.x0, tmf.dom)
    #     # ... in case `tmf` is a simple constant polynomial
    #     range_g = (g(tmf.pol))(tmf.dom-tmf.x0) + remainder(tmf)
    #     return TaylorModelN(range_g, _order, tmf.x0, tmf.dom)
    # else
    #     v = get_variables(T, _order)
    #     any( tmf.pol .== v ) && _rpaar(g, tmf.x0, tmf.dom, _order)
    # end

    f_pol = tmf.pol
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    I = tmf.dom

    # Range of tmf including remainder (Δf)
    range_tmf = f_pol(I-x0) + Δf

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    # Note that tmg is a TaylorModel1 !!
    tmg = _rpaar(g, f_pol0, range_tmf, _order)

    # Use original independent variable
    # tm1 = copy(tmf)
    # tm1.pol[0] = zero(f_pol0)
    tm1 = tmf - f_pol0   # OVER-ESTIMATION; IMPROVE
    tmres = tmg( tm1 )

    # Final remainder
    Δ = remainder(tmres) + remainder(tmg)
    return TaylorModelN(tmres.pol, Δ, x0, I)
end


"""
   rpa(g::Function, tmf::RTaylorModel1)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function rpa(g::Function, tmf::RTaylorModel1)
    _order = get_order(tmf)

    # # Avoid overestimations:
    # if tmf == RTaylorModel1(_order, tmf.x0, tmf.dom)
    #     # ... if `tmf` is the independent variable
    #     return _rparr(g, tmf.x0, tmf.dom, _order)
    # elseif tmf == RTaylorModel1(constant_term(tmf.pol), _order, tmf.x0, tmf.dom)
    #     # ... in case `tmf` is a simple constant polynomial
    #     range_g = bound_taylor1(g(tmf.pol), tmf.dom-tmf.x0) + remainder(tmf)
    #     return RTaylorModel1(range_g, _order, tmf.x0, tmf.dom)
    # end

    f_pol = tmf.pol
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    I = tmf.dom

    # Range of tmf including remainder (Δf)
    # range_tmf = bound_taylor1(f_pol, I-x0) + Δf * (I-x0)^(_order+1)
    range_tmf = f_pol(I-x0) + Δf * (I-x0)^(_order+1)

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    tmg = _rparr(g, f_pol0, range_tmf, _order)
    # tm1 = copy(tmf)
    # tm1.pol[0] = zero(f_pol0)
    tm1 = tmf - f_pol0   # OVER-ESTIMATION; IMPROVE
    tmres = tmg( tm1 )

    tmn = RTaylorModel1(Taylor1(copy(tm1.pol.coeffs)), tm1.rem, tm1.x0, tm1.dom)
    for i = 1:_order
        tmn = tmn * tm1
    end
    Δ = remainder(tmres) + remainder(tmn) * remainder(tmg)
    return RTaylorModel1(tmres.pol, Δ, x0, I)
end



"""
    fp_rpa(tm::TaylorModel1{Interval{T},T})
    fp_rpa(tm::RTaylorModel1{Interval{T},T})

Convert a `tm` TaylorModel1 to a TaylorModel1 whose polynomial coefficients
are `Float64`. The accumulated error is added to the remainder. The
mid point of the expansion interval is preferentially rounded down
if it is not an exactly representable value.

"""
function fp_rpa(tm::TaylorModel1{Interval{T},T}) where {T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    I = tm.dom
    order = get_order(tm)
    # ξ0 = mid(x0, α_mid)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    @inbounds for ind=order:-1:0
        t[ind] = mid(fT[ind], α_mid)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    δ = b(I-x0)
    Δ = Δ + δ
    return TaylorModel1(t, Δ, x0, I)
end

fp_rpa(tm::TaylorModel1{T, T}) where {T} = tm

function fp_rpa(tm::RTaylorModel1{Interval{T},T}) where {T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    I = tm.dom
    order = get_order(tm)
    # ξ0 = mid(x0, α_mid)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    @inbounds for ind=order:-1:0
        t[ind] = mid(fT[ind], α_mid)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    δ = b(I-x0)
    # Is the following correct for RTaylorModel1?
    Δ = Δ + δ
    return RTaylorModel1(t, Δ, x0, I)
end

fp_rpa(tm::RTaylorModel1{T, T}) where {T} = tm

function fp_rpa(tm::TaylorModelN{N,Interval{T},T}) where {N,T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    I = tm.dom
    order = get_order(tm)

    b = zero(fT)
    t = TaylorN([HomogeneousPolynomial(zeros(T,N))], order)
    for ind = 0:order
        @inbounds for homPol in 1:length(fT[ind])
            t[ind][homPol] = mid(fT[ind][homPol])
            b[ind][homPol] = fT[ind][homPol] - t[ind][homPol]
        end
    end
    Δ = Δ + b(I-x0)
    return TaylorModelN(t, Δ, x0, I)
end

fp_rpa(tm::TaylorModelN{N, T, T}) where {N, T} = tm

# Elementary functions
fnlist = (:inv, :sqrt, :exp, :log, :sin, :cos, :tan,
    :asin, :acos, :atan, :sinh, :cosh, :tanh)

for fn in fnlist
    for TM in tupleTMs
        @eval $fn(tm::$TM) = rpa($fn, tm)
    end
    @eval $fn(tm::TaylorModelN) = rpa($fn, tm)
end
