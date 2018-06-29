# rpa_functions.jl

# α_mid=0.484375 is used to preferentially round-down the mid point
# when the mid point is not exactly representable. Use for conversion
const α_mid = 0.484375 # corrsponds to 31/64; to preferentiably round down


"""
   _rpaar(f::Function, x0::Interval, ii::Interval, _order::Int)

Rigurous polynomial approximation (RPA) with absolute remainder
for the function `f` on the interval `ii`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
by `boundarem`(@ref)
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function _rpaar(f::Function, x0::Interval{T}, ii::Interval{T}, _order::Int) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( ii+Taylor1(Interval{T}, _order+1) )
    Δ = boundarem(f, polf, polfI, x0, ii)
    return TM1AbsRem(polf, Δ, x0, ii)
end


"""
   _rparr(f::Function, x0::Interval, ii::Interval, _order::Int)

Rigurous polynomial approximation (RPA) with relative remainder
for the function `f` on the interval `ii`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
by `boundrrem`(@ref)
exploiting monotonicity if possible, otherwise, it uses the Lagrange
coefficient.

"""
function _rparr(f::Function, x0::Interval{T}, ii::Interval{T}, _order::Int) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( ii+Taylor1(Interval{T}, _order+2) )
    Δ = boundrrem(f, polf, polfI, x0, ii)
    return TM1RelRem(polf, Δ, x0, ii)
end


"""
   rpa(g::Function, tmf::TM1AbsRem)
   rpa(g::Function, tmf::TMNAbsRem)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function rpa(g::Function, tmf::TM1AbsRem)
    _order = get_order(tmf)

    # Avoid overestimations:
    if tmf == TM1AbsRem(_order, tmf.x0, tmf.iI)
        # ... if `tmf` is the independent variable
        return _rpaar(g, tmf.x0, tmf.iI, _order)
    elseif tmf == TM1AbsRem(constant_term(tmf.pol), _order, tmf.x0, tmf.iI)
        # ... in case `tmf` is a simple constant polynomial
        range_g = bound_taylor1(g(tmf.pol), tmf.iI-tmf.x0) + remainder(tmf)
        return TM1AbsRem(range_g, _order, tmf.x0, tmf.iI)
    end

    f_pol = tmf.pol
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    iI = tmf.iI

    # Range of tmf including remainder (Δf)
    range_tmf = bound_taylor1(f_pol, iI-x0) + Δf

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    tmg = _rpaar(g, f_pol0, range_tmf, _order)

    # Use original independent variable
    tm1 = tmf - f_pol0
    tmres = tmg( tm1 )

    # Final remainder
    Δ = remainder(tmres) + remainder(tmg)
    return TM1AbsRem(tmres.pol, Δ, x0, iI)
end

function rpa(g::Function, tmf::TMNAbsRem{N,T,S}) where {N,T,S}
    _order = get_order(tmf)

    # # Avoid overestimations
    # if tmf == TMNAbsRem(constant_term(tmf.pol), _order, tmf.x0, tmf.iI)
    #     # ... in case `tmf` is a simple constant polynomial
    #     range_g = (g(tmf.pol))(tmf.iI-tmf.x0) + remainder(tmf)
    #     return TMNAbsRem(range_g, _order, tmf.x0, tmf.iI)
    # else
    #     v = get_variables(T, _order)
    #     any( tmf.pol .== v ) && _rpaar(g, tmf.x0, tmf.iI, _order)
    # end

    f_pol = tmf.pol
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    iI = tmf.iI

    # Range of tmf including remainder (Δf)
    range_tmf = f_pol(iI-x0) + Δf

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    # Note that tmg is a TM1AbsRem !!
    tmg = _rpaar(g, f_pol0, range_tmf, _order)

    # Use original independent variable
    tm1 = tmf - f_pol0
    tmres = tmg( tm1 )

    # Final remainder
    Δ = remainder(tmres) + remainder(tmg)
    return TMNAbsRem(tmres.pol, Δ, x0, iI)
end


"""
   rpa(g::Function, tmf::TM1RelRem)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function rpa(g::Function, tmf::TM1RelRem)
    _order = get_order(tmf)

    # Avoid overestimations:
    if tmf == TM1RelRem(_order, tmf.x0, tmf.iI)
        # ... if `tmf` is the independent variable
        return _rparr(g, tmf.x0, tmf.iI, _order)
    elseif tmf == TM1RelRem(constant_term(tmf.pol), _order, tmf.x0, tmf.iI)
        # ... in case `tmf` is a simple constant polynomial
        range_g = bound_taylor1(g(tmf.pol), tmf.iI-tmf.x0) + remainder(tmf)
        return TM1RelRem(range_g, _order, tmf.x0, tmf.iI)
    end

    f_pol = tmf.pol
    Δf = remainder(tmf)
    x0 = tmf.x0
    iI = tmf.iI

    # Range of tmf including remainder (Δf)
    range_tmf = bound_taylor1(f_pol, iI-x0) + Δf * (iI-x0)^(_order+1)

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    tmg = _rparr(g, constant_term(f_pol), range_tmf, _order)
    tm1 = tmf - constant_term(f_pol)
    tmres = tmg( tm1 )

    tmn = TM1RelRem(Taylor1(copy(tm1.pol.coeffs)), tm1.rem, tm1.x0, tm1.iI)
    for i = 1:_order
        tmn = tmn * tm1
    end
    Δ = remainder(tmres) + remainder(tmn) * remainder(tmg)
    return TM1RelRem(tmres.pol, Δ, x0, iI)
end



"""
    fp_rpa(tm::TM1AbsRem{Interval{T},T})
    fp_rpa(tm::TM1RelRem{Interval{T},T})

Convert a `tm` TaylorModel to a TaylorModel whose polynomial coefficients
of type `T<:Real`. The accumulated error is added to the remainder. The
mid point of the expansion interval is preferentially rounded down
if it is not an exactly representable value.

"""
function fp_rpa(tm::TM1AbsRem{Interval{T},T}) where {T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    ii = tm.iI
    order = get_order(tm)
    ξ0 = mid(x0, α_mid)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    for ind=0:order
        t[ind] = mid(fT[ind], α_mid)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    δ = b(ii-x0)
    Δ = Δ + δ
    return TM1AbsRem(t, Δ, x0, ii)
end

function fp_rpa(tm::TM1RelRem{Interval{T},T}) where {T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    ii = tm.iI
    order = get_order(tm)
    ξ0 = mid(x0, α_mid)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    for ind=0:order
        t[ind] = mid(fT[ind], α_mid)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    δ = b(ii-x0)
    # Is the following correct for TM1RelRem?
    Δ = Δ + δ
    return TM1RelRem(t, Δ, x0, ii)
end


# Elementary functions
fnlist = (:inv, :sqrt, :exp, :log, :sin, :cos, :tan,
    :asin, :acos, :atan, :sinh, :cosh, :tanh)

for fn in fnlist
    for TM in tupleTMs
        @eval $fn(tm::$TM) = rpa($fn, tm)
    end
    @eval $fn(tm::TMNAbsRem) = rpa($fn, tm)
end
