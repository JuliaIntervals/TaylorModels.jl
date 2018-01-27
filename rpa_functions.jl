# rpa_functions.jl

"""
   rpa(f::Function, x0::Interval, ii::Interval, _order::Int)

Rigurous polynomial approximation (RPA) for the function `f` on
the interval `ii`,  using a Taylor expansion around the *interval* `x0`
of order `_order`. The bound is computed exploiting monotonicity if
possible, otherwise, it uses Laplace bound.

"""
function rpa(f::Function, x0::Interval{T}, ii::Interval{T}, _order::Int) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( ii+Taylor1(Interval{T}, _order+1) )
    Δ = bound_arem(f, polf, polfI[_order+1], x0, ii)
    return TMAbsRem(polf, Δ, x0, ii)
end

function rpa(g::Function, tmf::TMAbsRem)
    _order = get_order(tmf)

    # Do not overestimate if `tmf` is the independent variable
    tmf == TMAbsRem(_order, tmf.x0, tmf.iI) && return rpa(g, tmf.x0, tmf.iI, _order)

    f_pol = tmf.pol
    Δf = remainder(tmf)
    x0 = tmf.x0
    iI = tmf.iI

    # Range of tmf including remainder (Δf)
    range_tmf = bound_taylor1(f_pol, iI-x0) + Δf

    # Compute RPA for `g`, around f_pol[0], over range_tmf
    tmg = rpa(g, f_pol[0], range_tmf, _order)
    tm1 = tmf - f_pol[0]
    tmres = tmg( tm1 )
    Δ = remainder(tmres) + remainder(tmg)
    return TMAbsRem(tmres.pol, Δ, x0, iI)
end


# evaluate, and function-like evaluation for TMAbsRem
function evaluate(tmg::TMAbsRem, tmf::TMAbsRem)
    _order = get_order(tmf)
    @assert _order == get_order(tmg)

    # R = typeof(tmg.pol[0] * tmf.pol[0])
    tmres = TMAbsRem(tmg.pol[_order], _order, tmf.x0, tmf.iI)
    @inbounds for k = _order-1:-1:0
        tmres = tmres * tmf
        tmres = tmres + TMAbsRem(tmg.pol[k], _order, tmf.x0, tmf.iI)
    end
    return tmres
end

(tm::TMAbsRem)(x::TMAbsRem) = evaluate(tm, x)



"""
    rpafp(tm::TMAbsRem{T,S})

Convert a `tm::TMAbsRem{T,S}` to a T-type RPA. It returns the `Taylor1{T}`
polynomial, the accumulated aremor `Δ::Interval{S}`, and `ξ0` which is
the mid point about which the expansion is obtained. If ξ0 is not exactly
representable, it returns *preferentiably* a rounded-down value.
This function is primarily used for plotting.

"""
rpafp(tm::TMAbsRem{T}) where {T} = _rpafp(tm.pol, tm.arem, tm.x0, tm.iI)

function _rpafp(fT::Taylor1{Interval{T}}, Δ::Interval{S},
        x0::Interval{S}, ii::Interval{S}) where {T,S}
    order = fT.order
    # α=0.484375 is used to get preferentially round-down of the mid point
    # when the mid point is not exactly representable
    α = 0.484375
    ξ0 = mid(x0, α)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    for ind=0:order
        t[ind] = mid(fT[ind], α)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    # if norm(b, Inf) < eps()
    #     δ = b(ii-x0)
    # else
    #     δ = bound_taylor1(b, ii-x0)   ## Used alone may yield problems
    # end
    δ = b(ii-x0)
    Δ = Δ + δ
    # return TMAbsRem(t, Δ, interval(ξ0), ii)
    return t, Δ, ξ0
end


# Elementary functions
fnlist = (:inv, :sqrt, :exp, :log, :sin, :cos, :tan,
    :asin, :acos, :atan, :sinh, :cosh, :tanh)

for fn in fnlist
    @eval $fn(tm::TMAbsRem) = rpa($fn, tm)
end
