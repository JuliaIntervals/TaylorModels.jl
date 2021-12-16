# rpa_functions.jl

# α_mid=0.484375 is used to preferentially round-down the mid point
# when the mid point is not exactly representable.
const α_mid = 0.484375 # == 31/64

"""
```
   _rpa(::Type{TaylorModel1}, f::Function, x0::Interval, I::Interval, _order::Integer)
   _rpa(::Type{RTaylorModel1}, f::Function, x0::Interval, I::Interval, _order::Integer)
```

Rigurous polynomial approximation (RPA) with absolute/relative remainder
for the function `f` on the interval `I`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
by `bound_remainder`(@ref)
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

""" _rpa

for TM in tupleTMs
    @eval begin
        function _rpa(::Type{$TM}, f::Function, x0::Interval{T}, I::Interval{T},
                _order::Integer) where {T}

            polf  = f( Taylor1([x0,one(x0)], _order) )
            polfI = f( Taylor1([I,one(I)], _order+1+($TM==RTaylorModel1) ) )
            Δ = bound_remainder($TM, f, polf, polfI, x0, I)
            return $TM(polf, Δ, x0, I)
        end

        function _rpa(::Type{$TM}, f::Function, x0::TaylorN{T}, I::Interval{T},
                      _order::Integer) where {T}

            polf = f(Taylor1([x0, one(x0)], _order))
            polfI = f( Taylor1([I,one(I)], _order+1+($TM==RTaylorModel1)))
            x0I = Interval(constant_term(x0))
            Δ = bound_remainder($TM, f, polf, polfI, x0I, I)
            return $TM(polf, Δ, x0I, I)
        end

        function _rpa(::Type{$TM}, f::Function, x0::T, I::Interval{T},
                _order::Integer) where {T}

            polf  = f( Taylor1([x0,one(x0)], _order) )
            polfI = f( Taylor1([I,one(I)], _order+1+($TM==RTaylorModel1)) )
            x0I = Interval(x0)
            Δ = bound_remainder($TM, f, polf, polfI, x0I, I)
            return $TM(polf, Δ, x0I, I)
        end
    end
end

function _rpa(::Type{TaylorModel1}, f::Function, x0::TaylorModelN, I::Interval,
        _order::Integer)
    polf  = f( Taylor1([x0, one(x0)], _order) )
    polfI = f( Taylor1([I, one(I)], _order+1) )
    Δ = bound_remainder(TaylorModel1, f, polf, polfI, I, I)
    return TaylorModel1(polf, Δ, x0(x0.x0), I)
end


"""
   rpa(g::Function, tmf::TaylorModel1)
   rpa(g::Function, tmf::RTaylorModel1)
   rpa(g::Function, tmf::TaylorModelN)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute/relative remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

""" rpa

for TM in tupleTMs
    @eval begin
        function rpa(g::Function, tmf::$TM)
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
            f_pol = polynomial(tmf)
            f_pol0 = constant_term(f_pol)
            Δf = remainder(tmf)
            x0 = tmf.x0
            I = domain(tmf)

            # Range of tmf including remainder (Δf)
            if $TM == TaylorModel1
                # range_tmf = bound_taylor1(f_pol, I-x0) + Δf
                range_tmf = f_pol(I-x0) + Δf
            else
                # range_tmf = bound_taylor1(f_pol, I-x0) + Δf * (I-x0)^(_order+1)
                range_tmf = f_pol(I-x0) + Δf * (I-x0)^(_order+1)
            end

            # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
            tmg = _rpa($TM, g, f_pol0, range_tmf, _order)
            if tmf == $TM(_order, x0, I)  # indep variable
                return tmg
            end

            # Use original independent variable
            tm1 = tmf - f_pol0   # OVER-ESTIMATION; IMPROVE
            tmres = tmg( tm1 )

            # Final remainder
            if $TM == TaylorModel1
                Δ = remainder(tmres) + remainder(tmg)
            else
                tmn = RTaylorModel1(Taylor1(copy(tm1.pol.coeffs)), tm1.rem, x0, I)
                for i = 1:_order
                    tmn = tmn * tmf
                end
                Δ = remainder(tmres) + remainder(tmn) * remainder(tmg)
            end

            return $TM(tmres.pol, Δ, x0, I)
        end
    end
end

function rpa(g::Function, tmf::TaylorModel1{TaylorN{T}}) where {T}
    _order = get_order(tmf)
    f_pol = polynomial(tmf)
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    I = domain(tmf)
    range_tmf = f_pol(I-x0) + Δf # TaylorN{Interval{...}}
    N = get_numvars()
    symIbox = IntervalBox(-1 .. 1, Val(N))
    interval_range_tmf = range_tmf(symIbox)
    tmg = _rpa(TaylorModel1, g, f_pol0, interval_range_tmf, _order)
    tm1 = tmf - f_pol0
    tmres = tmg(tm1)
    Δ = remainder(tmres) + remainder(tmg)
    return TaylorModel1(tmres.pol, Δ, x0, I)
end

function rpa(g::Function, tmf::RTaylorModel1{TaylorN{T}}) where {T}
    _order = get_order(tmf)
    f_pol = polynomial(tmf)
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    I = domain(tmf)
    range_tmf = f_pol(I-x0) + Δf # TaylorN{Interval{...}}
    N = get_numvars()
    symIbox = IntervalBox(-1 .. 1, Val(N))
    interval_range_tmf = range_tmf(symIbox)
    tmg = _rpa(RTaylorModel1, g, f_pol0, interval_range_tmf, _order)
    tm1 = tmf - f_pol0
    tmres = tmg(tm1)
    tmn = RTaylorModel1(Taylor1(copy(tm1.pol.coeffs)), tm1.rem, x0, I)
    for i = 1:_order
        tmn = tmn * tmf
    end
    Δ = remainder(tmres) + remainder(tmn) * remainder(tmg)
    return RTaylorModel1(tmres.pol, Δ, x0, I)
end

function rpa(g::Function, tmf::TaylorModel1{TaylorModelN{N,S,T},T}) where {N, T<:Real, S<:Real}
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

    f_pol = polynomial(tmf)
    f_pol0 = constant_term(f_pol)
    Δf = remainder(tmf)
    x0 = tmf.x0
    I = domain(tmf)

    # Range of tmf including remainder (Δf)
    # range_tmf = bound_taylor1(f_pol, I-x0) + Δf
    range_tmf = f_pol(I-x0) + Δf  # TaylorModelN
    interval_range_tmf = range_tmf(range_tmf.dom-range_tmf.x0)

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    tmg = _rpa(TaylorModel1, g, f_pol0, interval_range_tmf, _order)

    # Use original independent variable
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
    I = domain(tmf)

    # Range of tmf including remainder (Δf)
    range_tmf = f_pol(I-x0) + Δf

    # Compute RPA for `g`, around constant_term(f_pol), over range_tmf
    # Note that tmg is a TaylorModel1 !!
    tmg = _rpa(TaylorModel1, g, f_pol0, range_tmf, _order)

    # Use original independent variable
    tm1 = tmf - f_pol0   # OVER-ESTIMATION; IMPROVE
    tmres = tmg( tm1 )

    # Final remainder
    Δ = remainder(tmres) + remainder(tmg)
    return TaylorModelN(tmres.pol, Δ, x0, I)
end


"""
    fp_rpa(tm::TaylorModel1{Interval{T},T})
    fp_rpa(tm::RTaylorModel1{Interval{T},T})

Convert a `tm` TaylorModel1 to a TaylorModel1 whose polynomial coefficients
are `Float64`. The accumulated error is added to the remainder. The
mid point of the expansion interval is preferentially rounded down
if it is not an exactly representable value.

""" fp_rpa

for TM in tupleTMs
    @eval begin
        fp_rpa(tm::$TM{T, T}) where {T} = tm

        function fp_rpa(tm::$TM{Interval{T},T}) where {T}
            fT = polynomial(tm)
            Δ = remainder(tm)
            x0 = tm.x0
            I = domain(tm)
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
            return $TM(t, Δ, x0, I)
        end
    end
end

fp_rpa(tm::TaylorModelN{N, T, T}) where {N, T} = tm

function fp_rpa(tm::TaylorModelN{N,Interval{T},T}) where {N,T}
    fT = polynomial(tm)
    Δ = remainder(tm)
    x0 = tm.x0
    I = domain(tm)
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


# Elementary functions
fnlist = (:inv, :sqrt, :exp, :log, :sin, :cos, :tan,
    :asin, :acos, :atan, :sinh, :cosh, :tanh)

for fn in fnlist
    for TM in tupleTMs
        @eval $fn(tm::$TM) = rpa($fn, tm)
    end
    @eval $fn(tm::TaylorModelN) = rpa($fn, tm)
end
