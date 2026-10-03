# evaluate.jl

# evaluate and _evaluate functions


# Interal: Split an accumulated Interval{S} into a stored coefficient of type T,
# T being S or Interval{S}, and the residual interval (such that
# acc ⊆ coef + residual) that must go into the remainder.
_split(::Type{S}, acc::Interval{S}) where {S<:AbstractFloat} =
    (m = mid(acc); (m, acc - interval(m)))
_split(::Type{Interval{S}}, acc::Interval{S}) where {S<:AbstractFloat} =
    (acc, zero(Interval{S}))


for TM in tupleTMs
    # Evaluates the $TM at an interval `a`; the computation includes the remainder
    @eval function evaluate(tm::$TM{T,S}, a::Interval) where {T,S}
        @assert iscontained(a, tm)
        if $(TM) == TaylorModel1
            Δ = remainder(tm)
        else
            _order = TS.order(tm)
            Δ = remainder(tm) * Base.literal_pow(^, a, Val(_order+1))
        end
        return tm.pol(a) + Δ
    end

    # Real point: enclose it first, so that the rounding of the evaluation of
    # the polynomial and remainder are enclosed.
    @eval evaluate(tm::$TM{T,S}, a::Real) where {T,S} = evaluate(tm, interval(S, a))

    # Other arguments: unchanged behaviour
    @eval function evaluate(tm::$TM{T,S}, a) where {T,S}
        @assert iscontained(a, tm)
        if $(TM) == TaylorModel1
            Δ = remainder(tm)
        else
            _order = TS.order(tm)
            Δ = remainder(tm) * Base.literal_pow(^, a, Val(_order+1))
        end
        return tm.pol(a) + Δ
    end

    @eval (tm::$TM{T,S})(a) where {T,S} = evaluate(tm, a)

    @eval evaluate(tm::Vector{$TM{T,S}}, a) where {T,S} = evaluate.(tm, a)

    # Evaluate the $TM{TaylorN} by assuming the TaylorN vars are properly symmetrized,
    # and thus `a` is contained in the corresponding [-1,1] box; this is not checked.
    @eval function evaluate(tm::$TM{TaylorN{T},S}, a::AbstractVector{R}) where
            {T<:NumberNotSeries,S,R}
        @assert length(a) == get_numvars()
        pol = tm.pol(a)
        return $(Symbol(:unsafe_, TM))(pol, tm.rem, tm.x0, tm.dom)
    end

    # _evaluate corresponds to composition: substitute tmf into tmg
    # It **does not** include the remainder
    @eval function _evaluate(tmg::$TM, tmf::$TM)
        _order = TS.order(tmf)
        @assert _order == TS.order(tmg)

        tmres = zero(tmg.pol[_order]) * tmf
        tmres = tmres + tmg.pol[_order]
        @inbounds for k = _order-1:-1:0
            tmres = tmres * tmf
            tmres = tmres + tmg.pol[k]
        end

        # Returned result does not include remainder of tmg
        return tmres
    end

    @eval (tm::$TM)(x::$TM) = _evaluate(tm, x)
end



# Substitute a TaylorModelN into a TM1; it **does not** include the remainder
function _evaluate(tmg::TaylorModel1{T,S}, tmf::TaylorModelN{N,T,S}) where{N,T,S}
    _order = TS.order(tmf)
    @assert _order == TS.order(tmg)
    tmres = TaylorModelN(space(tmf), zero(constant_term(tmg.pol)), _order,
        expansion_point(tmf), domain(tmf))
    @inbounds for k = _order:-1:0
        tmres = tmres * tmf
        tmres = tmres + tmg.pol[k]
    end
    # Returned result does not include remainder of tmg
    return tmres
end

(tm::TaylorModel1)(x::TaylorModelN) = _evaluate(tm, x)


# Evaluates the TMN on an interval, or array with proper dimension;
# the computation includes the remainder
function evaluate(tm::TaylorModelN{N,T,S}, a::AbstractVector{Interval{S}}) where {N,T,S}
    @assert iscontained(a, tm)
    Δ = remainder(tm)
    return tm.pol(a) + Δ
end

(tm::TaylorModelN{N,T,S})(a::AbstractVector{Interval{S}}) where {N,T,S} = evaluate(tm, a)

function evaluate(tm::TaylorModelN{N,T,S}, a::AbstractVector{R}) where {N,T,S,R<:Real}
    @assert iscontained(a, tm)
    return evaluate(tm, interval.(Ref(S), a))
end

(tm::TaylorModelN{N,T,S})(a::AbstractVector{R}) where {N,T,S,R} = evaluate(tm, a)

evaluate(tm::Vector{TaylorModelN{N,T,S}}, a::AbstractVector{Interval{S}}) where {N,T,S} =
    interval.( tm[i](a) for i in eachindex(tm) )


for D in (:S, :(Interval{S}))
    @eval begin
        evaluate(a::Taylor1{TaylorModelN{N,T,S}}, dx::$D) where {N,T,S<:AbstractFloat} =
            _horner!(deepcopy(a[0]), a, dx, _aux_horner(a[0]))
        evaluate(tm::TaylorModel1{TaylorModelN{N,T,S},S}, dx::$D) where {N,T,S<:AbstractFloat} =
            __evaluate!(deepcopy(tm[0]), tm, dx)
    end
end
function evaluate(a::Taylor1{TaylorModelN{N,T,S}}, v::AbstractVector{R}) where {N,T,S,R}
    suma = Taylor1(zero(a[0])(v), TS.order(a))
    for k in eachindex(a)
        suma[k] = a[k](v)
    end
    return suma
end
function evaluate(tm::TaylorModel1{TaylorModelN{N,T,S},S}, v::AbstractVector) where {N,T,S}
    suma = Taylor1(zero(tm[0])(v), TS.order(tm))
    for k in eachindex(suma)
        suma[k] = tm[k](v)
    end
    return unsafe_TaylorModel1(suma, remainder(tm), expansion_point(tm), domain(tm))
end

function _evaluate(tm::TaylorModelN{N,T,S},
        dx::AbstractVector{TaylorModelN{N,T,S}}) where {N,T,S}
    @assert N == length(dx)
    res = evaluate(polynomial(tm), dx)
    res.rem += tm.rem
    return res
end


"""
    __evaluate!(tmn, tm::TaylorModel1{TaylorModelN}, dx[, aux]) -> tmn

In-place rigorous evaluation: `tmn ⊇ tm(dx)` (all remainders included), `dx`
an `S` or an `Interval{S}` within the centred domain of `tm`. `aux` is an
optional `TaylorN{Interval{S}}` workspace (see `_aux_horner`); without it, one
is allocated.
"""
function __evaluate!(tmn::TaylorModelN{N,T,S},
        tm::TaylorModel1{TaylorModelN{N,T,S},S}, dx::Union{S,Interval{S}},
        aux::TaylorN{Interval{S}}) where {N,T,S}
    @assert issubset_interval(interval(dx), centered_dom(tm))
    _horner!(tmn, polynomial(tm), dx, aux)
    tmn.rem += remainder(tm)
    return tmn
end
__evaluate!(tmn::TaylorModelN{N,T,S},
        tm::TaylorModel1{TaylorModelN{N,T,S},S},
        dx::Union{S,Interval{S}}) where {N,T,S} =
    __evaluate!(tmn, tm, dx, _aux_horner(tm[0]))



"""
    _horner!(tmn, a::Taylor1{TaylorModelN}, dx, aux) -> tmn

Rigorous `tmn ⊇ Σ_k a[k] dx^k`, including the remainders of the TMN
coefficients, for `dx` an `S` or an `Interval{S}` (a displacement of the
Taylor1 variable). Horner's scheme runs coefficientwise in `Interval{S}`,
accumulated in the workspace `aux::TaylorN{Interval{S}}` (same structure as the
coefficients of `a`); coefficients are stored in `tmn` with `_split`, the
residuals are kept in `aux` and bounded over `centered_dom(tmn)` (as in
`fp_rpa`). No allocations besides the evaluation of the residual polynomial.
`tmn` provides expansion point and domain; it must not alias any `a[k]`, and
`aux` must not alias anything.
"""
function _horner!(tmn::TaylorModelN{N,T,S}, a::Taylor1{TaylorModelN{N,T,S}},
        dx::Union{S,Interval{S}}, aux::TaylorN{Interval{S}}) where {N,T,S}
    zI = zero(Interval{S})
    dxI = interval(dx)
    @inbounds for q in eachindex(aux.coeffs)
        for h in eachindex(aux.coeffs[q].coeffs)
            aux.coeffs[q].coeffs[h] = zI
        end
    end
    rem = zI
    @inbounds for k in reverse(eachindex(a))
        rem = rem * dxI + remainder(a[k])
        ck = a[k].pol
        for q in eachindex(aux.coeffs)
            for h in eachindex(aux.coeffs[q].coeffs)
                aux.coeffs[q].coeffs[h] = aux.coeffs[q].coeffs[h] * dxI +
                    interval(ck.coeffs[q].coeffs[h])
            end
        end
    end
    @inbounds for q in eachindex(aux.coeffs), h in eachindex(aux.coeffs[q].coeffs)
        tmn.pol.coeffs[q].coeffs[h], aux.coeffs[q].coeffs[h] =
            _split(T, aux.coeffs[q].coeffs[h]) # aux ← residuals
    end
    tmn.rem = rem + aux(centered_dom(tmn))
    return tmn
end

# Interval workspace with the structure of the coefficients of `a`
_aux_horner(a::TaylorModelN{N,T,S}) where {N,T,S} =
    interval(zero(S)) * zero(polynomial(a))     # TaylorN{Interval{S}}
