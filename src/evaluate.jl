# evaluate.jl

# evaluate and _evaluate functions

for TM in tupleTMs
    # Evaluates the TM1{T,S} at a::T; the computation includes the remainder
    @eval function evaluate(tm::$TM{T,S}, a) where {T,S}
        @assert iscontained(a, tm)

        _order = get_order(tm)

        if $(TM) == TaylorModel1
            Δ = remainder(tm)
        else
            Δ = remainder(tm) * a^(_order + 1)
        end

        return tm.pol(a) + Δ
    end

    @eval (tm::$TM{T,S})(a) where {T,S} = evaluate(tm, a)

    @eval evaluate(tm::Vector{$TM{T,S}}, a) where {T,S} = evaluate.(tm, a)

    # Evaluate the $TM{TaylorN} by assuming the TaylorN vars are properly symmetrized,
    # and thus `a` is contained in the corresponding [-1,1] box; this is not checked.
    @eval evaluate(tm::$TM{TaylorN{T},S}, a::IntervalBox) where {T<:NumberNotSeries,S} =
        evaluate(tm, Vector(a.v))
    @eval function evaluate(tm::$TM{TaylorN{T},S}, a::AbstractVector{R}) where {T<:NumberNotSeries,S,R}
        @assert length(a) == get_numvars()
        pol = tm.pol(a)
        return $TM(pol, tm.rem, one(pol[0])*tm.x0, tm.dom)
    end

    # _evaluate corresponds to composition: substitute tmf into tmg
    # It **does not** include the remainder
    @eval function _evaluate(tmg::$TM, tmf::$TM)
        _order = get_order(tmf)
        @assert _order == get_order(tmg)

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
    _order = get_order(tmf)
    @assert _order == get_order(tmg)

    tmres = TaylorModelN(zero(constant_term(tmg.pol)), _order,
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
function evaluate(tm::TaylorModelN{N,T,S}, a::IntervalBox{N,S}) where {N,T,S}
    @assert iscontained(a, tm)
    _order = get_order(tm)

    Δ = remainder(tm)

    return tm.pol(a) + Δ
end

(tm::TaylorModelN{N,T,S})(a::IntervalBox{N,S}) where {N,T,S} = evaluate(tm, a)

function evaluate(tm::TaylorModelN{N,T,S}, a::AbstractVector{R}) where {N,T,S,R}
    @assert iscontained(a, tm)
    _order = get_order(tm)

    Δ = remainder(tm)

    return tm.pol(a) + Δ
end

(tm::TaylorModelN{N,T,S})(a::AbstractVector{R}) where {N,T,S,R} = evaluate(tm, a)

evaluate(tm::Vector{TaylorModelN{N,T,S}}, a::IntervalBox{N,S}) where {N,T,S} =
    IntervalBox( [ tm[i](a) for i in eachindex(tm) ] )


function evaluate(a::Taylor1{TaylorModelN{N,T,S}}, dx::T) where {N, T<:Real, S<:Real}
    @inbounds suma = a[end]*one(dx)
    @inbounds for k in get_order(a)-1:-1:0
        suma = suma*dx + a[k]
    end
    return suma
end
