# evaluate.jl

# evaluate and _evaluate functions

for TM in tupleTMs
    # Evaluates the TM1{T,S} at a::T; the computation includes the remainder
    @eval function evaluate(tm::$TM{T,S}, a) where {T,S}
        _order = get_order(tm)

        if $(TM) == TM1AbsRem
            Δ = tm.rem
        else
            Δ = tm.rem * a^(_order + 1)
        end

        return tm.pol(a) + Δ
    end

    @eval (tm::$TM{T,S})(a) where {T,S} = evaluate(tm, a)


    # _evaluate corresponds to composition: substitute tmf into tmg
    # It **does not** include the remainder
    @eval function _evaluate(tmg::$TM, tmf::$TM)
        _order = get_order(tmf)
        @assert _order == get_order(tmg)

        tmres = $TM(zero(constant_term(tmg.pol)), _order, tmf.x0, tmf.iI)
        @inbounds for k = _order:-1:0
            tmres = tmres * tmf
            tmres = tmres + tmg.pol[k]
        end

        # Returned result does not include remainder of tmg
        return tmres
    end

    @eval (tm::$TM)(x::$TM) = _evaluate(tm, x)

end



# Substitute a TMNAbsRem into a TM1; it **does not** include the remainder
function _evaluate(tmg::TM1AbsRem, tmf::TMNAbsRem)
    _order = get_order(tmf)
    @assert _order == get_order(tmg)

    tmres = TMNAbsRem(zero(constant_term(tmg.pol)), _order, tmf.x0, tmf.iI)
    @inbounds for k = _order:-1:0
        tmres = tmres * tmf
        tmres = tmres + tmg.pol[k]
    end

    # Returned result does not include remainder of tmg
    return tmres
end

(tm::TM1AbsRem)(x::TMNAbsRem) = _evaluate(tm, x)



# Evaluates the TMN on an interval, or array with proper dimension;
# the computation includes the remainder
function evaluate(tm::TMNAbsRem{N,T,S}, a::IntervalBox{N,S}) where {N,T,S}
    _order = get_order(tm)

    Δ = tm.rem

    return tm.pol(a...) + Δ
end

(tm::TMNAbsRem{N,T,S})(a::IntervalBox{N,S}) where {N,T,S} = evaluate(tm, a)

function evaluate(tm::TMNAbsRem{N,T,S}, a::Array{R,1}) where {N,T,S,R}
    _order = get_order(tm)

    Δ = tm.rem

    return tm.pol(a...) + Δ
end

(tm::TMNAbsRem{N,T,S})(a::Array{R,1}) where {N,T,S,R} = evaluate(tm, a)
