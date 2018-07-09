# promotion.jl
#

# Promotion
for TM in tupleTMs
    # @eval promote_type(::Type{$TM{T,S}}, ::Type{T}) where {T,S} = $TM{T,S}
    # @eval promote_type(::Type{$TM{T,S}}, ::Type{S}) where {T,S} = $TM{T,S}

    @eval function promote(a::$TM{T,S}, b::T) where {T,S}
        (a, $(TM)(Taylor1([b], get_order(a)), zero(a.rem), a.x0, a.iI))
    end
    @eval function promote(a::$TM{T,S}, b::S) where {T,S}
        (a, $TM(Taylor1([convert(T, b)], get_order(a)), zero(a.rem), a.x0, a.iI))
    end
    @eval function promote(a::$TM{TMNAbsRem{N,T,S},S}, b::T) where {N,T,S}
        tmN = TMNAbsRem(b, get_order(a.pol[0]), a.pol[0].x0, a.pol[0].iI)
        (a, $TM(Taylor1([tmN], get_order(a)), zero(a.rem), a.x0, a.iI))
    end
end

# promote_type(::Type{TMNAbsRem{N,T,S}}, ::Type{T}) where {N,T,S} = TMNAbsRem{N,T,S}
# promote_type(::Type{TMNAbsRem{N,T,S}}, ::Type{S}) where {N,T,S} = TMNAbsRem{N,T,S}

function promote(a::TMNAbsRem{N,T,S}, b::T) where {N,T,S}
    (a, TMNAbsRem(b, get_order(a), a.x0, a.iI))
end
function promote(a::TMNAbsRem{N,T,S}, b::S) where {N,T,S}
    (a, TMNAbsRem(b, get_order(a), a.x0, a.iI))
end
