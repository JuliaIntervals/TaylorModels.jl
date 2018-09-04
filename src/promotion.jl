# promotion.jl
#

# Promotion
promote(a::TaylorModelN{N,Interval{T},S}, b::S) where {N, T<:Real, S<:Real} =
    (a, TaylorModelN(IntervalArithmetic.atomic(Interval{T}, b), get_order(a), a.x0, a.I))
promote(b::S, a::TaylorModelN{N,Interval{T},S}) where {N, T<:Real, S<:Real} =
    reverse( promote(a, b) )
#
promote(a::TaylorModelN{N,Interval{T},S}, b::R) where {N, T<:Real, S<:Real, R<:Real} =
    promote(a, IntervalArithmetic.atomic(Interval{T}, b))
promote(b::R, a::TaylorModelN{N,Interval{T},S}) where {N, T<:Real, S<:Real, R<:Real} =
    reverse( promote(a, b) )
#
promote(a::TaylorModelN{N,Interval{T},S}, b::Interval{T}) where {N, T<:Real, S<:Real} =
    (a, TaylorModelN(b, get_order(a), a.x0, a.I))
promote(b::Interval{T}, a::TaylorModelN{N,Interval{T},S}) where {N, T<:Real, S<:Real} =
    reverse( promote(a, b) )
#
function promote(a::TaylorModelN{N,T,S}, b::TaylorN{R}) where {N, T, S, R}
    RR = promote_type(T,R)
    aa = TaylorModelN(convert(TaylorN{RR},a.pol), a.rem, a.x0, a.I)
    bb = TaylorModelN(convert(TaylorN{RR},b), 0..0, a.x0, a.I)
    (aa, bb)
end
promote(b::TaylorN{R}, a::TaylorModelN{N,T,S}) where {N, T, S, R} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,Interval{T},S}}, b::Taylor1{Interval{T}}) where
        {N, T<:Real, S<:Real}

    orderTMN = get_order(a[0])
    bTN = Array{TaylorModelN{N,Interval{T},S}}(get_order(a)+1)
    for ord = 1:length(bTN)
        bTN[ord] = b.coeffs[ord]*TaylorModelN(1..1, orderTMN, a[0].x0, a[0].I)
    end
    (a, Taylor1(bTN))
end
promote(b::Taylor1{Interval{T}}, a::Taylor1{TaylorModelN{N,Interval{T},S}}) where
    {N, T<:Real, S<:Real} = reverse( promote(a, b) )

# function promote(a::Taylor1{TaylorModelN{N,T,S}}, b::Taylor1{S}) where {N,T<:Real,S<:Real}
#     (a, Taylor1(TaylorModelN(interval(b), get_order(a), a.x0, a.I), b.order))
# end
# promote(b::Taylor1{Interval{T}}, a::Taylor1{TaylorModelN{N,Interval{T},S}}) where
#     {N, T<:Real, S<:Real} = promote(a, b)

for TM in tupleTMs
    @eval promote(a::$TM{T,S}, b::T) where {T,S} =
        (a, $(TM)(Taylor1([b], get_order(a)), zero(a.rem), a.x0, a.I))
    @eval promote(b::T, a::$TM{T,S}) where {T,S} = reverse( promote(a,b) )
    #
    @eval promote(a::$TM{T,S}, b::S) where {T,S} =
        (a, $TM(Taylor1([convert(T, b)], get_order(a)), zero(a.rem), a.x0, a.I))
    @eval promote(b::S, a::$TM{T,S}) where {T,S} = reverse( promote(a,b) )
    #
    @eval promote(a::$TM{T,S}, b::R) where {T,S,R} = promote(a, convert(S, b))
    @eval promote(b::R, a::$TM{T,S}) where {T,S,R} = reverse( promote(a,b) )
    #
    @eval function promote(a::$TM{TaylorModelN{N,T,S},S}, b::T) where {N,T,S}
        tmN = TaylorModelN(b, get_order(a.pol[0]), a.pol[0].x0, a.pol[0].I)
        (a, $TM(Taylor1([tmN], get_order(a)), zero(a.rem), a.x0, a.I))
    end
    @eval promote(b::T, a::$TM{TaylorModelN{N,T,S},S}) where {N,T,S} = reverse( promote(a,b) )
end
