# promotion.jl
#

# # Promotion
function promote(a::TaylorModelN{N,T,S}, b::R) where {N, T<:Real, S<:Real, R<:Real}
    orderTMN = get_order(a[0])
    apol, bb = promote(a.pol, b)
    return (TaylorModelN(apol, a.rem, a.x0, a.dom), TaylorModelN(bb, zero(a.rem), a.x0, a.dom))
end
promote(b::R, a::TaylorModelN{N,T,S}) where {N, T<:Real, S<:Real, R<:Real} =
    reverse(promote(a,b))

#
function promote(a::TaylorModelN{N,T,S}, b::TaylorN{R}) where {N, T, S, R}
    RR = promote_type(T,R)
    aa = TaylorModelN(convert(TaylorN{RR},a.pol), a.rem, a.x0, a.dom)
    bb = TaylorModelN(convert(TaylorN{RR},b), 0..0, a.x0, a.dom)
    (aa, bb)
end
promote(b::TaylorN{R}, a::TaylorModelN{N,T,S}) where {N, T, S, R} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,Interval{T},S}}, b::Taylor1{Interval{T}}) where
        {N, T<:Real, S<:Real}

    orderTMN = get_order(a[0])
    bTN = Array{TaylorModelN{N,Interval{T},S}}(undef, get_order(a)+1)
    unoTM = TaylorModelN(Interval(T(1), T(1)), orderTMN, a[0].x0, a[0].dom)
    @inbounds for ord = 1:length(bTN)
        bTN[ord] = b.coeffs[ord] * unoTM
    end
    (a, Taylor1(bTN))
end
promote(b::Taylor1{Interval{T}}, a::Taylor1{TaylorModelN{N,Interval{T},S}}) where
    {N, T<:Real, S<:Real} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,Interval{T},S}}, b::T) where
        {N, T<:Real, S<:Real}

    orderTMN = get_order(a[0])
    bTN = Array{TaylorModelN{N,Interval{T},S}}(undef, get_order(a)+1)
    unoTM = TaylorModelN(Interval(one(T), one(T)), orderTMN, a[0].x0, a[0].dom)
    bTN[1] = b * unoTM
    @inbounds for ord = 2:length(bTN)
        bTN[ord] = zero(b)*unoTM
    end
    (a, Taylor1(bTN))
end
promote(b::T, a::Taylor1{TaylorModelN{N,Interval{T},S}}) where
        {N, T<:Real, S<:Real} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,T,S}}, b::R) where
        {N, T<:Real, S<:Real, R<:Real}

    orderTMN = get_order(a[0])
    TT = promote_type(T,R)
    aTN = Array{TaylorModelN{N,TT,S}}(undef, get_order(a)+1)
    bTN = Array{TaylorModelN{N,TT,S}}(undef, get_order(a)+1)
    unoTT = one(TT)
    zeroTT = zero(TT)
    @inbounds for ord = 1:length(aTN)
        aTN[ord] = unoTT * a.coeffs[ord]
        bTN[ord] = zeroTT * a.coeffs[ord]
    end
    bTN[1] += b
    (Taylor1(aTN), Taylor1(bTN))
end
promote(b::R, a::Taylor1{TaylorModelN{N,T,S}}) where
        {N, T<:Real, S<:Real, R<:Real} = reverse( promote(a, b) )

# function promote(a::Taylor1{TaylorModelN{N,T,S}}, b::Taylor1{S}) where {N,T<:Real,S<:Real}
#     (a, Taylor1(TaylorModelN(interval(b), get_order(a), a.x0, a.dom), b.order))
# end
# promote(b::Taylor1{Interval{T}}, a::Taylor1{TaylorModelN{N,Interval{T},S}}) where
#     {N, T<:Real, S<:Real} = promote(a, b)

for TM in tupleTMs
    @eval promote(a::$TM{T,S}, b::T) where {T,S} =
        (a, $(TM)(Taylor1([b], get_order(a)), zero(a.rem), a.x0, a.dom))
    @eval promote(b::T, a::$TM{T,S}) where {T,S} = reverse( promote(a,b) )
    #
    @eval promote(a::$TM{T,S}, b::S) where {T,S} =
        (a, $TM(Taylor1([convert(T, b)], get_order(a)), zero(a.rem), a.x0, a.dom))
    @eval promote(b::S, a::$TM{T,S}) where {T,S} = reverse( promote(a,b) )
    #
    @eval promote(a::$TM{T,S}, b::R) where {T,S,R} = promote(a, convert(S, b))
    @eval promote(b::R, a::$TM{T,S}) where {T,S,R} = reverse( promote(a,b) )
    #
    @eval function promote(a::$TM{TaylorModelN{N,T,S},S}, b::T) where {N,T,S}
        tmN = TaylorModelN(b, get_order(a.pol[0]), a.pol[0].x0, a.pol[0].dom)
        (a, $TM(Taylor1([tmN], get_order(a)), zero(a.rem), a.x0, a.dom))
    end
    @eval promote(b::T, a::$TM{TaylorModelN{N,T,S},S}) where {N,T,S} = reverse( promote(a,b) )
end
