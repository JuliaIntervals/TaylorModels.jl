# promotion.jl
#

# # Promotion
function promote(a::TaylorModelN{T,S}, b::R) where {T<:Real, S<:Real, R<:Real}
    apol, bb = promote(a.pol, b)
    a_rem = remainder(a)
    return (TaylorModelN(apol, a_rem, expansion_point(a), domain(a)),
        TaylorModelN(bb, zero(a_rem), expansion_point(a), domain(a)))
end
# promote(b::R, a::TaylorModelN{T,S}) where { T<:Real, S<:Real, R<:Real} =
#     reverse(promote(a,b))

#
function promote(a::TaylorModelN{N,T,S}, b::TaylorN{R}) where {N, T, S, R}
    RR = promote_type(T,R)
    aa = TaylorModelN(convert(TaylorN{RR},a.pol), remainder(a), expansion_point(a),
        domain(a))
    bb = TaylorModelN(convert(TaylorN{RR},b), 0..0, expansion_point(a), domain(a))
    return (aa, bb)
end
# promote(b::TaylorN{R}, a::TaylorModelN{T,S}) where { T, S, R} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,Interval{T},S}}, b::Taylor1{Interval{T}}) where
        {N,T<:Real, S<:Real}
    orderTMN = get_order(a[0])
    bTN = Array{TaylorModelN{N,Interval{T},S}}(undef, get_order(a)+1)
    unoTM = TaylorModelN(interval(one(T)), orderTMN, expansion_point(a[0]), domain(a[0]))
    @inbounds for ord in eachindex(bTN)
        bTN[ord] = b.coeffs[ord] * unoTM
    end
    return (a, Taylor1(bTN))
end
# promote(b::Taylor1{Interval{T}}, a::Taylor1{TaylorModelN{Interval{T},S}}) where
#     { T<:Real, S<:Real} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,Interval{T},S}}, b::T) where
        {N, T<:Real, S<:Real}
    orderTMN = get_order(a[0])
    bTN = Array{TaylorModelN{N,Interval{T},S}}(undef, get_order(a)+1)
    unoTM = TaylorModelN(interval(one(T)), orderTMN, expansion_point(a[0]), domain(a[0]))
    bTN[1] = b * unoTM
    @inbounds for ord = 2:length(bTN)
        bTN[ord] = zero(b)*unoTM
    end
    return (a, Taylor1(bTN))
end
# promote(b::T, a::Taylor1{TaylorModelN{Interval{T},S}}) where
#         { T<:Real, S<:Real} = reverse( promote(a, b) )

function promote(a::Taylor1{TaylorModelN{N,T,S}}, b::R) where
        {N, T<:Real, S<:Real, R<:Real}
    # orderTMN = get_order(a[0])
    TT = promote_type(T,R)
    aTN = Array{TaylorModelN{TT,S}}(undef, get_order(a)+1)
    bTN = Array{TaylorModelN{TT,S}}(undef, get_order(a)+1)
    unoTT = one(TT)
    zeroTT = zero(TT)
    @inbounds for ord = 1:length(aTN)
        aTN[ord] = unoTT * a.coeffs[ord]
        bTN[ord] = zeroTT * a.coeffs[ord]
    end
    bTN[1] += b
    return (Taylor1(aTN), Taylor1(bTN))
end
# promote(b::R, a::Taylor1{TaylorModelN{T,S}}) where
#         { T<:Real, S<:Real, R<:Real} = reverse( promote(a, b) )


for TM in tupleTMs
    @eval promote(a::$TM{T,S}, b::T) where {T,S<:Real} =
        (a, $(TM)(Taylor1([b], get_order(a)), zero(remainder(a)), expansion_point(a), domain(a)))
    @eval promote(b::T, a::$TM{T,S}) where {T,S<:Real} = reverse( promote(a,b) )
    #
    @eval promote(a::$TM{T,S}, b::S) where {T,S} =
        (a, $TM(Taylor1([convert(T, b)], get_order(a)), zero(remainder(a)), expansion_point(a), domain(a)))
    # @eval promote(b::S, a::$TM{T,S}) where {T,S} = reverse( promote(a,b) )
    #
    # @eval promote(a::$TM{T,S}, b::R) where {T,S<:Real,R} = promote(a, convert(S, b))
    # @eval promote(b::R, a::$TM{T,S}) where {T,S,R} = reverse( promote(a,b) )
    #
    @eval function promote(a::$TM{TaylorModelN{N,T,S},S}, b::T) where {N,T,S}
        a_pol0 = a.pol[0]
        tmN = TaylorModelN(b, get_order(a_pol0), expansion_point(a_pol0), domain(a_pol0))
        return (a, $TM(Taylor1([tmN], get_order(a)), zero(remainder(a)), expansion_point(a), domain(a)))
    end
    # @eval promote(b::T, a::$TM{TaylorModelN{N,T,S},S}) where {N,T,S} = reverse( promote(a,b) )
end


function Base.convert(::Type{TaylorModel1{TaylorN{T}, S}}, tm::TaylorModel1{TaylorModelN{N,T,S}, S}) where {N,T,S}
    order = get_order(tm)
    pol = Taylor1(zero(polynomial(tm[0])), order)
    for k in eachindex(polynomial(tm))
        pol[k] = polynomial(tm[k])
    end
    rem = total_remainder(tm)
    return TaylorModel1(pol, rem, expansion_point(tm), domain(tm))
end

function Base.convert(::Type{TaylorModel1{TaylorModelN{N,T,S}, S}}, tm::TaylorModel1{TaylorN{T}, S}) where {N,T,S}
    z = interval(zero(S))
    symIbox = symmetric_box(S)
    zSB = interval.(mid.(symIbox))
    order = get_order(tm)
    pol = Taylor1(TaylorModelN.(tm[:], z, Ref(zSB), Ref(symIbox)), order)
    return TaylorModel1(pol, remainder(tm), expansion_point(tm), domain(tm))
end

"""
    total_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S})

Computes de total reminder of a `TaylorModel1{TaylorModelN}` by
computing the polynomial associated to the `TaylorModelN` reminders,
evaluated in the centered domain, and adding the remainder of the
`TaylorModel1`.
"""
function total_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S}) where {N,T,S}
    order = get_order(tm)
    polI = Taylor1(zero(remainder(tm[0])), order)
    for k in eachindex(polynomial(tm))
        polI[k] = remainder(tm[k])
    end
    return polI(centered_dom(tm))+remainder(tm)
end

"""
    shift_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S})

Returns a `TaylorModel1{TaylorModelN}` with null remainder for
the `TaylorModelN` part, and the total remainder in the `TaylorModel1`.
"""
function shift_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S}) where {N,T,S}
    rem = total_remainder(tm)
    z = interval(0.0)
    for k in eachindex(polynomial(tm))
        tm[k] = TaylorModelN(tm[k], z)
    end
    return TaylorModel1(tm, rem)
end
