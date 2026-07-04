# Mutating functions
#

function TS.identity!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        k::Int) where {T,S,N}
    TS.identity!(c.coeffs[k+1], a.coeffs[k+1])
    return nothing
end
function TS.identity!(c::TaylorModelN{N,T,S}, a::TaylorModelN{N,T,S}) where {T,S,N}
    for k in eachindex(c)
        TS.identity!(polynomial(c), polynomial(a), k)
    end
    c.rem = a.rem
    return nothing
end


function TS.zero!(a::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {T,S,N}
    TS.zero!(a.coeffs[k+1])
    return nothing
end
function TS.zero!(a::TaylorModelN)
    for k in eachindex(a)
        TS.zero!(polynomial(a), k)
    end
    a.rem = zero(a.rem)
    return nothing
end
function TS.zero!(a::TaylorModelN, k::Int)
    TS.zero!(polynomial(a), k)
    return nothing
end

function TS.one!(a::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {T,S,N}
    if k == 0
        TS.one!(a.coeffs[1])
    else
        TS.zero!(a.coeffs[k+1])
    end
    return nothing
end
function TS.one!(a::TaylorModelN)
    for k in eachindex(polynomial(a))
        TS.zero!(polynomial(a), k)
    end
    TS.one!(polynomial(a), 0)
    a.rem = zero(a.rem)
    return nothing
end
function TS.one!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        k::Int) where {T,S,N}
    if k == 0
        c.coeffs[1] = one(a.coeffs[1])
    else
        TS.zero!(c.coeffs[k+1])
    end
    return nothing
end

#
function TS.add!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        k::Int) where {N,T,S}
    TS.identity(c.coeffs[k+1], a.coeffs[k+1])
    return nothing
end
# function TS.add!(c::TaylorModelN{N,T,S}, a::TaylorModelN{N,T,S}) where {N,T,S}
#     for k in eachindex(c)
#         TS.add!(polynomial(c), polynomial(a), k)
#     end
#     c.rem = a.rem
#     return nothing
# end
function TS.add!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    c.coeffs[k+1] = a.coeffs[k+1] + b.coeffs[k+1]
    return nothing
end
# function TS.add!(c::TaylorModelN{N,T,S}, a::TaylorModelN{N,T,S}) where {N,T,S}
#     for k in eachindex(c)
#         TS.add!(polynomial(c), polynomial(a), polynomial(b), k)
#     end
#     c.rem = a.rem + b.rem
#     return nothing
# end

function TS.subst!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        k::Int) where {N,T,S}
    TS.identity!(c.coeffs[k+1], - a.coeffs[k+1])
    return nothing
end
# function TS.subst!(c::TaylorModelN{N,T,S}, a::TaylorModelN{N,T,S}) where {N,T,S}
#     for k in eachindex(c)
#         TS.subst!(polynomial(c), polynomial(a), k)
#     end
#     c.rem = -a.rem
#     return nothing
# end
function TS.subst!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    c.coeffs[k+1] = a.coeffs[k+1] - b.coeffs[k+1]
    return nothing
end
# function TS.subst!(c::TaylorModelN{N,T,S}, a::TaylorModelN{N,T,S},
#         b::Taylor1{TaylorModelN{N,T,S}}) where {N,T,S}
#     for k in eachindex(c)
#         TS.subst!(polynomial(c), polynomial(a), polynomial(b), k)
#     end
#     c.rem = a.rem - b.rem
#     return nothing
# end

function TS.mul!(res::Taylor1{TaylorModelN{N,T,S}}, a::NumberNotSeries,
        b::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    res.coeffs[k+1] = a * b.coeffs[k+1]
    return nothing
end

TS.mul!(res::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
    b::NumberNotSeries, k::Int) where {N,T,S} = TS.mul!(res, b, a, k)

function TS.mul!(res::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::Taylor1{TaylorModelN{N,T,S}}, ordT::Int) where {N,T,S}
    # Sanity
    TS.zero!(res, ordT)
    res_k = res.coeffs[ordT+1]
    for k in 0:ordT
        res_k += a.coeffs[k+1] * b.coeffs[ordT+1-k]
    end
    res[ordT] = res_k
    return nothing
end


function TS.div!(res::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::NumberNotSeries, k::Int) where {N,T,S}
    res.coeffs[k+1] = a.coeffs[k+1] / b
    return nothing
end

function TS.div!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    # order and coefficient of first factorized term
    # ordfact, cdivfact = divfactorization(a, b)
    anz = findfirst(a)
    bnz = findfirst(b)
    anz = anz ≥ 0 ? anz : TS.order(a)
    bnz = bnz ≥ 0 ? bnz : TS.order(a)
    ordfact = min(anz, bnz)
    # Is the polynomial factorizable?
    iszero(b[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )
    TS.zero!(c, k)
    if k == 0
        c.coeffs[1] = a.coeffs[ordfact+1] / b.coeffs[ordfact+1]
        # TS.div!(c[0], a[ordfact], b[ordfact])
        return nothing
    end
    imin = max(0, k+ordfact-TS.order(b))
    tmp = c.coeffs[imin+1] * b.coeffs[k+1+ordfact-imin]
    # TS.mul!(c[k], c[imin], b[k+ordfact-imin])
    for i = imin+1:k-1
        tmp += c.coeffs[i+1] * b.coeffs[k+1+ordfact-i]
        # TS.mul!(c[k], c[i], b[k+ordfact-i])
    end
    if k+ordfact ≤ TS.order(b)
        tmp = (a.coeffs[k+1+ordfact]-tmp)
        # for l in eachindex(c[k])
        #     TS.subst!(c[k], a[k+ordfact], c[k], l)
        # end
        c.coeffs[k+1] = tmp / b.coeffs[ordfact+1]
        # TS.div!(c[k], b[ordfact])
    else
        # c[k] = (-c[k]) / b[ordfact]
        c.coeffs[k+1] = tmp
        TS.div_scalar!(c.coeffs[k+1], -1, b.coeffs[ordfact+1])
    end
    return nothing
end

#
function TS.pow!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        tmp::Taylor1{TaylorModelN{N,T,S}}, r::U, k::Int) where {N,T,S,U}
    (r == 0) && return TS.one!(c, a, k)
    (r == 1) && return TS.identity!(c, a, k)
    (r == 2) && return TS.sqr!(c, a, constant_term(tmp), k)
    (r == 0.5) && return TS.sqrt!(c, a, tmp, k)
    # Sanity
    TS.zero!(c, k)
    # First non-zero coefficient
    l0 = findfirst(a)
    l0 < 0 && return nothing
    # # Index of first non-zero coefficient of the result; must be integer
    # !isinteger(r*l0) && throw(DomainError(a,
    #     """The 0-th order Taylor1 coefficient must be non-zero
    #     to raise the Taylor1 polynomial to a non-integer exponent."""))
    # lnull = trunc(Int, r*l0 )
    # kprime = k-lnull
    # (kprime < 0 || lnull > TS.order(a)) && return nothing
    # # Relevant for positive integer r, to avoid round-off errors
    # isinteger(r) && r > 0 && (k > r*findlast(a)) && return nothing
    # # First non-zero coeff
    # if k == lnull
    #     @inbounds c[k] = (a[l0])^float(r)
    #     return nothing
    # end
    # # The recursion formula
    # if l0+kprime ≤ TS.order(a)
    #     @inbounds c[k] = r * kprime * c[lnull] * a[l0+kprime]
    # end
    # for i = 1:k-lnull-1
    #     ((i+lnull) > TS.order(a) || (l0+kprime-i > TS.order(a))) && continue
    #     aaux = r*(kprime-i) - i
    #     @inbounds c[k] += aaux * c[i+lnull] * a[l0+kprime-i]
    # end
    # @inbounds c[k] = c[k] / (kprime * a[l0])
    # return nothing
    # #
    if k == l0
        c[0] = ( a[l0] )^r
        return nothing
    end
    tmp[k-l0] = zero(c[k-l0])
    for i = 1:k-l0-1
        aaux = r*(k-i) - i
        tmp[k-l0] += aaux * a[k-i] * c[i]
    end
    aaux = k - l0*(r+1)
    c[k-l0] = tmp[k-l0] / (aaux * a[l0])
    return nothing
end

function TS.sqr!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}, tmp::TaylorModelN{N,T,S},
        k::Int) where {N,T,S}
    if k == 0
        c[0] = a[0]^2
        return nothing
    end
    kodd = k%2
    kend = (k - 2 + kodd) >> 1
    TS.zero!(tmp.pol)
    for i = 0:kend
        tmp += a[i] * a[k-i]
    end
    c[k] = 2 * tmp
    kodd == 1 && return nothing
    c[k] += a[k>>1]^2
    return nothing
end

function TS.sqr_orderzero!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}) where {N,T,S}
    c[0] = a[0]^2
    return nothing
end

function TS.sqrt!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        ::Taylor1{TaylorModelN{N,T,S}}, k::Int, k0::Int=0) where {N,T,S}
    if k == k0
        c[k] = sqrt(a[2*k0])
        return nothing
    end
    kodd = (k - k0)%2
    kend = (k - k0 - 2 + kodd) >> 1
    tmp = zero(c[k])
    for i = k0+1:k0+kend
        tmp += c[i] * c[k+k0-i]
    end
    tmp = a[k+k0] - 2*tmp
    if kodd == 0
        tmp = tmp - (c[kend+k0+1])^2
    end
    c[k] = tmp / (2*c[k0])
    return nothing
end

#
function TS.exp!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        c[0] = exp(constant_term(a))
        return nothing
    end
    tmp = zero(c[k])
    for i = 0:k-1
        tmp += (k-i) * a[k-i] * c[i]
    end
    c[k] = tmp / k
    return nothing
end

function TS.log!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        c[0] = log(constant_term(a))
        return nothing
    end
    tmp = zero(c[k])
    for i = 1:k-1
        tmp += (k-i) * a[i] * c[k-i]
    end
    c[k] = (a[k] - tmp/k) / constant_term(a)
    return nothing
end

function TS.sincos!(s::Taylor1{TaylorModelN{N,T,S}},
        c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        a0 = constant_term(a)
        s[0], c[0] = sin( a0 ), cos( a0 )
        return nothing
    end
    tmps = zero(s[k])
    tmpc = zero(c[k])
    for i = 1:k
        x = i * a[i]
        tmps += x * c[k-i]
        tmpc -= x * s[k-i]
    end
    s[k] = tmps / k
    c[k] = tmpc / k
    return nothing
end

function TS.tan!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        c2::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        aux = tan( constant_term(a) )
        c[0] = aux
        c2[0] = aux^2
        return nothing
    end
    tmp = zero(c[k])
    for i = 0:k-1
        tmp += (k-i) * a[k-i] * c2[i]
    end
    c[k] = a[k] + tmp/k
    # c2 = c^2
    TS.sqr!(c2, c, tmp, k)
    return nothing
end

function TS.asin!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        r::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        a0 = constant_term(a)
        c[0] = asin( a0 )
        r[0] = sqrt( 1 - a0^2)
        return nothing
    end
    tmp = zero(c[k])
    for i in 1:k-1
        tmp += (k-i) * r[i] * c[k-i]
    end
    # r = sqrt(1-a^2)
    TS.sqrt!(r, 1-a^2, k)
    c[k] = (a[k] - tmp/k) / constant_term(r)
    return nothing
end

function TS.acos!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        r::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        a0 = constant_term(a)
        c[0] = acos( a0 )
        r[0] = sqrt( 1 - a0^2)
        return nothing
    end
    TS.asin!(c, -a, r, k)
    return nothing
end

function TS.atan!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        r::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        a0 = constant_term(a)
        c[0] = atan( a0 )
        r[0] = 1 + a0^2
        return nothing
    end
    tmp = zero(c[k])
    for i in 1:k-1
        tmp += (k-i) * r[i] * c[k-i]
    end
    c[k] = (a[k] - tmp/k) / constant_term(r)
    # r = a^2
    TS.sqr!(r, a, tmp, k)
    return nothing
end

function TS.sinhcosh!(s::Taylor1{TaylorModelN{N,T,S}},
        c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        s[0] = sinh( constant_term(a) )
        c[0] = cosh( constant_term(a) )
        return nothing
    end
    tmps = zero(s[k])
    tmpc = zero(c[k])
    for i = 1:k
        x = i * a[i]
        tmps += x * c[k-i]
        tmpc += x * s[k-i]
    end
    s[k] = tmps / k
    c[k] = tmpc / k
    return nothing
end

function TS.tanh!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        c2::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    if k == 0
        aux = tanh( constant_term(a) )
        c[0] = aux
        c2[0] = aux^2
        return nothing
    end
    tmp = zero(c[k])
    for i = 0:k-1
        tmp += (k-i) * a[k-i] * c2[i]
    end
    c[k] = a[k] - tmp/k
    # c2 = c^2
    TS.sqr!(c2, c, tmp, k)
    return nothing
end
