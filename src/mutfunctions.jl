# Mutating functions
#

function TS.zero!(a::TaylorModelN)
    for k in eachindex(a)
        TS.zero!(a, k)
    end
    return nothing
end
function TS.zero!(a::TaylorModelN, k::Int)
    TS.zero!(a[k])
    return nothing
end

#
function TS.pow!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}},
        tmp::Taylor1{TaylorModelN{N,T,S}}, r::S, k::Int, l0::Int=0) where {N,T,S}
    if k == l0
        c[0] = ( a[l0] )^r
        return nothing
    end
    tmp[k-l0] = zero(c[k-l0])
    for i = 0:k-l0-1
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
    aux = zero(a[0])
    @inbounds for ord in eachindex(c[0])
        sqr!(c[0], a[0], aux, ord)
    end
    return nothing
end

function TS.sqrt!(c::Taylor1{TaylorModelN{N,T,S}},
        a::Taylor1{TaylorModelN{N,T,S}}, k::Int, k0::Int=0) where {N,T,S}
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
