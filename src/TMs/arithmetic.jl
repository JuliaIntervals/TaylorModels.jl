# arithmetic.jl

# Addition, substraction and other functions
for TM in tupleTMs
    @eval begin
        zero(a::$TM) = $TM(zero(a.pol), zero(a.rem), a.x0, a.iI)
        one(a::$TM) = $TM(one(a.pol), zero(a.rem), a.x0, a.iI)

        # iszero(a::$TM) = iszero(a.pol) && iszero(zero(a.rem))

        findfirst(a::$TM) = findfirst(a.pol)

        ==(a::$TM, b::$TM) =
            a.pol == b.pol && a.rem == b.rem && a.x0 == b.x0 && a.iI == b.iI


        # Addition
        +(a::$TM) = $TM(a.pol, a.rem, a.x0, a.iI)

        function +(a::$TM, b::$TM)
            @assert a.x0 == b.x0 && a.iI == b.iI
            return $TM(a.pol+b.pol, a.rem+b.rem, a.x0, a.iI)
        end

        +(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol+b, a.rem, a.x0, a.iI)

        +(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b+a.pol, a.rem, a.x0, a.iI)


        # Substraction
        -(a::$TM) = $TM(-a.pol, -a.rem, a.x0, a.iI)

        function -(a::$TM, b::$TM)
            @assert a.x0 == b.x0 && a.iI == b.iI
            return $TM(a.pol-b.pol, a.rem-b.rem, a.x0, a.iI)
        end

        -(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol-b, a.rem, a.x0, a.iI)

        -(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b-a.pol, -a.rem, a.x0, a.iI)


        # Basic division
        function basediv(a::$TM, b::$TM)
            invb = rpa(x->inv(x), b)
            return a * invb
        end


        # Multiplication by numbers
        *(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol*b, b*a.rem, a.x0, a.iI)

        *(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(a.pol*b, b*a.rem, a.x0, a.iI)


        # Division by numbers
        /(a::$TM, b::T) where {T<:NumberNotSeries} = a * inv(b)

        /(b::T, a::$TM) where {T<:NumberNotSeries} = b * inv(a)


        # Power
        ^(a::$TM, r) = rpa(x->x^r, a)

        ^(a::$TM, n::Integer) = rpa(x->x^n, a)
    end
end

# This requires current master of TaylorSeries
# # In TaylorSeries v0.7.3, a^n, a::Taylor1 and n::Integer,
# # in general use pow!, which yields [-∞,∞] if the interval
# # contains zero; the following uses power_by_squaring
# ^(a::Taylor1{Interval{T}}, n::Integer) where {T} = Base.power_by_squaring(a, n)
# ^(a::TaylorN{Interval{T}}, n::Integer) where {T} = Base.power_by_squaring(a, n)


# Multiplication
function *(a::TM1AbsRem, b::TM1AbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), 2*order)
    bext = Taylor1(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product of polynomials
    res[0:order] .= zero(eltype(res))
    aux = a.iI - a.x0
    Δnegl = res(aux)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    return TM1AbsRem(bext, Δ, a.x0, a.iI)
end

function *(a::TM1RelRem, b::TM1RelRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), 2*order)
    bext = Taylor1(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product (properly factorized)
    res = Taylor1(copy(res.coeffs[order+2:2*order+1]), order-1)
    aux = a.iI - a.x0
    Δnegl = res(aux)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    V = aux^(order+1)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V

    return TM1RelRem(bext, Δ, a.x0, a.iI)
end


# Division
function /(a::TM1AbsRem, b::TM1AbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI
    return basediv(a, b)
end

function /(a::TM1RelRem, b::TM1RelRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # DetermineRootOrderUpperBound seems equivalent (optimized?) to `findfirst`
    bk = findfirst(b)
    bk ≤ 0 && return basediv(a, b)

    # In 2.3.12, `a` and `b` are now extended in order by `bk`,
    # but how can we do this without knowing the explicit function
    # that created them?
    #
    # Below we reduce the original order by bk.
    #
    order = get_order(a)
    ared = reducetoorder(
        TM1RelRem(Taylor1(a.pol.coeffs[bk+1:order+1]), a.rem, a.x0, a.iI), order-bk)
    order = get_order(b)
    bred = reducetoorder(
        TM1RelRem(Taylor1(b.pol.coeffs[bk+1:order+1]), b.rem, b.x0, b.iI), order-bk)

    return basediv( ared, bred )
end


"""
    reducetoorder(a::TM1RelRem, m::Int)

From `a:TM1RelRem`, it returns a Taylor Model of relative
remainder of order `m`.

"""
function reducetoorder(a::TM1RelRem, m::Int)
    order = get_order(a)
    @assert order ≥ m ≥ 0

    bpol = Taylor1(copy(a.pol.coeffs))
    zz = zero(bpol[0])
    for ind in 0:m
        bpol[ind] = zz
    end
    bf = bpol(a.iI-a.x0)
    Δ = bf + a.rem * (a.iI-a.x0)^(order-m)
    return TM1RelRem( Taylor1(copy(a.pol.coeffs[1:m+1])), Δ, a.x0, a.iI )
end


# Same as above, for TMNAbsRem
zero(a::TMNAbsRem) = TMNAbsRem(zero(a.pol), zero(a.rem), a.x0, a.iI)
one(a::TMNAbsRem) = TMNAbsRem(one(a.pol), zero(a.rem), a.x0, a.iI)

# iszero(a::TMNAbsRem) = iszero(a.pol) && iszero(zero(a.rem))

# findfirst(a::TMNAbsRem) = findfirst(a.pol)

==(a::TMNAbsRem, b::TMNAbsRem) =
    a.pol == b.pol && a.rem == b.rem && a.x0 == b.x0 && a.iI == b.iI


# Addition and substraction
for op in (:+, :-)
    @eval begin
        $(op)(a::TMNAbsRem) = TMNAbsRem($(op)(a.pol), $(op)(a.rem), a.x0, a.iI)

        function $(op)(a::TMNAbsRem, b::TMNAbsRem)
            @assert a.x0 == b.x0 && a.iI == b.iI
            return TMNAbsRem($(op)(a.pol,b.pol), $(op)(a.rem,b.rem), a.x0, a.iI)
        end

        $(op)(a::TMNAbsRem, b::T) where {T<:NumberNotSeries} =
            TMNAbsRem($(op)(a.pol, b), a.rem, a.x0, a.iI)

        $(op)(b::T, a::TMNAbsRem) where {T<:NumberNotSeries} =
            TMNAbsRem($(op)(b, a.pol), $(op)(a.rem), a.x0, a.iI)
    end
end


# Multiplication
function *(a::TMNAbsRem, b::TMNAbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    @assert 2*order ≤ get_order()
    aext = TaylorN(copy(a.pol.coeffs), 2*order)
    bext = TaylorN(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Returned polynomial
    bext = TaylorN( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product of polynomials
    res[0:order] .= zero(eltype(res))
    aux = a.iI - a.x0
    Δnegl = res(aux)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    return TMNAbsRem(bext, Δ, a.x0, a.iI)
end


# Multiplication by numbers
function *(b::T, a::TMNAbsRem) where {T<:NumberNotSeries}
    pol = a.pol * b
    rem = b * a.rem
    return TMNAbsRem(pol, rem, a.x0, a.iI)
end

*(a::TMNAbsRem, b::T) where {T<:NumberNotSeries} = b * a


# Basic division
function basediv(a::TMNAbsRem, b::TMNAbsRem)
    invb = rpa(x->inv(x), b)
    return a * invb
end


# Division by numbers
/(a::TMNAbsRem, b::T) where {T<:Number} = a * inv(b)
/(b::T, a::TMNAbsRem) where {T<:NumberNotSeries} = b * inv(a)


# Power
^(a::TMNAbsRem, r) = rpa(x->x^r, a)
^(a::TMNAbsRem, n::Integer) = rpa(x->x^n, a)
