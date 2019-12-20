# arithmetic.jl

# Addition, substraction and other functions
for TM in tupleTMs
    @eval begin
        tmdata(f::$TM) = (f.x0, f.dom)

        zero(a::$TM) = $TM(zero(a.pol), zero(a.rem), a.x0, a.dom)
        one(a::$TM) = $TM(one(a.pol), zero(a.rem), a.x0, a.dom)

        # iszero(a::$TM) = iszero(a.pol) && iszero(zero(a.rem))

        findfirst(a::$TM) = findfirst(a.pol)

        ==(a::$TM, b::$TM) =
            a.pol == b.pol && a.rem == b.rem && a.x0 == b.x0 && a.dom == b.dom


        # Addition
        +(a::$TM) = $TM(a.pol, a.rem, a.x0, a.dom)

        function +(a::$TM, b::$TM)
            a, b = fixorder(a, b)
            return $TM(a.pol+b.pol, a.rem+b.rem, a.x0, a.dom)
        end

        +(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol+b, a.rem, a.x0, a.dom)

        +(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b+a.pol, a.rem, a.x0, a.dom)


        # Substraction
        -(a::$TM) = $TM(-a.pol, -a.rem, a.x0, a.dom)

        function -(a::$TM, b::$TM)
            a, b = fixorder(a, b)
            return $TM(a.pol-b.pol, a.rem-b.rem, a.x0, a.dom)
        end

        -(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol-b, a.rem, a.x0, a.dom)

        -(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b-a.pol, -a.rem, a.x0, a.dom)


        # Basic division
        function basediv(a::$TM, b::$TM)
            invb = rpa(x->inv(x), b)
            return a * invb
        end


        # Multiplication by numbers
        *(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol*b, b*a.rem, a.x0, a.dom)

        *(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(a.pol*b, b*a.rem, a.x0, a.dom)


        # Division by numbers
        /(a::$TM, b::T) where {T<:NumberNotSeries} = a * inv(b)

        /(b::T, a::$TM) where {T<:NumberNotSeries} = b * inv(a)


        # Power
        ^(a::$TM, r) = rpa(x->x^r, a)

        ^(a::$TM, n::Integer) = rpa(x->x^n, a)
    end
end


# Multiplication
function *(a::TaylorModel1, b::TaylorModel1)
    @assert tmdata(a) == tmdata(b)
    order, ordermax = minmax(get_order(a), get_order(b))

    # Polynomial product with extended order
    aext = Taylor1(copy(a.pol.coeffs), order+ordermax)
    bext = Taylor1(copy(b.pol.coeffs), order+ordermax)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product of polynomials
    aux = a.dom - a.x0
    Δnegl = bound_truncation(TaylorModel1, res, aux, order)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    return TaylorModel1(bext, Δ, a.x0, a.dom)
end

function *(a::TaylorModel1{TaylorModelN{N,T,S},S},
        b::TaylorModel1{TaylorModelN{N,T,S},S}) where{N,T,S}
    @assert tmdata(a) == tmdata(b)
    order, ordermax = minmax(get_order(a), get_order(b))

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), order+ordermax)
    bext = Taylor1(copy(b.pol.coeffs), order+ordermax)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product of polynomials
    aux = a.dom - a.x0
    Δnegl = bound_truncation(TaylorModel1, res, aux, order)
    # Δnegl = res(aux)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    # Evaluate at the TMN centered domain
    auxN = a[0].dom - a[0].x0
    ΔN = Δ(auxN)

    return TaylorModel1(bext, ΔN, a.x0, a.dom)
end

function *(a::RTaylorModel1, b::RTaylorModel1)
    @assert tmdata(a) == tmdata(b)
    order, ordermax = minmax(get_order(a), get_order(b))

    # Polynomial product with extended order
    aext = Taylor1(copy(a.pol.coeffs), order+ordermax)
    bext = Taylor1(copy(b.pol.coeffs), order+ordermax)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product (properly factorized)
    aux = a.dom - a.x0
    Δnegl = bound_truncation(RTaylorModel1, res, aux, order)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    V = aux^(order+1)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V

    return RTaylorModel1(bext, Δ, a.x0, a.dom)
end


# Division
function /(a::TaylorModel1, b::TaylorModel1)
    @assert tmdata(a) == tmdata(b)
    return basediv(a, b)
end

function /(a::RTaylorModel1, b::RTaylorModel1)
    @assert tmdata(a) == tmdata(b)

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
    ared = truncate_taylormodel(
        RTaylorModel1(Taylor1(a.pol.coeffs[bk+1:order+1]), a.rem, a.x0, a.dom), order-bk)
    order = get_order(b)
    bred = truncate_taylormodel(
        RTaylorModel1(Taylor1(b.pol.coeffs[bk+1:order+1]), b.rem, b.x0, b.dom), order-bk)

    return basediv( ared, bred )
end


"""
    truncate_taylormodel(a::RTaylorModel1, m::Integer)

Truncates `a::RTaylorModel1` to order `m`.

"""
function truncate_taylormodel(a::RTaylorModel1, m::Integer)
    order = get_order(a)
    m ≥ order && return a

    apol = Taylor1(copy(a.pol.coeffs[1:m+1]))
    bpol = Taylor1(copy(a.pol.coeffs))
    aux = a.dom - a.x0
    Δnegl = bound_truncation(RTaylorModel1, bpol, aux, m)
    Δ = Δnegl + a.rem * (aux)^(order-m)
    return RTaylorModel1( apol, Δ, a.x0, a.dom )
end


# Same as above, for TaylorModelN
tmdata(f::TaylorModelN) = (f.x0, f.dom)
zero(a::TaylorModelN) = TaylorModelN(zero(a.pol), zero(a.rem), a.x0, a.dom)
one(a::TaylorModelN) = TaylorModelN(one(a.pol), zero(a.rem), a.x0, a.dom)

# iszero(a::TaylorModelN) = iszero(a.pol) && iszero(zero(a.rem))

findfirst(a::TaylorModelN) = findfirst(a.pol)

==(a::TaylorModelN, b::TaylorModelN) =
    a.pol == b.pol && a.rem == b.rem && a.x0 == b.x0 && a.dom == b.dom


# Addition and substraction
for op in (:+, :-)
    @eval begin
        $(op)(a::TaylorModelN) = TaylorModelN($(op)(a.pol), $(op)(a.rem), a.x0, a.dom)

        function $(op)(a::TaylorModelN, b::TaylorModelN)
            a, b = fixorder(a, b)
            return TaylorModelN($(op)(a.pol,b.pol), $(op)(a.rem,b.rem), a.x0, a.dom)
        end

        $(op)(a::TaylorModelN, b::T) where {T<:NumberNotSeries} =
            TaylorModelN($(op)(a.pol, b), a.rem, a.x0, a.dom)

        $(op)(b::T, a::TaylorModelN) where {T<:NumberNotSeries} =
            TaylorModelN($(op)(b, a.pol), $(op)(a.rem), a.x0, a.dom)
    end
end


# Multiplication
function *(a::TaylorModelN, b::TaylorModelN)
    @assert tmdata(a) == tmdata(b)
    order, ordermax = minmax(get_order(a), get_order(b))

    # Polynomial product with extended order
    @assert order+ordermax ≤ get_order()
    aext = TaylorN(copy(a.pol.coeffs), order+ordermax)
    bext = TaylorN(copy(b.pol.coeffs), order+ordermax)
    res = aext * bext

    # Returned polynomial
    bext = TaylorN( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product of polynomials
    aux = a.dom - a.x0
    Δnegl = bound_truncation(TaylorModelN, res, aux, order)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    return TaylorModelN(bext, Δ, a.x0, a.dom)
end


# Multiplication by numbers
function *(b::T, a::TaylorModelN) where {T<:NumberNotSeries}
    pol = a.pol * b
    rem = b * a.rem
    return TaylorModelN(pol, rem, a.x0, a.dom)
end

*(a::TaylorModelN, b::T) where {T<:NumberNotSeries} = b * a


# Basic division
function basediv(a::TaylorModelN, b::TaylorModelN)
    invb = rpa(x->inv(x), b)
    return a * invb
end
# /(a::TaylorModelN, b::TaylorModelN) = basediv(a, b)


# Division by numbers
/(a::TaylorModelN, b::T) where {T<:Number} = a * inv(b)
/(b::T, a::TaylorModelN) where {T<:NumberNotSeries} = b * inv(a)


# Power
^(a::TaylorModelN, r) = rpa(x->x^r, a)
^(a::TaylorModelN, n::Integer) = rpa(x->x^n, a)
