# arithmetic.jl

# Addition, substraction and other functions
for TM in tupleTMs
    @eval begin
        tmdata(f::$TM) = (f.x0, f.I)

        zero(a::$TM) = $TM(zero(a.pol), zero(a.rem), a.x0, a.I)
        one(a::$TM) = $TM(one(a.pol), zero(a.rem), a.x0, a.I)

        # iszero(a::$TM) = iszero(a.pol) && iszero(zero(a.rem))

        findfirst(a::$TM) = findfirst(a.pol)

        ==(a::$TM, b::$TM) =
            a.pol == b.pol && a.rem == b.rem && a.x0 == b.x0 && a.I == b.I


        # Addition
        +(a::$TM) = $TM(a.pol, a.rem, a.x0, a.I)

        function +(a::$TM, b::$TM)
            @assert tmdata(a) == tmdata(b)
            return $TM(a.pol+b.pol, a.rem+b.rem, a.x0, a.I)
        end

        +(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol+b, a.rem, a.x0, a.I)

        +(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b+a.pol, a.rem, a.x0, a.I)


        # Substraction
        -(a::$TM) = $TM(-a.pol, -a.rem, a.x0, a.I)

        function -(a::$TM, b::$TM)
            @assert tmdata(a) == tmdata(b)
            return $TM(a.pol-b.pol, a.rem-b.rem, a.x0, a.I)
        end

        -(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol-b, a.rem, a.x0, a.I)

        -(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b-a.pol, -a.rem, a.x0, a.I)


        # Basic division
        function basediv(a::$TM, b::$TM)
            invb = rpa(x->inv(x), b)
            return a * invb
        end


        # Multiplication by numbers
        *(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol*b, b*a.rem, a.x0, a.I)

        *(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(a.pol*b, b*a.rem, a.x0, a.I)


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

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), 2*order)
    bext = Taylor1(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product of polynomials
    res[0:order] .= zero(eltype(res))
    aux = a.I - a.x0
    Δnegl = res(aux)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    return TaylorModel1(bext, Δ, a.x0, a.I)
end

function *(a::RTaylorModel1, b::RTaylorModel1)
    @assert tmdata(a) == tmdata(b)

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), 2*order)
    bext = Taylor1(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Returned polynomial
    bext = Taylor1( copy(res.coeffs[1:order+1]) )

    # Bound for the neglected part of the product (properly factorized)
    res = Taylor1(copy(res.coeffs[order+2:2*order+1]), order-1)
    aux = a.I - a.x0
    Δnegl = res(aux)

    # Remainder of the product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    V = aux^(order+1)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V

    return RTaylorModel1(bext, Δ, a.x0, a.I)

end


# Division
function /(a::TaylorModel1, b::TaylorModel1)
    @assert a.x0 == b.x0 && a.I == b.I
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
    ared = reducetoorder(
        RTaylorModel1(Taylor1(a.pol.coeffs[bk+1:order+1]), a.rem, a.x0, a.I), order-bk)
    order = get_order(b)
    bred = reducetoorder(
        RTaylorModel1(Taylor1(b.pol.coeffs[bk+1:order+1]), b.rem, b.x0, b.I), order-bk)

    return basediv( ared, bred )
end


"""
    reducetoorder(a::RTaylorModel1, m::Integer)

From `a::RTaylorModel1`, it returns a RTaylorModel1 of order `m`.

"""
function reducetoorder(a::RTaylorModel1, m::Integer)
    order = get_order(a)
    @assert order ≥ m ≥ 0

    bpol = Taylor1(copy(a.pol.coeffs))
    zz = zero(bpol[0])
    for ind in 0:m
        bpol[ind] = zz
    end
    bf = bpol(a.I-a.x0)
    Δ = bf + a.rem * (a.I-a.x0)^(order-m)
    return RTaylorModel1( Taylor1(copy(a.pol.coeffs[1:m+1])), Δ, a.x0, a.I )
end


# Same as above, for TaylorModelN
zero(a::TaylorModelN) = TaylorModelN(zero(a.pol), zero(a.rem), a.x0, a.I)
one(a::TaylorModelN) = TaylorModelN(one(a.pol), zero(a.rem), a.x0, a.I)

# iszero(a::TaylorModelN) = iszero(a.pol) && iszero(zero(a.rem))

findfirst(a::TaylorModelN) = findfirst(a.pol)

==(a::TaylorModelN, b::TaylorModelN) =
    a.pol == b.pol && a.rem == b.rem && a.x0 == b.x0 && a.I == b.I


# Addition and substraction
for op in (:+, :-)
    @eval begin
        $(op)(a::TaylorModelN) = TaylorModelN($(op)(a.pol), $(op)(a.rem), a.x0, a.I)

        function $(op)(a::TaylorModelN, b::TaylorModelN)
            @assert a.x0 == b.x0 && a.I == b.I
            return TaylorModelN($(op)(a.pol,b.pol), $(op)(a.rem,b.rem), a.x0, a.I)
        end

        $(op)(a::TaylorModelN, b::T) where {T<:NumberNotSeries} =
            TaylorModelN($(op)(a.pol, b), a.rem, a.x0, a.I)

        $(op)(b::T, a::TaylorModelN) where {T<:NumberNotSeries} =
            TaylorModelN($(op)(b, a.pol), $(op)(a.rem), a.x0, a.I)
    end
end


# Multiplication
function *(a::TaylorModelN, b::TaylorModelN)
    @assert a.x0 == b.x0 && a.I == b.I
    aux = a.I - a.x0

    # Some definitions related to the order of some polynomials
    orderTS = get_order()  # maximum order, fixed by `TaylorSeries._params_TaylorN_`
    order_a = get_order(a)
    order_b = get_order(b)
    order_prod = 2 * max(order_a, order_b)

    # The returned polynomial has no terms beyond `orderTS`
    if order_prod ≤ orderTS
        apol = TaylorN(a.pol.coeffs, order_prod)
        bpol = TaylorN(b.pol.coeffs, order_prod)
        res = apol * bpol
        Δa = a.pol(aux)
        Δb = b.pol(aux)
        Δ = Δb * a.rem + Δa * b.rem + a.rem * b.rem
        return TaylorModelN(res, Δ, a.x0, a.I)
    end

    # The resulting product has a part whose order is larger than `orderTS`;
    # a bound for this polynomial has to be included
    apol = TaylorN(a.pol.coeffs, orderTS)
    bpol = TaylorN(b.pol.coeffs, orderTS)
    res = apol * bpol    # Trunctated polynomial to order `orderTS`
    Δa = apol(aux)
    Δb = bpol(aux)
    Δ = Δb * a.rem + Δa * b.rem + a.rem * b.rem
    N = get_numvars()

    # We evaluate *explicitely* each term of the product of coefficients at `aux`
    Δnegl = zero(Δ)
    for order = orderTS+1:order_prod
        orderini = order - orderTS
        for inda = orderini:order_a
            indb = order - inda
            # Δnegl += evaluate(apol[inda], aux) * evaluate(bpol[indb], aux)

            # Coefficient tables of corresponding `HomogeneousPolynomial`s
            @inbounds cta = TaylorSeries.coeff_table[inda+1]
            @inbounds ctb = TaylorSeries.coeff_table[indb+1]

            for i in eachindex(cta)
                @inbounds a_coeff = apol.coeffs[inda+1][i]
                iszero(a_coeff) && continue
                for j in eachindex(ctb)
                    @inbounds b_coeff = bpol.coeffs[indb+1][j]
                    iszero(b_coeff) && continue

                    # Evaluation of monomials
                    tmp = one(Δ)
                    for n in 1:N
                        @inbounds tmp *= aux[n]^(cta[i][n]+ctb[j][n])
                    end
                    Δnegl += a_coeff * b_coeff * tmp
                end
            end
        end
    end

    return TaylorModelN(res, Δ+Δnegl, a.x0, a.I)
end


# Multiplication by numbers
function *(b::T, a::TaylorModelN) where {T<:NumberNotSeries}
    pol = a.pol * b
    rem = b * a.rem
    return TaylorModelN(pol, rem, a.x0, a.I)
end

*(a::TaylorModelN, b::T) where {T<:NumberNotSeries} = b * a


# Basic division
function basediv(a::TaylorModelN, b::TaylorModelN)
    invb = rpa(x->inv(x), b)
    return a * invb
end


# Division by numbers
/(a::TaylorModelN, b::T) where {T<:Number} = a * inv(b)
/(b::T, a::TaylorModelN) where {T<:NumberNotSeries} = b * inv(a)


# Power
^(a::TaylorModelN, r) = rpa(x->x^r, a)
^(a::TaylorModelN, n::Integer) = rpa(x->x^n, a)
