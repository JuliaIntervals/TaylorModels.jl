# arithmetic.jl

for TM in tupleTMs
    @eval begin
        zero(a::$TM) = $TM(zero(a.pol), zero(a.rem), a.x0, a.iI)

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

        +(a::$TM, b::T) where {T} = $TM(a.pol+b, a.rem, a.x0, a.iI)
        +(b::T, a::$TM) where {T} = $TM(b+a.pol, a.rem, a.x0, a.iI)


        # Substraction
        -(a::$TM) = $TM(-a.pol, -a.rem, a.x0, a.iI)

        function -(a::$TM, b::$TM)
            @assert a.x0 == b.x0 && a.iI == b.iI
            return $TM(a.pol-b.pol, a.rem-b.rem, a.x0, a.iI)
        end

        -(a::$TM, b::T) where {T} = $TM(a.pol-b, a.rem, a.x0, a.iI)
        -(b::T, a::$TM) where {T} = $TM(b-a.pol, -a.rem, a.x0, a.iI)


        # Base div
        function basediv(a::$TM, b::$TM)
            invb = rpa(x->inv(x), b)
            return a * invb
        end


        # Multiplication by numbers
        *(a::$TM, b::T) where {T} = $TM(a.pol*b, b*a.rem, a.x0, a.iI)
        *(b::T, a::$TM) where {T} = $TM(a.pol*b, b*a.rem, a.x0, a.iI)


        # Division by numbers
        /(a::$TM, b::T) where {T} = a * inv(b)
        /(b::T, a::$TM) where {T} = b * inv(a)


        # Power
        ^(a::$TM, r) = rpa(x->x^r, a)
        ^(a::$TM, n::Integer) = rpa(x->x^n, a)
    end
end


# Multiplication
function *(a::TMAbsRem, b::TMAbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), 2*order)
    bext = Taylor1(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Neglected polynomial resulting from the product
    polnegl = Taylor1(zero(eltype(res)), 2*order)
    polnegl.coeffs[order+1:2order] .= res.coeffs[order+1:2order]

    # Remainder of the product
    Δnegl = polnegl(a.iI-a.x0)
    Δa = a.pol(a.iI-a.x0)
    Δb = b.pol(a.iI-a.x0)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    # Returned polynomial
    polret = Taylor1( copy(res.coeffs[1:order+1]) )

    return TMAbsRem(polret, Δ, a.x0, a.iI)
end
function *(a::TMRelRem, b::TMRelRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # Polynomial product with extended order
    order = max(get_order(a), get_order(b))
    aext = Taylor1(copy(a.pol.coeffs), 2*order)
    bext = Taylor1(copy(b.pol.coeffs), 2*order)
    res = aext * bext

    # Neglected polynomial resulting from the product
    polnegl = Taylor1(zero(eltype(res)), 2*order)
    polnegl.coeffs[order+1:2order] .= res.coeffs[order+1:2order]

    # Remainder of the product
    Δnegl = polnegl(a.iI-a.x0)
    Δa = a.pol(a.iI-a.x0)
    Δb = b.pol(a.iI-a.x0)
    V = (a.iI-a.x0)^(order+1)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V

    # Returned polynomial
    polret = Taylor1( copy(res.coeffs[1:order+1]) )

    return TMRelRem(polret, Δ, a.x0, a.iI)
end


# Division
function /(a::TMAbsRem, b::TMAbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI
    return basediv(a,b)
end
function /(a::TMAbsRem, b::TMAbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI

    # DetermineRootOrderUpperBound seems equivalent (optimized?) to `findfirst`
    bk = findfirst(b)
    bk ≤ 0 && return basediv(a,b)


    # TMRelRem
    # TMReduceOrder
    # TMBaseDiv
    invb = rpa(x->inv(x), b)
    return a * invb
end
