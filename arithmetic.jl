# arithmetic.jl

zero(a::TMAbsRem) = TMAbsRem(zero(a.pol), zero(a.arem), a.x0, a.iI)

# iszero(a::TMAbsRem) = iszero(a.pol) && iszero(zero(a.arem))

==(a::TMAbsRem, b::TMAbsRem) =
    a.pol == b.pol && a.arem == b.arem && a.x0 == b.x0 && a.iI == b.iI


# Addition
+(a::TMAbsRem) = a

function +(a::TMAbsRem, b::TMAbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI
    return TMAbsRem(a.pol+b.pol, a.arem+b.arem, a.x0, a.iI)
end

+(a::TMAbsRem, b::T) where {T} = TMAbsRem(a.pol+b, a.arem, a.x0, a.iI)

+(b::T, a::TMAbsRem) where {T} = TMAbsRem(b+a.pol, a.arem, a.x0, a.iI)


# Substraction
-(a::TMAbsRem) = TMAbsRem(-a.pol, -a.arem, a.x0, a.iI)

function -(a::TMAbsRem, b::TMAbsRem)
    @assert a.x0 == b.x0 && a.iI == b.iI
    return TMAbsRem(a.pol-b.pol, a.arem-b.arem, a.x0, a.iI)
end

-(a::TMAbsRem, b::T) where {T} = TMAbsRem(a.pol-b, a.arem, a.x0, a.iI)

-(b::T, a::TMAbsRem) where {T} = TMAbsRem(b-a.pol, -a.arem, a.x0, a.iI)


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
    Δ = Δnegl + Δb * a.arem + Δa * b.arem + a.arem * b.arem

    # Returned polynomial
    polret = Taylor1( copy(res.coeffs[1:order+1]) )

    return TMAbsRem(polret, Δ, a.x0, a.iI)
end

*(a::TMAbsRem, b::T) where {T} = TMAbsRem(a.pol*b, b*a.arem, a.x0, a.iI)

*(b::T, a::TMAbsRem) where {T} = a * b


# Division
function /(a::TMAbsRem, b::TMAbsRem)
    invb = rpa(x->inv(x), b)
    return a * invb
end

/(a::TMAbsRem, b::T) where {T} = a * inv(b)

/(b::T, a::TMAbsRem) where {T} = b * inv(a)


# Power
^(a::TMAbsRem, r) = rpa(x->x^r, a)
^(a::TMAbsRem, n::Integer) = rpa(x->x^n, a)
