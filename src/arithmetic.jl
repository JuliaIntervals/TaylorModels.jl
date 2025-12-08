# arithmetic.jl

# Addition, substraction and other functions
for TM in tupleTMs
    @eval begin
        zero(a::$TM{T,S}) where {T,S} =
            $TM(zero(a.pol), zero(remainder(a)), expansion_point(a), domain(a))
        one(a::$TM{T,S}) where {T,S} =
            $TM(one(a.pol), zero(remainder(a)), expansion_point(a), domain(a))

        # iszero(a::$TM) = iszero(a.pol) && iszero(zero(remainder(a)))

        findfirst(a::$TM{T,S}) where {T,S} = findfirst(a.pol)

        ==(a::$TM, b::$TM) = all(isequal_interval.(tmdata(a), tmdata(b))) &&
            isequal_interval(remainder(a), remainder(b)) && a.pol == b.pol


        # Addition
        +(a::$TM{T,S}) where {T,S} = $TM(a.pol, remainder(a), expansion_point(a), domain(a))

        function +(a::$TM, b::$TM)
            a, b = fixorder(a, b)
            return $TM(a.pol+b.pol, remainder(a)+remainder(b), expansion_point(a), domain(a))
        end

        +(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol+b, remainder(a),
            expansion_point(a), domain(a))

        +(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b+a.pol, remainder(a),
            expansion_point(a), domain(a))


        # Substraction
        -(a::$TM{T,S}) where {T,S} = $TM(-a.pol, -remainder(a), expansion_point(a), domain(a))

        function -(a::$TM, b::$TM)
            a, b = fixorder(a, b)
            return $TM(a.pol-b.pol, remainder(a)-remainder(b), expansion_point(a), domain(a))
        end

        -(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol-b, remainder(a),
            expansion_point(a), domain(a))

        -(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b-a.pol, -remainder(a),
            expansion_point(a), domain(a))


        # Basic division
        function basediv(a::$TM, b::$TM)
            # invb = rpa(x->inv(x), b)
            invb = inv(b)
            return a * invb
        end


        # Multiplication by numbers
        *(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol*b, b*remainder(a),
            expansion_point(a), domain(a))

        *(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(a.pol*b, b*remainder(a),
            expansion_point(a), domain(a))

        # Multiplication
        function *(a::$TM{T,S}, b::$TM{T,S}) where {T,S}
            @assert all(isequal_interval.(tmdata(a), tmdata(b)))
            # Returned polynomial
            res = a.pol * b.pol
            # Remainder
            Δ = remainder_product(a, b, get_order(res))
            return $TM(res, Δ, expansion_point(a), domain(a))
        end

        # Division by numbers
        /(a::$TM, b::T) where {T<:NumberNotSeries} = a * inv(b)

        /(b::T, a::$TM) where {T<:NumberNotSeries} = b * inv(a)

        # Power
        # Square
        function TS.square(a::$TM{T,S}) where {T,S}
            # Returned polynomial
            res = TS.square(a.pol)
            # Remainder
            Δ = remainder_square(a, get_order(res))
            return $TM(res, Δ, expansion_point(a), domain(a))
        end

        function ^(a::$TM{T,S}, r::NumberNotSeries) where {T,S}
            isinteger(r) && return a^Int(r)
            return rpa(x->x^r, a)
        end

        function ^(a::$TM{T,S}, n::Integer) where {T,S}
            iszero(n) && return one(a)
            n == 1 && return a
            n == 2 && return TS.square(a)
            n > 0 && return Base.power_by_squaring(a, n)
            return rpa(x->x^n, a)
        end

        function Base.power_by_squaring(x::$TM, p::Integer)
            @assert p ≥ 0
            (p == 0) && return one(x)
            (p == 1) && return $TM(
                Taylor1(x.pol.coeffs[:], x.pol.order), x.rem, x.x0, x.dom)
            (p == 2) && return TS.square(x)
            (p == 3) && return x*TS.square(x)
            t = trailing_zeros(p) + 1
            p >>= t
            while (t -= 1) > 0
                x = TS.square(x)
            end
            y = x
            while p > 0
                t = trailing_zeros(p) + 1
                p >>= t
                while (t -= 1) ≥ 0
                    x = TS.square(x)
                end
                y *= x
            end
            return y
        end

    end
end


# Remainder of the product: checks the three (algebraically equivalent) forms
# to evaluate the remainder, and chooses the best (cf. Bunger 2020)
# TaylorModel1
function remainder_product(a::TaylorModel1{T,S}, b::TaylorModel1{T,S},
        order::Int) where {T<:NumberNotSeries,S}
    aux = centered_dom(a)
    if iseven(order)
        # Remaining terms of the product as a Taylor1 pol
        rnegl = _neglected_polynomial(a, b, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux)
    else
        # Remaining terms of the product as a factored Taylor1 pol
        rnegl = _neglected_polynomial_factored(a, b, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux) * aux^interval(order+1)
    end
    # Remainder of product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl + Δb * a_rem + Δa * b_rem + a_rem * b_rem
    Δ1 = Δnegl + Δb * a_rem + (Δa + a_rem) * b_rem
    Δ2 = Δnegl + Δa * b_rem + (Δb + b_rem) * a_rem
    return _intersect_reminders((Δ1, Δ2, Δ))
end
function remainder_product(a::TaylorModel1{TaylorN{T},S},
        b::TaylorModel1{TaylorN{T},S}, order::Int) where {T,S}
    auxT = centered_dom(a)
    if iseven(order)
        # Remaining terms of the product as a Taylor1 pol
        rnegl = _neglected_polynomial(a, b, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(auxT)
    else
        # Remaining terms of the product as a factored Taylor1 pol
        rnegl = _neglected_polynomial_factored(a, b, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(auxT) * auxT^interval(order+1)
    end
    # An N-dimensional symmetrical IntervalBox is assumed
    # to bound the TaylorN part
    auxQ = symmetric_box(S)
    #
    Δa = a.pol(auxT)(auxQ)
    Δb = b.pol(auxT)(auxQ)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl(auxQ) + Δb * a_rem + Δa * b_rem + a_rem * b_rem
    Δ1 = Δnegl(auxQ) + Δb * a_rem + (Δa + a_rem) * b_rem
    Δ2 = Δnegl(auxQ) + Δa * b_rem + (Δb + b_rem) * a_rem
    return _intersect_reminders((Δ1, Δ2, Δ))
end
function remainder_product(a::TaylorModel1{TaylorModelN{N,T,S},S},
        b::TaylorModel1{TaylorModelN{N,T,S},S}, order::Int) where {N,T,S}
    aux = centered_dom(a)
    if iseven(order)
        # Remaining terms of the product as a Taylor1 pol
        rnegl = _neglected_polynomial(a, b, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux)
    else
        # Remaining terms of the product as a factored Taylor1 pol
        rnegl = _neglected_polynomial_factored(a, b, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux) * aux^interval(order+1)
    end
    #
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl + Δb * a_rem + Δa * b_rem + a_rem * b_rem
    # Evaluate at the TMN centered domain
    auxN = centered_dom(a[0])
    ΔN = Δ(auxN)
    return ΔN
end
# RTaylorModel1
function remainder_product(a::RTaylorModel1{T,S}, b::RTaylorModel1{T,S},
        order::Int) where {T<:NumberNotSeries,S}
    # Remaining terms of the product as a Taylor1 pol
    rnegl = _neglected_polynomial_factored(a, b, order)
    # Remainder of the neglected terms
    aux = centered_dom(a)
    Δnegl = rnegl(aux)
    # Remainder of product
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    V = aux^interval(order+1)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V
    Δ1 = Δnegl + Δb * a_rem + (Δa + a_rem * V) * b_rem
    Δ2 = Δnegl + Δa * b_rem + (Δb + b_rem * V) * a_rem
    return _intersect_reminders((Δ1, Δ2, Δ))
end
function remainder_product(a::RTaylorModel1{TaylorN{T},S},
        b::RTaylorModel1{TaylorN{T},S}, order::Int) where {T,S}
    # Remaining terms of the product as a Taylor1 pol
    rnegl = _neglected_polynomial_factored(a, b, order)
    # Remainder of the neglected terms
    aux = centered_dom(a)
    Δnegl = rnegl(aux)
    # An N-dimensional symmetrical IntervalBox is assumed
    # to bound the TaylorN part
    auxQ = symmetric_box(S)
    Δa = a.pol(aux)(auxQ)
    Δb = b.pol(aux)(auxQ)
    V = aux^interval(order+1)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl(auxQ) + Δb * a_rem + Δa * b_rem + a_rem * b_rem * V
    Δ1 = Δnegl(auxQ) + Δb * a_rem + (Δa + a_rem * V) * b_rem
    Δ2 = Δnegl(auxQ) + Δa * b_rem + (Δb + b_rem * V) * a_rem
    return _intersect_reminders((Δ1, Δ2, Δ))
end

# TaylorModel1
function remainder_square(a::TaylorModel1{T,S}, order::Int) where
        {T<:NumberNotSeries,S}
    aux = centered_dom(a)
    if iseven(order)
        # Remaining terms of the product as a Taylor1 pol
        rnegl = _neglected_polynomial(a, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux)
    else
        # Remaining terms of the product as a factored Taylor1 pol
        rnegl = _neglected_polynomial_factored(a, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux) * aux^interval(order+1)
    end
    Δa = a.pol(aux)
    a_rem = remainder(a)
    Δ = Δnegl + interval(2) * Δa * a_rem + a_rem^interval(2)
    Δ1 = Δnegl + (interval(2) * Δa + a_rem) * a_rem
    return _intersect_reminders((Δ1, Δ))
end
function remainder_square(a::TaylorModel1{TaylorN{T}, S},
        order::Int) where {T<:NumberNotSeries,S}
    auxT = centered_dom(a)
    if iseven(order)
        # Remaining terms of the product as a Taylor1 pol
        rnegl = _neglected_polynomial(a, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(auxT)
    else
        # Remaining terms of the product as a factored Taylor1 pol
        rnegl = _neglected_polynomial_factored(a, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(auxT) * auxT^interval(order+1)
    end
    #
    # An N-dimensional symmetrical IntervalBox is assumed
    # to bound the TaylorN part
    auxQ = symmetric_box(S)
    Δa = a.pol(auxT)(auxQ)
    a_rem = remainder(a)
    Δ = Δnegl(auxQ) + interval(2) * Δa * a_rem + a_rem^interval(2)
    Δ1 = Δnegl(auxQ) + (interval(2) * Δa + a_rem) * a_rem
    return _intersect_reminders((Δ1, Δ))
end
function remainder_square(a::TaylorModel1{TaylorModelN{N,T,S},S},
        order::Int) where {N,T<:NumberNotSeries,S}
    aux = centered_dom(a)
    if iseven(order)
        # Remaining terms of the product as a Taylor1 pol
        rnegl = _neglected_polynomial(a, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux)
    else
        # Remaining terms of the product as a factored Taylor1 pol
        rnegl = _neglected_polynomial_factored(a, order)
        # Remainder of the neglected terms
        Δnegl = rnegl(aux) * aux^interval(order+1)
    end
    #
    Δa = a.pol(aux)
    a_rem = remainder(a)
    Δ = Δnegl + interval(2) * Δa * a_rem + a_rem^interval(2)
    # Evaluate at the TMN centered domain
    auxN = centered_dom(a[0])
    return Δ(auxN)
end
# RTaylorModel1
function remainder_square(a::RTaylorModel1{T,S}, order::Int) where
        {T<:NumberNotSeries,S}
    # Remaining terms of the product as a Taylor1 pol
    rnegl = _neglected_polynomial_factored(a, order)
    # Remainder of the neglected terms
    aux = centered_dom(a)
    Δnegl = rnegl(aux)
    #
    Δa = a.pol(aux)
    V = aux^interval(order+1)
    a_rem = remainder(a)
    Δ = Δnegl + interval(2) * Δa * a.rem + a.rem^interval(2) * V
    Δ1 = Δnegl + (interval(2) * Δa + a_rem * V) * a_rem
    return _intersect_reminders((Δ1, Δ))
end
function remainder_square(a::RTaylorModel1{TaylorN{T},S}, order::Int) where
        {T<:NumberNotSeries,S}
    # Remaining terms of the product as a Taylor1 pol
    rnegl = _neglected_polynomial_factored(a, order)
    # Remainder of the neglected terms
    aux = centered_dom(a)
    Δnegl = rnegl(aux)
    # An N-dimensional symmetrical IntervalBox is assumed
    # to bound the TaylorN part
    auxQ = symmetric_box(S)
    Δa = a.pol(aux)(auxQ)
    V = aux^interval(order+1)
    a_rem = remainder(a)
    Δ = Δnegl(auxQ) + interval(2) * Δa * a_rem + a_rem^interval(2) * V
    Δ1 = Δnegl(auxQ) + (interval(2) * Δa + a_rem * V) * a_rem
    return _intersect_reminders((Δ1, Δ))
end

# Neglected polynomial for the product and square
function _neglected_polynomial(a::TaylorModel1{T,S}, b::TaylorModel1{T,S},
        order::Int) where {T,S}
    a_order = a.pol.order
    b_order = b.pol.order
    rnegl_order = a_order + b_order
    # Remaining terms of the product as reduced Taylor1
    z = zero(a.pol[0])
    rnegl = Taylor1(z, rnegl_order)
    for k in order+1:rnegl_order
        tmp = z
        @inbounds for i = 0:k
            (i > a_order || k-i > b_order) && continue
            tmp += a.pol[i] * b.pol[k-i]
        end
        rnegl[k] = tmp
    end
    return rnegl
end
function _neglected_polynomial(a::TaylorModel1{T,S}, order::Int) where {T,S}
    z = zero(a.pol[0])
    rnegl = Taylor1(z, 2*order)
    for k in order+1:2*order
        tmp = z
        @inbounds for i = 0:k
            (i > order || k-i > order) && continue
            tmp += a.pol[i] * a.pol[k-i]
        end
        rnegl[k] = tmp
    end
    return rnegl
end
for TM in tupleTMs
    @eval function _neglected_polynomial_factored(a::$TM{T,S}, b::$TM{T,S},
            order::Int) where {T,S}
        a_order = a.pol.order
        b_order = b.pol.order
        rnegl_order = a_order + b_order
        # Remaining terms of the product as reduced Taylor1 (factored polynomial)
        z = zero(a.pol[0])
        rnegl = Taylor1(z, rnegl_order-order)
        for k in order+1:rnegl_order
            l = k-order-1
            tmp = z
            @inbounds for i = 0:k
                (i > a_order || k-i > b_order) && continue
                tmp += a.pol[i] * b.pol[k-i]
            end
            rnegl[l] = tmp
        end
        return rnegl
    end
    @eval function _neglected_polynomial_factored(a::$TM{T,S},
            order::Int) where {T,S}
        # Remaining terms of the product as reduced Taylor1 (factored polynomial)
        z = zero(a.pol[0])
        rnegl = Taylor1(z, order)
        for k in order+1:2*order
            tmp = z
            @inbounds for i = 0:k
                (i > order || k-i > order) && continue
                tmp += a.pol[i] * a.pol[k-i]
            end
            rnegl[k-order-1] = tmp
        end
        return rnegl
    end
end

function _intersect_reminders(Δs::NTuple{N,Interval{T}}) where{N,T}
    rems = (Δs..., intersect_interval(Δs...))
    _, i = findmin(diam, rems)
    return rems[i]
end


# Division
function /(a::TaylorModel1, b::TaylorModel1)
    @assert all(isequal_interval.(tmdata(a), tmdata(b)))
    return basediv(a, b)
end

function /(a::RTaylorModel1, b::RTaylorModel1)
    @assert all(isequal_interval.(tmdata(a), tmdata(b)))

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
        RTaylorModel1(Taylor1(a.pol.coeffs[bk+1:order+1]), remainder(a),
        expansion_point(a), domain(a)), order-bk)
    order = get_order(b)
    bred = truncate_taylormodel(
        RTaylorModel1(Taylor1(b.pol.coeffs[bk+1:order+1]), remainder(b),
        expansion_point(b), domain(b)), order-bk)

    return basediv( ared, bred )
end


"""
    truncate_taylormodel(a::RTaylorModel1, m::Integer)

Truncates `a::RTaylorModel1` to order `m`.

"""
function truncate_taylormodel(a::RTaylorModel1, m::Integer)
    order = get_order(a)
    m ≥ order && return a

    apol_coeffs = polynomial(a).coeffs
    apol = Taylor1(copy(apol_coeffs[1:m+1]))
    bpol = Taylor1(copy(apol_coeffs))
    aux = centered_dom(a)
    Δnegl = bound_truncation(RTaylorModel1, bpol, aux, m)
    Δ = Δnegl + remainder(a) * aux^interval(order-m)
    return RTaylorModel1( apol, Δ, expansion_point(a), domain(a) )
end


# Same as above, for TaylorModelN
zero(a::TaylorModelN) = TaylorModelN(zero(a.pol), zero(remainder(a)),
    expansion_point(a), domain(a))
one(a::TaylorModelN) = TaylorModelN(one(a.pol), zero(remainder(a)),
    expansion_point(a), domain(a))
zero(a::Taylor1{TaylorModelN{N,T,S}}) where {N,T,S} = Taylor1(zero(a[0]), a.order)

# iszero(a::TaylorModelN) = iszero(a.pol) && iszero(zero(remainder(a)))

findfirst(a::TaylorModelN) = findfirst(a.pol)

==(a::TaylorModelN, b::TaylorModelN) =
    a.pol == b.pol && isequal_interval(remainder(a), remainder(b)) &&
        all(isequal_interval.(tmdata(a), tmdata(b)))
        # expansion_point(a) == expansion_point(b) && domain(a) == domain(b)


# Addition and substraction
for op in (:+, :-)
    @eval begin
        $(op)(a::TaylorModelN) = TaylorModelN($(op)(a.pol), $(op)(remainder(a)),
            expansion_point(a), domain(a))

        function $(op)(a::TaylorModelN, b::TaylorModelN)
            a, b = fixorder(a, b)
            return TaylorModelN($(op)(a.pol, b.pol), $(op)(remainder(a), remainder(b)),
                expansion_point(a), domain(a))
        end

        $(op)(a::TaylorModelN, b::T) where {T<:NumberNotSeries} =
            TaylorModelN($(op)(a.pol, b), remainder(a), expansion_point(a), domain(a))

        $(op)(b::T, a::TaylorModelN) where {T<:NumberNotSeries} =
            TaylorModelN($(op)(b, a.pol), $(op)(remainder(a)), expansion_point(a), domain(a))
    end
end


# Multiplication
function *(a::TaylorModelN, b::TaylorModelN)
    @assert all(isequal_interval.(tmdata(a), tmdata(b)))
    @assert get_order(a)+get_order(b) ≤ get_order()
    # Returned polynomial
    res = a.pol * b.pol
    # Remainder
    Δ = remainder_product(a, b, get_order(res))
    return TaylorModelN(res, Δ, expansion_point(a), domain(a))
end

function remainder_product(a::TaylorModelN{N,T,S}, b::TaylorModelN{N,T,S},
        order::Int) where {N,T,S}
    # Remaining terms of the product as a Taylor1 pol
    suma = _neglected_polynomial(a, b, order)
    # Remainder of the neglected terms
    Δnegl = sum( suma ) # = sum( sort!(suma, by=abs2) )
    # Remainder of product
    aux = centered_dom(a)
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl + Δb * a_rem + Δa * b_rem + a_rem * b_rem
    Δ1 = Δnegl + Δb * a_rem + (Δa + a_rem) * b_rem
    Δ2 = Δnegl + Δa * b_rem + (Δb + b_rem) * a_rem
    return _intersect_reminders((Δ1, Δ2, Δ))
end
function _neglected_polynomial(a::TaylorModelN{N,T,S}, b::TaylorModelN{N,T,S},
        order::Int) where {N,T,S}
    # Remaining terms of the product
    a_order = a.pol.order
    b_order = b.pol.order
    rnegl_order = a_order + b_order
    aux = centered_dom(a)
    R = TS.numtype(a.pol[0]*b.pol[0])
    D = TS.numtype(domain(a))
    z = zero(R)
    suma = Vector{promote_type(R,D)}(undef, rnegl_order-order)
    for k in order+1:rnegl_order
        tmp = HomogeneousPolynomial(z, k)
        @inbounds for i = 0:k
            (i > a_order || k-i > b_order) && continue
            TS.mul!(tmp, a.pol[i], b.pol[k-i])
        end
        suma[k-order] = tmp(aux)
    end
    return suma
end

# Multiplication by numbers
function *(b::T, a::TaylorModelN) where {T<:NumberNotSeries}
    pol = a.pol * b
    rem = b * remainder(a)
    return TaylorModelN(pol, rem, expansion_point(a), domain(a))
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
# Square
function TS.square(a::TaylorModelN)
    @assert 2*get_order(a) ≤ get_order()
    res = TS.square(a.pol)
    Δ = remainder_square(a, get_order(res))
    return TaylorModelN(res, Δ, expansion_point(a), domain(a))
end

function remainder_square(a::TaylorModelN{N,T,S},
        order::Int) where {N,T,S}
    a_order = get_order(a)
    rnegl_order = 2*a_order
    aux = centered_dom(a)
    # Remaining terms of the product
    vv = Array{HomogeneousPolynomial{TS.numtype(a.pol)}}(undef, rnegl_order-order)
    suma = Array{promote_type(TS.numtype(a.pol),
                    TS.numtype(domain(a)))}(undef, rnegl_order-order)
    for k in order+1:rnegl_order
        vv[k-order] = HomogeneousPolynomial(zero(TS.numtype(a.pol)), k)
        @inbounds for i = 0:k
            (i > a_order || k-i > a_order) && continue
            TS.mul!(vv[k-order], a.pol[i], a.pol[k-i])
        end
        suma[k-order] = vv[k-order](aux)
    end
    # Bound for the neglected part of the product of polynomials
    Δnegl = sum( suma ) # = sum( sort!(suma, by=abs2) )
    Δa = a.pol(aux)
    a_rem = remainder(a)
    Δ = Δnegl + interval(2) * Δa * a_rem + a_rem^interval(2)
    Δ1 = Δnegl + (interval(2) * Δa + a_rem) * a_rem
    return _intersect_reminders((Δ1, Δ))
end

function ^(a::TaylorModelN{N,T,S}, r::Number) where {N,T,S}
    isinteger(r) && return a^Int(r)
    return rpa(x->x^r, a)
end
function ^(a::TaylorModelN{N,T,S}, n::Integer) where {N,T,S}
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return TS.square(a)
    n > 0 && return Base.power_by_squaring(a, n)
    return rpa(x->x^n, a)
end

function Base.power_by_squaring(x::TaylorModelN, p::Integer)
    @assert p ≥ 0
    (p == 0) && return one(x)
    (p == 1) && return TaylorModelN(
        Taylor1(x.pol.coeffs[:], x.pol.order), x.rem, x.x0[:], x.dom[:])
    (p == 2) && return TS.square(x)
    (p == 3) && return x*TS.square(x)
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x = TS.square(x)
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) ≥ 0
            x = TS.square(x)
        end
        y *= x
    end
    return y
end

# Operations involving TaylorModel1{TaylorModelN}
for f in (:+, :-)
    @eval function ($f)(a::TaylorModel1{TaylorModelN{N,R,S},S},
            b::TaylorModel1{TaylorModelN{N,R1,S1},S1}) where {N,R,S,R1,S1}
        a, b = fixorder(a,b)
        c = Taylor1(zero(a[0])+zero(b[0]), get_order(a))
        for i in eachindex(a)
            c[i] = ($f)(a[i], b[i])
        end
        return TaylorModel1(c, ($f)(a.rem, b.rem), expansion_point(a), domain(a))
    end
    @eval function ($f)(b::T, a::TaylorModel1{TaylorModelN{N,R,S}}) where
            {N,R,S,T<:NumberNotSeries}
        c = Taylor1(zero(($f)(b, a[0])), get_order(a))
        for i in eachindex(a)
            c[i] = ($f)(b, a[i])
        end
        return TaylorModel1(c, ($f)(zero(interval(b)),a.rem), expansion_point(a), domain(a))
    end
    @eval function ($f)(a::TaylorModel1{TaylorModelN{N,R,S}}, b::T) where
            {N,R,S,T<:NumberNotSeries}
        c = Taylor1(zero(($f)(a[0]), b), get_order(a))
        for i in eachindex(a)
            c[i] = ($f)(a[i], b)
        end
        return TaylorModel1(c, ($f)(a.rem), expansion_point(a), domain(a))
    end
end

function *(a::TaylorModel1{TaylorModelN{N,T,S},S}, b::TaylorModelN{N,T,S}) where {N,T,S}
    res = polynomial(a)*b
    return TaylorModel1(res, remainder(a), expansion_point(a), domain(a))
end
*(b::TaylorModelN{N,T,S}, a::TaylorModel1{TaylorModelN{N,T,S},S}) where {N,T,S} = a * b

# Multiplication by numbers
*(a::Taylor1{TaylorModelN{N,R,S}}, b::T) where {N,R,S,T<:NumberNotSeries} = Taylor1(a.coeffs .* b)
*(b::T, a::Taylor1{TaylorModelN{N,R,S}}) where {N,R,S,T<:NumberNotSeries} = a * b
*(a::Taylor1{TaylorModelN{N,S,R}}, b::TaylorModelN{N,S,R}) where {N,S,R} =
    Taylor1(a.coeffs .* b)
*(b::TaylorModelN{N,S,R}, a::Taylor1{TaylorModelN{N,S,R}}) where {N,S,R} = a * b

#
function TS.mul!(res::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::Taylor1{TaylorModelN{N,T,S}}, ordT::Int) where {N,T,S}
    # Sanity
    TS.zero!(res[ordT])
    tmp = res[ordT]
    for k in 0:ordT
        # res[ordT] += a[k] * b[ordT-k]
        tmp += a[k] * b[ordT-k]
        # for ordQ in eachindex(a[ordT])
        #     TS.mul!(res[ordT], a[k], b[ordT-k], ordQ)
        # end
    end
    res[ordT] = tmp
    return nothing
end

function TS.div!(c::Taylor1{TaylorModelN{N,T,S}}, a::Taylor1{TaylorModelN{N,T,S}},
        b::Taylor1{TaylorModelN{N,T,S}}, k::Int) where {N,T,S}
    # order and coefficient of first factorized term
    # ordfact, cdivfact = divfactorization(a, b)
    anz = findfirst(a)
    bnz = findfirst(b)
    anz = anz ≥ 0 ? anz : a.order
    bnz = bnz ≥ 0 ? bnz : a.order
    ordfact = min(anz, bnz)
    # Is the polynomial factorizable?
    iszero(b[ordfact]) && throw( ArgumentError(
        """Division does not define a Taylor1 polynomial;
        order k=$(ordfact) => coeff[$(ordfact)]=$(cdivfact).""") )
    TS.zero!(c, k)
    if k == 0
        c[0] = a[ordfact]/b[ordfact]
        # TS.div!(c[0], a[ordfact], b[ordfact])
        return nothing
    end
    imin = max(0, k+ordfact-b.order)
    tmp = c[imin] * b[k+ordfact-imin]
    # TS.mul!(c[k], c[imin], b[k+ordfact-imin])
    for i = imin+1:k-1
        tmp += c[i] * b[k+ordfact-i]
        # TS.mul!(c[k], c[i], b[k+ordfact-i])
    end
    if k+ordfact ≤ b.order
        tmp = (a[k+ordfact]-tmp)
        # for l in eachindex(c[k])
        #     TS.subst!(c[k], a[k+ordfact], c[k], l)
        # end
        c[k] = tmp / b[ordfact]
        # TS.div!(c[k], b[ordfact])
    else
        # c[k] = (-c[k]) / b[ordfact]
        c[k] = tmp
        TS.div_scalar!(c[k], -1, b[ordfact])
    end
    return nothing
end