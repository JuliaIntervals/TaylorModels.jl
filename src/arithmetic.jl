# arithmetic.jl

# Addition, substraction and other functions
for TM in tupleTMs
    @eval begin
        zero(a::$TM) = $TM(zero(a.pol), zero(remainder(a)), expansion_point(a), domain(a))
        one(a::$TM) = $TM(one(a.pol), zero(remainder(a)), expansion_point(a), domain(a))

        # iszero(a::$TM) = iszero(a.pol) && iszero(zero(remainder(a)))

        findfirst(a::$TM) = findfirst(a.pol)

        ==(a::$TM, b::$TM) =
            a.pol == b.pol &&  isequal_interval(remainder(a), remainder(b)) &&
                all(isequal_interval.(tmdata(a), tmdata(b)))
                # expansion_point(a) == expansion_point(b) && domain(a) == domain(b)


        # Addition
        +(a::$TM) = $TM(a.pol, remainder(a), expansion_point(a), domain(a))

        function +(a::$TM, b::$TM)
            a, b = fixorder(a, b)
            return $TM(a.pol+b.pol, remainder(a)+remainder(b), expansion_point(a), domain(a))
        end

        +(a::$TM, b::T) where {T<:NumberNotSeries} = $TM(a.pol+b, remainder(a),
            expansion_point(a), domain(a))

        +(b::T, a::$TM) where {T<:NumberNotSeries} = $TM(b+a.pol, remainder(a),
            expansion_point(a), domain(a))


        # Substraction
        -(a::$TM) = $TM(-a.pol, -remainder(a), expansion_point(a), domain(a))

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
        function *(a::$TM, b::$TM)
            @assert all(isequal_interval.(tmdata(a), tmdata(b)))
            a_order = get_order(a)
            b_order = get_order(b)
            rnegl_order = a_order + b_order
            aux = centered_dom(a)

            # Returned polynomial
            a_pol = polynomial(a)
            b_pol = polynomial(b)
            res = a_pol * b_pol
            order = get_order(res)

            # Remainder of the product
            if $TM == TaylorModel1

                # Remaining terms of the product as reduced Taylor1 (factored polynomial)
                rnegl = Taylor1(zero(res[0]), rnegl_order)
                for k in order+1:rnegl_order
                    @inbounds for i = 0:k
                        (i > a_order || k-i > b_order) && continue
                        rnegl[k] += a_pol[i] * b_pol[k-i]
                    end
                end

                # Bound for the neglected part of the product of polynomials
                Δnegl = rnegl(aux)
                Δ = remainder_product(a, b, aux, Δnegl)
            else

                # Remaining terms of the product as reduced Taylor1 (factored polynomial)
                rnegl = Taylor1(zero(res[0]), rnegl_order-order)
                for k in order+1:rnegl_order
                    @inbounds for i = 0:k
                        (i > a_order || k-i > b_order) && continue
                        rnegl[k-order-1] += a_pol[i] * b_pol[k-i]
                    end
                end
                Δnegl = rnegl(aux)
                Δ = remainder_product(a, b, aux, Δnegl, order)
            end

            return $TM(res, Δ, expansion_point(a), domain(a))
        end

        # Division by numbers
        /(a::$TM, b::T) where {T<:NumberNotSeries} = a * inv(b)

        /(b::T, a::$TM) where {T<:NumberNotSeries} = b * inv(a)

        # Power
        function ^(a::$TM, r::Number)
            r == zero(r) && return one(a)
            r == 1 && return a
            r == 2 && return a*a
            return rpa(x->x^r, a)
        end

        function ^(a::$TM, n::Integer)
            n == 0 && return one(a)
            n == 1 && return a
            n == 2 && return a*a
            return rpa(x->x^n, a)
        end
    end
end

# Remainder of the product
# TaylorModel1
function remainder_product(a, b, aux, Δnegl)
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl + Δb * a_rem + Δa * b_rem + a_rem * b_rem
    return Δ
end
function remainder_product(a::TaylorModel1{TaylorN{T}, S},
                           b::TaylorModel1{TaylorN{T}, S},
                           auxT, Δnegl) where {T, S}
    # An N-dimensional symmetrical IntervalBox is assumed
    # to bound the TaylorN part
    auxQ = symmetric_box(S)
    Δa = a.pol(auxT)(auxQ)
    Δb = b.pol(auxT)(auxQ)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl(auxQ) + Δb * a_rem + Δa * b_rem + a_rem * b_rem
    return Δ
end
function remainder_product(a::TaylorModel1{TaylorModelN{T,S},S},
        b::TaylorModel1{TaylorModelN{T,S},S}, aux, Δnegl) where {T,S}
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
function remainder_product(a, b, aux, Δnegl, order)
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    V = aux^(order+1)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V
    return Δ
end
function remainder_product(a::RTaylorModel1{TaylorN{T},S}, b::RTaylorModel1{TaylorN{T},S},
                            aux, Δnegl, order) where {T, S}
    # N = get_numvars()
    # An N-dimensional symmetrical IntervalBox is assumed
    # to bound the TaylorN part
    auxQ = symmetric_box(T)
    Δa = a.pol(aux)(auxQ)
    Δb = b.pol(aux)(auxQ)
    V = aux^(order+1)
    a_rem = remainder(a)
    b_rem = remainder(b)
    Δ = Δnegl(auxQ) + Δb * a_rem + Δa * b_rem + a_rem * b_rem * V
    return Δ
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
    Δ = Δnegl + remainder(a) * (aux)^(order-m)
    return RTaylorModel1( apol, Δ, expansion_point(a), domain(a) )
end


# Same as above, for TaylorModelN
zero(a::TaylorModelN) = TaylorModelN(zero(a.pol), zero(remainder(a)),
    expansion_point(a), domain(a))
one(a::TaylorModelN) = TaylorModelN(one(a.pol), zero(remainder(a)),
    expansion_point(a), domain(a))

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

        $(op)(b::T, a::TaylorModelN) where {T<:NumberNotSeries} = TaylorModelN($(
            op)(b, a.pol), $(op)(remainder(a)), expansion_point(a), domain(a))
    end
end


# Multiplication
function *(a::TaylorModelN, b::TaylorModelN)
    @assert all(isequal_interval.(tmdata(a), tmdata(b)))
    a_order = get_order(a)
    b_order = get_order(b)
    rnegl_order = a_order + b_order
    @assert rnegl_order ≤ get_order()
    aux = centered_dom(a)

    # Returned polynomial
    a_pol = polynomial(a)
    b_pol = polynomial(b)
    res = a_pol * b_pol
    order = get_order(res)

    # Remaing terms of the product
    vv = Array{HomogeneousPolynomial{TS.numtype(res)}}(undef, rnegl_order-order)
    suma = Array{promote_type(TS.numtype(res),
                    TS.numtype(domain(a)))}(undef, rnegl_order-order)
    for k in order+1:rnegl_order
        vv[k-order] = HomogeneousPolynomial(zero(TS.numtype(res)), k)
        @inbounds for i = 0:k
            (i > a_order || k-i > b_order) && continue
            TaylorSeries.mul!(vv[k-order], a_pol[i], b_pol[k-i])
        end
        suma[k-order] = vv[k-order](aux)
    end

    # Bound for the neglected part of the product of polynomials
    Δnegl = sum( suma ) # = sum( sort!(suma, by=abs2) )
    Δ = remainder_product(a, b, aux, Δnegl)

    return TaylorModelN(res, Δ, expansion_point(a), domain(a))
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
function ^(a::TaylorModelN{T,S}, r::Number) where {T,S}
    r == 0 && return one(a)
    r == 1 && return a
    r == 2 && return a*a
    return rpa(x->x^r, a)
end
function ^(a::TaylorModelN{T,S}, n::Integer) where {T,S}
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return a*a
    return rpa(x->x^n, a)
end
