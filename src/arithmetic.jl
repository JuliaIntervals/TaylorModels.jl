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

        # Multiplication
        function *(a::$TM, b::$TM)
            @assert tmdata(a) == tmdata(b)
            a_order = get_order(a)
            b_order = get_order(b)
            rnegl_order = a_order + b_order
            aux = a.dom - a.x0

            # Returned polynomial
            res = a.pol * b.pol
            order = res.order

            # Remainder of the product
            if $TM == TaylorModel1

                # Remaing terms of the product as reduced Taylor1 (factored polynomial)
                rnegl = Taylor1(zero(res[0]), rnegl_order)
                for k in order+1:rnegl_order
                    @inbounds for i = 0:k
                        (i > a_order || k-i > b_order) && continue
                        rnegl[k] += a.pol[i] * b.pol[k-i]
                    end
                end

                # Bound for the neglected part of the product of polynomials
                Δnegl = rnegl(aux)
                Δ = remainder_product(a, b, aux, Δnegl)

            else

                # Remaing terms of the product as reduced Taylor1 (factored polynomial)
                rnegl = Taylor1(zero(res[0]), rnegl_order-order)
                for k in order+1:rnegl_order
                    @inbounds for i = 0:k
                        (i > a_order || k-i > b_order) && continue
                        rnegl[k-order-1] += a.pol[i] * b.pol[k-i]
                    end
                end
                Δnegl = rnegl(aux)
                Δ = remainder_product(a, b, aux, Δnegl, order)
            end

            return $TM(res, Δ, a.x0, a.dom)
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
function remainder_product(a, b, aux, Δnegl)
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem
    return Δ
end
function remainder_product(a::TaylorModel1{TaylorModelN{N,T,S},S},
        b::TaylorModel1{TaylorModelN{N,T,S},S}, aux, Δnegl) where {N,T,S}
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem

    # Evaluate at the TMN centered domain
    auxN = a[0].dom - a[0].x0
    ΔN = Δ(auxN)
    return ΔN
end
function remainder_product(a::RTaylorModel1, b::RTaylorModel1, aux, Δnegl, order)
    Δa = a.pol(aux)
    Δb = b.pol(aux)
    V = aux^(order+1)
    Δ = Δnegl + Δb * a.rem + Δa * b.rem + a.rem * b.rem * V
    return Δ
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
    a_order = get_order(a)
    b_order = get_order(b)
    rnegl_order = a_order + b_order
    @assert rnegl_order ≤ get_order()
    aux = a.dom - a.x0

    # Returned polynomial
    res = a.pol * b.pol
    order = res.order

    # Remaing terms of the product
    vv = Array{HomogeneousPolynomial{eltype(res)}}(undef, rnegl_order-order)
    suma = Array{promote_type(eltype(res), eltype(a.dom))}(undef, rnegl_order-order)
    for k in order+1:rnegl_order
        vv[k-order] = HomogeneousPolynomial(zero(eltype(res)), k)
        @inbounds for i = 0:k
            (i > a_order || k-i > b_order) && continue
            TaylorSeries.mul!(vv[k-order], a.pol[i], b.pol[k-i])
        end
        suma[k-order] = vv[k-order](aux)
    end

    # Bound for the neglected part of the product of polynomials
    Δnegl = sum( sort!(suma, by=abs2) )
    Δ = remainder_product(a, b, aux, Δnegl)

    return TaylorModelN(res, Δ, a.x0, a.dom)
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
function ^(a::TaylorModelN{N,T,S}, r::Number) where {N,T,S}
    r == 0 && return one(a)
    r == 1 && return a
    r == 2 && return a*a
    return rpa(x->x^r, a)
end
function ^(a::TaylorModelN{N,T,S}, n::Integer) where {N,T,S}
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return a*a
    return rpa(x->x^n, a)
end
