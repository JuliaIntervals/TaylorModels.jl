# auxiliary.jl

# getindex, fistindex, lastindex
for TM in (:TaylorModel1, :RTaylorModel1, :TaylorModelN)
    @eval begin
        copy(f::$TM) = $TM(copy(f.pol), f.rem, f.x0, f.dom)
        @inline firstindex(a::$TM) = 0
        @inline lastindex(a::$TM) = get_order(a)

        getindex(a::$TM, n::Integer) = a.pol[n]
        getindex(a::$TM, u::UnitRange) = a.pol[u]
        getindex(a::$TM, c::Colon) = a.pol[c]
        # getindex(a::$TM, u::StepRange{Int,Int}) = a.pol[u[:]]

        constant_term(a::$TM) = constant_term(polynomial(a))
        linear_polynomial(a::$TM) = linear_polynomial(polynomial(a))

        # norm
        norm(x::$TM, p::Real=2) = norm(polynomial(x), p)
    end
end

# setindex, iscontained
for TM in tupleTMs
    @eval begin
        setindex!(a::$TM{T,S}, x::T, n::Integer) where {T<:Number, S} = a.pol[n] = x
        setindex!(a::$TM{T,S}, x::T, u::UnitRange) where {T<:Number, S} = a.pol[u] .= x
        function setindex!(a::$TM{T,S}, x::Array{T,1}, u::UnitRange) where {T<:Number, S}
            @assert length(u) == length(x)
            for ind in eachindex(x)
                a.pol[u[ind]] = x[ind]
            end
        end
        setindex!(a::$TM{T,S}, x::T, c::Colon) where {T<:Number, S} = a[c] .= x
        setindex!(a::$TM{T,S}, x::Array{T,1}, c::Colon) where {T<:Number, S} = a[c] .= x

        iscontained(a, tm::$TM) = a ∈ domain(tm)-tm.x0
        iscontained(a::Interval, tm::$TM) = a ⊆ domain(tm)-tm.x0
    end
end
iscontained(a, tm::TaylorModelN) = a ∈ domain(tm)-tm.x0
iscontained(a::IntervalBox, tm::TaylorModelN) = a ⊆ domain(tm)-tm.x0


# fixorder and bound_truncation
for TM in tupleTMs
    @eval begin
        function fixorder(a::$TM, b::$TM)
            @assert tmdata(a) == tmdata(b)
            a.pol.order == b.pol.order && return a, b

            order = min(a.pol.order, b.pol.order)
            apol0, bpol0 = polynomial.((a, b))
            apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

            # Bound for the neglected part of the polynomial
            dom = a.dom - a.x0
            Δa = bound_truncation($TM, apol0, dom, order) + remainder(a)
            Δb = bound_truncation($TM, bpol0, dom, order) + remainder(b)

            return $TM(apol, Δa, a.x0, a.dom), $TM(bpol, Δb, b.x0, b.dom)
        end

        function bound_truncation(::Type{$TM}, a::Taylor1, aux::Interval,
                order::Int)
            order ≥ get_order(a) && return zero(aux)
            if $TM == TaylorModel1
                res = Taylor1(copy(a.coeffs))
                res[0:order] .= zero(res[0])
            else
                res = Taylor1(copy(a.coeffs[order+2:end]), get_order(a)-order )
            end
            return res(aux)
        end

    end
end

function fixorder(a::TaylorModelN, b::TaylorModelN)
    @assert tmdata(a) == tmdata(b)
    a.pol.order == b.pol.order && return a, b

    order = min(a.pol.order, b.pol.order)
    apol0, bpol0 = polynomial.((a, b))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Bound for the neglected part of the polynomial
    dom = a.dom - a.x0
    Δa = bound_truncation(TaylorModelN, apol0, dom, order) + remainder(a)
    Δb = bound_truncation(TaylorModelN, bpol0, dom, order) + remainder(b)

    return TaylorModelN(apol, Δa, a.x0, a.dom), TaylorModelN(bpol, Δb, b.x0, b.dom)
end

function bound_truncation(::Type{TaylorModelN}, a::TaylorN, aux::IntervalBox,
        order::Int)
    order ≥ get_order(a) && return zero(aux[1])
    res = deepcopy(a)
    res[0:order] .= zero(res[0])
    return res(aux)
end
