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
    end
end

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
    end
end


# fixorder
function fixorder(a::TaylorModel1, b::TaylorModel1)
    @assert tmdata(a) == tmdata(b)
    a.pol.order == b.pol.order && return a, b

    order = min(a.pol.order, b.pol.order)
    apol0, bpol0 = polynomial.((a, b))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Bound for the neglected part of the polynomial
    dom = a.dom - a.x0
    Δa = bound_truncation(TaylorModel1, apol0, dom, order) + remainder(a)
    Δb = bound_truncation(TaylorModel1, bpol0, dom, order) + remainder(b)

    return TaylorModel1(apol, Δa, a.x0, a.dom), TaylorModel1(bpol, Δb, b.x0, b.dom)
end
function fixorder(a::RTaylorModel1, b::RTaylorModel1)
    @assert tmdata(a) == tmdata(b)
    a.pol.order == b.pol.order && return a, b

    order = min(a.pol.order, b.pol.order)
    apol0, bpol0 = polynomial.((a, b))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Bound for the neglected part of the polynomial
    dom = a.dom - a.x0
    Δa = bound_truncation(RTaylorModel1, apol0, dom, order) + remainder(a)
    Δb = bound_truncation(RTaylorModel1, bpol0, dom, order) + remainder(b)

    return RTaylorModel1(apol, Δa, a.x0, a.dom), RTaylorModel1(bpol, Δb, b.x0, b.dom)
end
function fixorder(a::TaylorModelN, b::TaylorModelN)
    @assert tmdata(a) == tmdata(b)
    a.pol.order == b.pol.order && return a, b

    order = min(a.pol.order, b.pol.order)
    apol0, bpol0 = polynomial.((a, b))
    @show(typeof(apol0))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Bound for the neglected part of the polynomial
    dom = a.dom - a.x0
    Δa = bound_truncation(TaylorModelN, apol0, dom, order) + remainder(a)
    Δb = bound_truncation(TaylorModelN, bpol0, dom, order) + remainder(b)

    return TaylorModelN(apol, Δa, a.x0, a.dom), TaylorModelN(bpol, Δb, b.x0, b.dom)
end


# bound_truncation
function bound_truncation(::Type{TaylorModel1}, a::Taylor1, aux::Interval,
        order::Int)
    order ≥ get_order(a) && return zero(aux)
    res = Taylor1(copy(a.coeffs))
    res[0:order] .= zero(res[0])
    return res(aux)
end
function bound_truncation(::Type{RTaylorModel1}, a::Taylor1, aux::Interval,
        order::Int)
    order ≥ get_order(a) && return zero(aux)
    res = Taylor1(copy(a.coeffs[order+2:end]), get_order(a)-order )
    return res(aux)
end
function bound_truncation(::Type{TaylorModelN}, a::TaylorN, aux::IntervalBox,
        order::Int)
    order ≥ get_order(a) && return zero(aux[1])
    res = deepcopy(a)
    res[0:order] .= zero(res[0])
    return res(aux)
end
