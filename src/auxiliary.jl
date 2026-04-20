# auxiliary.jl

# getindex, fistindex, lastindex
for TM in (:TaylorModel1, :RTaylorModel1, :TaylorModelN)
    @eval begin
        tmdata(f::$TM) = (expansion_point(f), domain(f))

        copy(f::$TM) = $TM(copy(f.pol), remainder(f), expansion_point(f), domain(f))
        @inline firstindex(a::$TM) = 0
        @inline lastindex(a::$TM) = get_order(a)
        @inline Base.length(a::$TM) = get_order(a)+1
        @inline Base.iterate(a::$TM, state=0) = state ≥ lastindex(a) ? nothing : (a[state+1], state+1)
        @inline Base.eachindex(a::$TM) = firstindex(a):lastindex(a)

        getindex(a::$TM, n::Integer) = getindex(polynomial(a), n)
        getindex(a::$TM, u::UnitRange) = getindex(polynomial(a), u)
        getindex(a::$TM, c::Colon) = getindex(polynomial(a), c)
        # getindex(a::$TM, u::StepRange{Int,Int}) = getindex(polynomial(a), u)

        constant_term(a::$TM) = constant_term(polynomial(a))
        linear_polynomial(a::$TM) = linear_polynomial(polynomial(a))
        nonlinear_polynomial(a::$TM) = nonlinear_polynomial(polynomial(a))

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

        iscontained(a, tm::$TM) = in_interval(a, centered_dom(tm))
        iscontained(a::Interval, tm::$TM) = issubset_interval(a, centered_dom(tm))
    end
end
iscontained(a, tm::TaylorModelN) = all(in_interval.(a, centered_dom(tm)))
iscontained(a::AbstractVector{<:Interval}, tm::TaylorModelN) =
    all(issubset_interval.(a, centered_dom(tm)))
setindex!(a::TaylorModel1{TaylorModelN{N,T,S},S}, x::TaylorModelN{N,T,S},
    n::Int) where {N,T,S} = a.pol[n] = x
setindex!(a::Taylor1{TaylorModelN{N,T,S}}, x::TaylorModelN{N,T,S},
        n::Int) where {N,T,S} = setindex!(a.coeffs,
    TaylorModelN(TaylorN(x.pol.coeffs[:], get_order(x.pol)), x.rem, x.x0[:], x.dom[:]),
    n+1)


"""
    symmetric_box(N::Int, [Type{T} = Float64])
    symmetric_box(::Type{T})

Create the interval box [-1, 1]^N as a SVector, with elements of type T. If N is omitted,
it corresponds to `get_numvars()`.
"""
symmetric_box(N::Int) = fill(interval(-1.0, 1.0), SVector{N})
symmetric_box(N::Int, ::Type{S}) where {S<:IA.NumTypes} =
    fill(interval(-one(S), one(S)), SVector{N})
symmetric_box(::Type{T}) where {T<:IA.NumTypes} = symmetric_box(get_numvars(), T)


# fixorder and bound_truncation
for TM in tupleTMs
    @eval begin
        function fixorder(a::$TM, b::$TM)
            @assert all(isequal_interval.(tmdata(a), tmdata(b)))
            get_order(a) == get_order(b) && return a, b

            order = min(get_order(a), get_order(b))
            apol0, bpol0 = polynomial.((a, b))
            apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

            # Bound for the neglected part of the polynomial
            dom = centered_dom(a)
            Δa = bound_truncation($TM, apol0, dom, order) + remainder(a)
            Δb = bound_truncation($TM, bpol0, dom, order) + remainder(b)

            return $TM(apol, Δa, expansion_point(a), domain(a)),
                $TM(bpol, Δb, expansion_point(b), domain(b))
        end

        function bound_truncation(::Type{$TM}, a::Taylor1, aux::Interval,
                order::Int)
            order_a = get_order(a)
            order ≥ order_a && return zero(aux)
            if $TM == TaylorModel1
                res = Taylor1(copy(a.coeffs))
                res[0:order] .= zero(res[0])
            else
                res = Taylor1(a.coeffs[order+2:end], order_a-order )
            end
            return res(aux)
        end

    end
end

function bound_truncation(::Type{TaylorModel1}, a::Taylor1{TaylorN{T}}, aux::Interval,
        order::Int) where {T}
    order ≥ get_order(a) && return zero(aux)
    # Assumes that the domain for the TaylorN variables is the symmetric normalized box -1 .. 1
    symIbox = symmetric_box(numtype(aux))
    res = evaluate(a, symIbox)
    res[0:order] .= zero(res[0])
    return res(aux)
end


function fixorder(a::TaylorModelN, b::TaylorModelN)
    @assert all(isequal_interval.(tmdata(a), tmdata(b)))
    get_order(a) == get_order(b) && return a, b

    order = min(get_order(a), get_order(b))
    apol0, bpol0 = polynomial.((a, b))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Bound for the neglected part of the polynomial
    dom = centered_dom(a)
    Δa = bound_truncation(TaylorModelN, apol0, dom, order) + remainder(a)
    Δb = bound_truncation(TaylorModelN, bpol0, dom, order) + remainder(b)

    return TaylorModelN(apol, Δa, expansion_point(a), domain(a)),
        TaylorModelN(bpol, Δb, expansion_point(b), domain(b))
end

function bound_truncation(::Type{TaylorModelN}, a::TaylorN, aux::AbstractVector{<:Interval},
        order::Int)
    order ≥ get_order(a) && return zero(aux[1])
    res = TaylorN(a.coeffs[:], get_order(a))
    res[0:order] .= zero(res[0])
    return res(aux)
end


"""
    pol_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S}) :: TaylorModel1{Interval{S},S}

TaylorModel1 formed by the TaylorModelN remainders.
"""
function pol_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S}) where {N,T,S}
    order = get_order(tm)
    polI = Taylor1(zero(remainder(tm[0])), order)
    for k in eachindex(tm)
        polI[k] = remainder(tm[k])
    end
    return TaylorModel1(polI, remainder(tm), expansion_point(tm), domain(tm))
end


"""
    total_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S})

Computes de total reminder of a `TaylorModel1{TaylorModelN}` by
computing the polynomial associated to the `TaylorModelN` reminders,
evaluated in the centered domain, and adding the remainder of the
`TaylorModel1`.
"""
total_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S}) where {N,T,S} =
    evaluate(pol_remainder(tm), centered_dom(tm))


"""
    shift_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S})

Returns a `TaylorModel1{TaylorModelN}` with null remainder for
the `TaylorModelN` part, and the total remainder in the `TaylorModel1`.
"""
function shift_remainder(tm::TaylorModel1{TaylorModelN{N,T,S}, S}) where {N,T,S}
    rem = total_remainder(tm)
    z = interval(0.0)
    for k in eachindex(tm)
        tm[k] = TaylorModelN(tm[k], z)
    end
    return TaylorModel1(tm, rem)
end


"""
    pure_polynomial(tm)

Return the pure polynomial part of a `tm::TaylorModel1{TaylorModelN{N,T,S},S}`
as a `Taylor1{TaylorN{T}}`.
"""
function pure_polynomial(tm::TaylorModel1{TaylorModelN{N,T,S},S}) where {N,T,S}
    order = get_order(tm)
    vTN = Vector{TaylorN{T}}(undef, order+1)
    for ind in eachindex(tm)
        vTN[ind+1] = polynomial(tm[ind])
    end
    return Taylor1(vTN, order)
end
