# auxiliary.jl

# getindex, fistindex, lastindex
for TM in (:TaylorModel1, :RTaylorModel1, :TaylorModelN)
    @eval begin
        copy(f::$TM) = $TM(copy(f.pol), remainder(f), expansion_point(f), domain(f))
        @inline firstindex(a::$TM) = 0
        @inline lastindex(a::$TM) = get_order(a)

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
iscontained(a::Vector{Interval{S}}, tm::TaylorModelN) where {S} =
    all(issubset_interval.(a, centered_dom(tm)))

symmetric_box(T::Type{S}) where {S<:IntervalArithmetic.NumTypes} = [interval(-one(T), one(T)) for _ = 1:get_numvars()]
symmetric_box(T,N) = [interval(-one(T), one(T)) for _ = 1:get_numvars()]


# fixorder and bound_truncation
for TM in tupleTMs
    @eval begin
        function fixorder(a::$TM, b::$TM)
            @assert tmdata(a) == tmdata(b)
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

function bound_truncation(::Type{TaylorModel1}, a::Taylor1{TaylorN{T}}, aux::Interval,
        order::Int) where {T}
    order ≥ get_order(a) && return zero(aux)
    # Assumes that the domain for the TaylorN variables is the symmetric normalized box -1 .. 1
    symIbox = symmetric_box(T)
    res = Taylor1(evaluate.(a.coeffs, Ref(symIbox)))
    res[0:order] .= zero(res[0])
    return res(aux)
end


function fixorder(a::TaylorModelN, b::TaylorModelN)
    @assert tmdata(a) == tmdata(b)
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

function bound_truncation(::Type{TaylorModelN}, a::TaylorN, aux::Vector{Interval},
        order::Int)
    order ≥ get_order(a) && return zero(aux[1])
    res = deepcopy(a)
    res[0:order] .= zero(res[0])
    return res(aux)
end


# TMSol utilities
@inline firstindex(a::TMSol) = firstindex(a.time)
@inline lastindex(a::TMSol)  = lastindex(a.time)
@inline Base.length(a::TMSol) = length(a.time)
@inline Base.iterate(a::TMSol, state=0) = state ≥ lastindex(a) ? nothing : (a[state+1], state+1)
@inline Base.eachindex(a::TMSol) = firstindex(a):lastindex(a)

getindex(a::TMSol, n::Integer) = a.xTM[:,n]
getindex(a::TMSol, u::UnitRange) = a.xTM[:,u]
getindex(a::TMSol, c::Colon) = a.xTM[:,c]
getindex(a::TMSol, n::Integer, m::Integer) = a.xTM[m,n]
getindex(a::TMSol, c::Colon, m::Integer) = a.xTM[m,c]
getindex(a::TMSol, u::UnitRange, m::Integer) = a.xTM[m,u]
