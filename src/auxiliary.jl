# auxiliary.jl

# getindex
for TM in (:TaylorModel1, :RTaylorModel1, :TaylorModelN)
    @eval begin
        copy(f::$TM) = $TM(copy(f.pol), f.rem, f.x0, f.dom)

        getindex(a::$TM, n::Integer) = a.pol[n]
        getindex(a::$TM, u::UnitRange) = view(a.pol, u .+ 1 )
        getindex(a::$TM, c::Colon) = view(a.pol, c)

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
