# Submodule related to validated integration methods

module ValidatedInteg

using TaylorModels
using Reexport
@reexport using TaylorIntegration
using LinearAlgebra: cond, mul!
using Parameters
using StaticArrays

export shrink_wrapping!, absorb_remainder
export validated_integ, validated_integ2, validated_integ3
export iscontractive, picard_lindelof, picard_lindelof!

const TI = TaylorIntegration
const IA = IntervalArithmetic
const IANumTypes = IA.NumTypes

include("cache.jl")
include("validated_integ.jl")
include("validated_integ2.jl")
include("validated_integ3.jl")


"""
    iscontractive(Δ, Δx)

Checks if `Δ .⊆ Δx` is satisfied.
"""
iscontractive(Δ::Interval{T}, Δx::Interval{T}) where {T} =
    issubset_interval(Δ, Δx)
iscontractive(Δ::AbstractVector{Interval{T}}, Δx::AbstractVector{Interval{T}}) where
    {T} = all(iscontractive.(Δ[:], Δx[:]))


function normalize_taylorNs!(q0::Vector{TaylorN{T}}, x0::Vector{Interval{T}},
        orderQ::Int) where {T}
    @inbounds for ind in eachindex(q0)
        mx0 = mid(x0[ind])
        q0[ind] = mx0 + TaylorN(ind, order=orderQ) * radius(x0[ind])
    end
    return q0
end

end # module