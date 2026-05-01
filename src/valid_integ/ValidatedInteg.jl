# Submodule related to validated integration methods

module ValidatedInteg

using TaylorModels
using Reexport
@reexport using TaylorIntegration
using LinearAlgebra: cond, mul!
using Parameters
using StaticArrays

export TMSol, TMSol3, flowpipe, get_xTM
export shrink_wrapping!, absorb_remainder
export validated_integ, validated_integ2, validated_integ3
export iscontractive, picard_lindelof, picard_lindelof!
export mince_in_time

const TI = TaylorIntegration
const IA = IntervalArithmetic
const IANumTypes = IA.NumTypes

include("cache.jl")
include("TMSol.jl")
include("integ_utils.jl")
include("validated_integ.jl")
include("validated_integ2.jl")
include("validated_integ3.jl")
include("recipe.jl")


"""
    iscontractive(Δ, Δx)

Checks if `Δ .⊆ Δx` is satisfied.
"""
iscontractive(Δ::Interval{T}, Δx::Interval{T}) where {T} =
    issubset_interval(Δ, Δx)
iscontractive(Δ::AbstractVector{Interval{T}}, Δx::AbstractVector{Interval{T}}) where
    {T} = all(iscontractive.(Δ[:], Δx[:]))


"""
    normalize_taylorNs!

Convert the `x0::Vector{Interval{T}}` to a `q0::Vector{TaylorN{T}}` inplace. Evaluating
each component of `q0` in the `interval(-1,1)` yields the corresponding component
of `x0`.
"""
function normalize_taylorNs!(q0::Vector{TaylorN{T}}, x0::Vector{Interval{T}}) where {T}
    @inbounds for ind in eachindex(q0)
        q0[ind][0][1] = mid(x0[ind])
        q0[ind][1][ind] = radius(x0[ind])
    end
    return q0
end

end # module