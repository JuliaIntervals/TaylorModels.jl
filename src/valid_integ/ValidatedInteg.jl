# Submodule related to validated integration methods

module ValidatedInteg

using TaylorModels
using Reexport
@reexport using TaylorIntegration
using LinearAlgebra: cond, mul!
using Parameters

export shrink_wrapping!, validated_integ, validated_integ2

const TI = TaylorIntegration
const IA = IntervalArithmetic
const IANumTypes = IA.NumTypes

include("tweaksTI.jl")
include("validated_integ.jl")
include("validated_integ2.jl")

end
