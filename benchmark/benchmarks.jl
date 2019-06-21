using BenchmarkTools, TaylorModels

SUITE = BenchmarkGroup()

include("arithmetic.jl")
include("daisy.jl")
