module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
using IntervalRootFinding

const Interval = IntervalArithmetic.Interval
setformat(:full)

import Base: ==, +, -, *, /, ^,
    zero, #iszero,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

import TaylorSeries: get_order, evaluate

export TMAbsRem, TMRelRem, Interval, remainder,
    rpa, rpafp, bound_arem, bound_rrem

include("constructors.jl")
include("bounds.jl")
include("rpa_functions.jl")
include("arithmetic.jl")


end # module
