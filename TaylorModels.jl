module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic

const Interval = IntervalArithmetic.Interval
setformat(:full)

import Base: ==, +, -, *, /, ^,
    zero, #iszero,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

import TaylorSeries: get_order, evaluate

export TMAbsRem, Interval, rpa, fprpa, bound_arem, remainder

include("constructors.jl")
include("bounds.jl")
include("rpa_functions.jl")
include("arithmetic.jl")


end # module
