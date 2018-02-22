module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
using IntervalRootFinding
using RecipesBase


const Interval = IntervalArithmetic.Interval
setformat(:full)

import Base: setindex!,
    ==, +, -, *, /, ^,
    zero, findfirst, #iszero,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

import TaylorSeries: integrate, get_order, evaluate


export Taylor1Model, bound, make_Taylor_model, TMcomposition,
        taylor1_var, integrate, degree,
        calculate_set, Taylor_step

export TMAbsRem, TMRelRem, Interval, remainder,
    rpa, rpafp, boundarem, boundrrem


include("Taylor1/Taylor1.jl")
include("TMs/constructors.jl")
include("TMs/bounds.jl")
include("TMs/rpa_functions.jl")
include("TMs/arithmetic.jl")
include("TMs/recipe.jl")


end # module
