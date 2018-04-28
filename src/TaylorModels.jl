module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
using IntervalRootFinding
using RecipesBase


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

export TM1AbsRem, TM1RelRem, remainder,
    rpa, rpafp, boundarem, boundrrem


include("TMs/constructors.jl")
include("TMs/bounds.jl")
include("TMs/rpa_functions.jl")
include("TMs/arithmetic.jl")
include("TMs/integration.jl")
include("TMs/recipe.jl")

include("Taylor1/Taylor1.jl")


end # module
