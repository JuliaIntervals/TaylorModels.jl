module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
# Use current master (until next tag) of IntervalRootFinding
using IntervalRootFinding#; Pkg.checkout("IntervalRootFinding")
using RecipesBase


import Base: setindex!,
    ==, +, -, *, /, ^,
    zero, one, findfirst, #iszero,
    promote, show,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

import TaylorSeries: integrate, get_order, evaluate,
    pretty_print


export Taylor1Model, bound, make_Taylor_model, TMcomposition,
        taylor1_var, integrate, degree,
        calculate_set, Taylor_step

export TM1AbsRem, TM1RelRem, TMNAbsRem,
    remainder, polynomial,
    rpa, fp_rpa, boundarem, boundrrem

export TaylorNModel



include("TMs/constructors.jl")
include("TMs/promotion.jl")
include("TMs/bounds.jl")
include("TMs/evaluate.jl")
include("TMs/rpa_functions.jl")
include("TMs/arithmetic.jl")
include("TMs/integration.jl")
include("TMs/recipe.jl")
include("TMs/show.jl")

include("Taylor1/Taylor1.jl")
include("TaylorN/TaylorN.jl")


end # module
