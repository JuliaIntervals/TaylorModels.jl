module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
# Use current master (until next tag) of IntervalRootFinding
using IntervalRootFinding
using RecipesBase


import Base: setindex!, getindex, copy,
    ==, +, -, *, /, ^,
    zero, one, findfirst, #iszero,
    promote, show,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

import TaylorSeries: integrate, get_order, evaluate,
    pretty_print


# export Taylor1Model, bound, make_Taylor_model, TMcomposition,
#         taylor1_var, integrate, degree,
#         calculate_set, Taylor_step

export TaylorModel1, RTaylorModel1, TaylorModelN

export remainder, polynomial,
    rpa, fp_rpa, bound_absrem, bound_relrem


include("constructors.jl")
include("promotion.jl")
include("bounds.jl")
include("evaluate.jl")
include("rpa_functions.jl")
include("arithmetic.jl")
include("integration.jl")
include("recipe.jl")
include("show.jl")

# include("Taylor1/Taylor1.jl")
# include("TaylorN/TaylorN.jl")


end # module
