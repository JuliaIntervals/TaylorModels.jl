module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
@reexport using TaylorIntegration

using IntervalRootFinding
using LinearAlgebra: norm, mul!, cond, isposdef

using RecipesBase

using Markdown

import Base: setindex!, getindex, copy, firstindex, lastindex,
    ==, +, -, *, /, ^,
    zero, one, findfirst, #iszero,
    promote, show,
    inv, sqrt, exp, log, sin, cos, tan,
    asin, acos, atan, sinh, cosh, tanh

using TaylorSeries: derivative, ∇

import TaylorSeries: integrate, get_order, evaluate, pretty_print,
    constant_term, linear_polynomial, fixorder, get_numvars

import IntervalArithmetic: showfull

import LinearAlgebra: norm

# export Taylor1Model, bound, make_Taylor_model, TMcomposition,
#         taylor1_var, integrate, degree,
#         calculate_set, Taylor_step

export TaylorModel1, RTaylorModel1, TaylorModelN, TMSol

export remainder, polynomial, domain, expansion_point, flowpipe, get_xTM,
    rpa, fp_rpa, bound_remainder,
    validated_integ, validated_integ2

export linear_dominated_bounder, quadratic_fast_bounder

include("constructors.jl")
include("auxiliary.jl")
include("promotion.jl")
include("bounds.jl")
include("evaluate.jl")
include("rpa_functions.jl")
include("arithmetic.jl")
include("integration.jl")
include("show.jl")
include("validatedODEs.jl")
include("recipe.jl")

# include("Taylor1/Taylor1.jl")
# include("TaylorN/TaylorN.jl")


end # module
