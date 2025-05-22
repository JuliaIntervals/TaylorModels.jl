module TaylorModels

using Reexport
@reexport using TaylorSeries, IntervalArithmetic
# @reexport using TaylorIntegration

using IntervalRootFinding
using LinearAlgebra: norm, isposdef

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
    constant_term, linear_polynomial, nonlinear_polynomial,
    fixorder, get_numvars

import LinearAlgebra: norm

export TaylorModel1, RTaylorModel1, TaylorModelN, TMSol

export remainder, polynomial, domain, expansion_point, flowpipe, get_xTM,
    rpa, fp_rpa, bound_remainder, centered_dom, symmetric_box

export linear_dominated_bounder, quadratic_fast_bounder

setdisplay(:full)

include("constructors.jl")
include("auxiliary.jl")
include("promotion.jl")
include("bounds.jl")
include("evaluate.jl")
include("arithmetic.jl")
include("rpa_functions.jl")
include("integration.jl")
include("show.jl")
include("valid_integ/ValidatedInteg.jl")
include("recipe.jl")


end # module
