module TaylorModels

using IntervalArithmetic, TaylorSeries
using RecipesBase

const Interval = IntervalArithmetic.Interval
import TaylorSeries.integrate

import Base: exp, sin, inv, cos, identity, +, *, /, ^, -

export Taylor1Model, bound, make_Taylor_model, TMcomposition,
        taylor1_var, integrate, degree,
        calculate_set, Taylor_step

export TaylorNModel


import Base: setindex!

include("Taylor1/Taylor1.jl")
include("TaylorN/TaylorN.jl")

end
