
include("TaylorModels.jl")
using TaylorModels, IntervalArithmetic


t1 = make_Taylor_model(exp, 3, 0, -0.5..0.5)
t2 = make_Taylor_model(sin, 3, 0, -0.5..0.5)

t1 + t2
t1 * t2

t = taylor_var(3, 0, -0.5..0.5)
compute_bound(t)

t = t1 + t2
compute_bound(t)
