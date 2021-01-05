using TaylorModels

t = TaylorModel1(3, interval(0), -0.5..0.5)

@show texp = exp(t)

@show tsin = sin(t)

s = texp + tsin
p = texp * tsin

@show TaylorModels.bound_taylor1(t)
@show TaylorModels.bound_taylor1(p)
@show TaylorModels.bound_taylor1(s)

#= Output
texp = exp(t) =  [1, 1] + [1, 1] t + [0.5, 0.5] t² + [0.166666, 0.166667] t³ + [-0, 0.00288794]
tsin = sin(t) =  [1, 1] t + [-0.166667, -0.166666] t³ + [-0.00124851, 0.00124851]
TaylorModels.bound_taylor1(t) = [-0.5, 0.5]
TaylorModels.bound_taylor1(p) = [-0.291667, 0.791667]
TaylorModels.bound_taylor1(s) = [0.124999, 2.12501]
[0.124999, 2.12501]
=#
