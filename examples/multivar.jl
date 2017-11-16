using TaylorModels, IntervalArithmetic, TaylorSeries

a, b = set_variables("a b", order=3)

u0 = Taylor1([1 + a], 3)

t = TaylorModel(3, 0..0, -0.1..0.1, u0, 0..0, [-0.05..0.05, -0.05..05])

bound(t)


# SOLVE ODE  ẋ = x:

using TaylorModels, IntervalArithmetic, TaylorSeries

n = 10  # order
a, b = set_variables("a b", order=n)
bounds = [-0.05..0.05, -0.05..0.05]  # bounds on a and b

u0 = Taylor1([1 + a], n)   # initial condition as function of a, b

# Find (non-rigorous) Taylor series expansion by Picard iteration:
∫ = integrate
u = u0
for i in 1:n+1   # how many iterations are required?
    u = u0 + ∫(u)   # solving ẋ = x
end

@assert u == u0 + ∫(u)   # fixed point of non-rigorous Picard

t_interval = 0..0.1

# Prepare a TaylorModel with 0 error term to estimate size of required interval
t = TaylorModel(n, 0..0, t_interval, u, 0..0, bounds)

t_new = ∫(t)  # just so that t_new exists outside the loop


# Search for interval giving contraction, and hence rigorous Picard iteration

for i in 1:10
    t_new = ∫(t, u0[0])
    @show t.Δ, t_new.Δ, t_new.Δ ⊆ t.Δ

    (t_new.Δ ⊆ t.Δ) && break  # if this happens, then have found contraction

    Δ = 2*t_new.Δ  # otherwise keep increasing the size of the error bound interval
    t = TaylorModel(n, 0..0, t_interval, u, Δ, bounds)
end


# Iterate contraction until reach fixed point:

for i in 1:20
    t_new = ∫(t, u0[0])
    @show t.Δ, t_new.Δ, t_new.Δ ⊆ t.Δ, t_new.Δ == t.Δ

    Δ = t_new.Δ
    t = TaylorModel(n, 0..0, t_interval, u, Δ, bounds)
end

# Find initial set for next step:
tf = (t_interval.hi)..(t_interval.hi)
U1 = t_new(tf)  # start from initial


## VOLTERRA EQNS (Makino & Berz, suppression of wrapping effect)
# ẋ = 2x(1 - y)
# ẏ = -y(1 - x)

# Init conds:
# x0 ∈ 1 + [-0.05, 0.05]
# y0 ∈ 3 + [-0.05, 0.05]


using TaylorModels, IntervalArithmetic, TaylorSeries

n = 6  # order
a, b = set_variables("a b", order=n)

h = 0.01
bounds = (-h..h) * ones(2)  # bounds on a and b

u0 = Taylor1([(1..1) + (1..1)*a], n)   # initial condition as function of a, b
v0 = Taylor1([(3..3) + (3..3)*b], n)   # initial condition as function of a, b

∫ = integrate
u = u0
v = v0
u_new = u0   # so exist outside loop
v_new = v0

for i in 1:n+1   # how many iterations are required?
    u_new = u0 + ∫(   2  * u * (1 - v) )
    v_new = v0 + ∫(       -v * (1 - u) )

    u, v = u_new, v_new
end


t_interval = 0..0.1
uu = TaylorModel(n, 0..0, t_interval, u, 0..0, bounds)
vv = TaylorModel(n, 0..0, t_interval, v, 0..0, bounds)

uu_new = ∫( 2 * uu * (1 - vv), u0[0] )
vv_new = ∫(    -vv * (1 - uu), v0[0] )
