using TaylorModels, IntervalArithmetic, TaylorSeries

t, a, b = set_variables(Interval{Float64}, "t a b", order=3)

u00 = 1+a   # initial condition

II = -0.1..0.1
IIbox = IntervalBox(II, II, II)

u0 = TaylorNModel(3, IntervalBox(0..0, 0..0, 0..0), IIbox, u00, 0..0)


# SOLVE ODE  ẋ = f(x):

f(x) = x

n = 3
∫_dt = u -> integrate(u, 1)
u = u0
for i in 1:n+5   # how many iterations are required?
    u = u0 + ∫_dt(f(u))
end

u

# u == u0 + ∫(u,1)   # fixed point of non-rigorous Picard
# need to define equality

integrate(u, 1)

t_interval = 0..0.1

# Prepare a TaylorModel with 0 error term to estimate size of required interval
w = TaylorNModel(u.n, u.x0, u.I, u.p, 0..0)

w_new = ∫(w, 1)  # just so that t_new exists outside the loop

# Search for interval giving contraction, and hence rigorous Picard iteration

for i in 1:10
    w_new = u0 + ∫_dt(w)
    @show w.Δ, w_new.Δ, w_new.Δ ⊆ w.Δ

    (w_new.Δ ⊆ w.Δ) && break  # if this happens, then have found contraction

    Δ = 2*w_new.Δ  # otherwise keep increasing the size of the error bound interval
    w = TaylorNModel(w.n, w.x0, w.I, w.p, Δ)
end

w.Δ
w_new.Δ

w.p

evaluate # Iterate contraction until reach fixed point:

for i in 1:20
    w_new = u0 + ∫_dt(w)
    @show w.Δ, w_new.Δ, w_new.Δ ⊆ w.Δ, w_new.Δ == w.Δ

    Δ = w_new.Δ
    w = TaylorNModel(w.n, w.x0, w.I, w.p, Δ)
end

w
w.Δ
w.p([t_interval.hi..t_interval.hi, a, b]) + w.Δ

(1+a)*exp(0.1..0.1)

all( ((1+a)*exp(tt..tt))([IIbox...]) ⊆ ( w.p([tt..tt, IIbox[2], IIbox[3]]) + w.Δ ) for tt in 0:0.01:0.1 )

for tt in 0:0.01:0.1
    println(((1+a)*exp(tt..tt))([IIbox...]), "\t", ( w.p([tt..tt, IIbox[2], IIbox[3]]) + w.Δ ))
end

(1+a)*exp(0.01)
w.p([0.01, a, b]) + w.Δ
# Find initial set for next step:
tf = (t_interval.hi)..(t_interval.hi)d
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
v0 = Taylor1([(3..3) + (1..1)*b], n)   # initial condition as function of a, b

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

uu, vv = uu_new, vv_new

uu_new = ∫( 2 * uu * (1 - vv), u0[0] )
vv_new = ∫(    -vv * (1 - uu), v0[0] )

# repeat the following, doubling each time



uu = TaylorModel(n, 0..0, t_interval, u, uu_new.Δ, bounds)
vv = TaylorModel(n, 0..0, t_interval, u, vv_new.Δ, bounds)

uu_new = ∫( 2 * uu * (1 - vv), u0[0] )
vv_new = ∫(    -vv * (1 - uu), v0[0] )

while ! ((uu_new.Δ ⊆ uu.Δ) && (vv_new.Δ ⊆ vv.Δ))
    uu = TaylorModel(n, 0..0, t_interval, u, 2*uu_new.Δ, bounds)
    vv = TaylorModel(n, 0..0, t_interval, u, 2*vv_new.Δ, bounds)

    uu_new = ∫( 2 * uu * (1 - vv), u0[0] )
    vv_new = ∫(    -vv * (1 - uu), v0[0] )
end



# only need to bound the integral, not actually carry out
# the whole integral

# contract:

for i in 1:10
   uu, vv = uu_new, vv_new
   uu_new = ∫( 2 * uu * (1 - vv), u0[0] )
   vv_new = ∫(    -vv * (1 - uu), v0[0] )
   @show uu_new.Δ, vv_new.Δ
   @show uu_new.Δ ⊆ uu.Δ, vv_new.Δ ⊆ vv.Δ

end

U = uu_new(t_interval.hi)
V = vv_new(t_interval.hi)


function calculate_set(U::TaylorN, V::TaylorN, bounds)

    xs = []
    ys = []

    # iterate over 4 sides of initial box

    num_points = 50

    b = bounds[2].lo
    for a in linspace(bounds[1].lo, bounds[2].hi, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    a = bounds[1].hi
    for b in linspace(bounds[2].lo, bounds[2].hi, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    b = bounds[2].hi
    for a in linspace(bounds[2].hi, bounds[2].lo, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    a = bounds[1].lo
    for b in linspace(bounds[2].hi, bounds[2].lo, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    return mid.(xs), mid.(ys)


end
