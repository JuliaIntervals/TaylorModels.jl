module TaylorModels

using IntervalArithmetic, TaylorSeries
using RecipesBase

const Interval = IntervalArithmetic.Interval
import TaylorSeries.integrate

import Base: exp, sin, inv, cos, identity, +, *, /, ^, -

export TaylorModel, bound, make_Taylor_model, TMcomposition, taylor_var,
        integrate, degree


import Base: setindex!

degree(f::Taylor1) = f.order
setindex!(f::Taylor1, i, x) = f.coeffs[i+1] = x

abstract type AbstractTaylorModel end

# the bounds field contains interval bounds for variables in Taylor model whose
# cofficients are multivar polynomials
struct TaylorModel{T,S} <: AbstractTaylorModel
    n::Int      # degree
    x0::Interval{T}  # expansion point
    I::Interval{T}   # interval over which the Taylor model is valid
    p::Taylor1{S}  # Taylor Taylor1nomial
    Δ::Interval{T}   # interval remainder bound
    bounds::Vector{Interval{T}}
end


doc"""
Compute a rigorous bound for a polynomial over an interval.
"""
function bound(f::TaylorModel)
    x0, I, p = f.x0, f.I, f.p

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + bound(p[i], f.bounds)
    end

    return B
end

function bound(p::Taylor1, x0, I)

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + p[i]
    end

    return B
end

function bound(p::Taylor1, x0, I, bounds)

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + bound(p[i], bounds)
    end

    return B
end


bound(f::TaylorN, bounds) = evaluate(f, bounds)  # can replace by better polynomial bounder


doc"""
Compute a rigorous bound for a TaylorModel.
"""
# bound(f::TaylorModel) = bound(f.p, f.x0, f.I)

bound(x::Interval, bounds) = x

# TaylorModel for a constant:
# TaylorModel(n::Int, x0, I, c::T) where {T<:AbstractFloat} = TaylorModel{T}(n, x0, I, Taylor1{Interval{T}}(c), Interval{T}(0), [])
#
# TaylorModel(n::Int, x0, I) = TaylorModel{Float64}(n, x0, I, Taylor1{Interval{Float64}}([0.0, 1.0]), Interval{Float64}(0.0))
#

#TaylorModel(n::Int, x0::Interval{T}, I::Interval{T}, p::Taylor1{S}, Δ::Interval{T}) where {T, S} = TaylorModel{T, S}(n, x0, I, p, Δ)

# TaylorModel for a variable:
taylor_var(n::Int, x0, I) = TaylorModel(n, Interval(x0), I, Taylor1{Interval{Float64}}(Interval{Float64}[0.0, 1.0], n), Interval{Float64}(0.0), Interval{Float64}[])

# assumes f and g are expansions around the same point x0 with the same order

import Base.copy
copy(f::TaylorModel) = TaylorModel(f.n, f.x0, f.I, copy(f.p), f.Δ, f.bounds)

function +(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    return TaylorModel(f.n, f.x0, f.I, f.p + g.p, f.Δ + g.Δ, f.bounds)
end

function -(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    return TaylorModel(f.n, f.x0, f.I, f.p - g.p, f.Δ - g.Δ, f.bounds)
end

function +(α::Real, f::TaylorModel)
    g = copy(f)

    g.p[0] += (α..α)

    return g
end

-(f::TaylorModel) = (-1 * f)
-(α::Real, f::TaylorModel) = α + (-f)


## Taylor coeffcients of different functions
## TODO: Change nth_deriv -> taylor_coeff
nth_deriv(::typeof(exp), i, x) = exp(x) / factorial(i)

function nth_deriv(::typeof(sin), i, x)
    j = i % 4

    (j == 0) && (return sin(x) / factorial(i))
    (j == 1) && (return cos(x) / factorial(i))
    (j == 2) && (return -sin(x) / factorial(i))
    (j == 3) && (return -cos(x) / factorial(i))
end

nth_deriv(::typeof(inv), i, x) = (-1)^i / (x^(i+1))

doc"""
Make a TaylorModel for a given function over a given domain.
"""
function make_Taylor_model(f, n, x0, I::Interval{T}, bounds) where T

    a = zeros(typeof(I), n+1)

    for i in 0:n
        a[i+1] = nth_deriv(f, i, x0)
    end

    p = Taylor1(a)

    Γ = nth_deriv(f, n+1, I)

    if sup(Γ) ≤ 0 || inf(Γ) ≥ 0

        lo = Interval(inf(I))
        hi = Interval(sup(I))

        Δlo = f(lo) - bound(p, x0, lo)
        Δhi = f(hi) - bound(p, x0, hi)
        Δx0 = f(x0) - bound(p, x0, x0)
        Δ = Interval(min(inf(Δlo), inf(Δx0), inf(Δhi)), max(sup(Δlo), sup(Δx0), sup(Δhi)))

    else
        V = (I - x0)^(n+1)
        Δ = V * Γ
    end

    return TaylorModel(n, x0, I, p, Δ, bounds)
end

*(α::Real, f::TaylorModel) = TaylorModel(f.n, f.x0, f.I, Taylor1(α*f.p), (α..α)*f.Δ, f.bounds)




function *(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I
    @assert f.bounds == g.bounds

    n, x0, I = f.n, f.x0, f.I


    c = zeros(eltype(g.p), 2n+1)

    a = f.p
    b = g.p

    for i in 0:n
        for j in 0:n
            c[i+j+1] += a[i] * b[j]
        end
    end

    d = zeros(eltype(g.p), 2n+1)

    d[n+2:end] = c[n+2:end]

    B = bound(Taylor1(d), x0, I, f.bounds)
    Bf = bound(a, x0, I, f.bounds)
    Bg = bound(b, x0, I, f.bounds)   # assuming g.bounds == f.bounds

    Δ = B + (f.Δ * Bg) + (g.Δ * Bf) + (f.Δ * g.Δ)

    return TaylorModel(n, f.x0, f.I, Taylor1(c[1:n+1]), Δ, f.bounds)

end


import Base.zero
zero(::Type{TaylorModel{T}}, n, x0, I::Interval{T}, bounds) where {T<:AbstractFloat} = TaylorModel(n, x0, I, Taylor1{Interval{Float64}}(zeros(n+1)), Interval{T}(0), bounds)


doc"""
Evaluate a polynomial of a TaylorModel.
`b` are the coefficients of the Taylor1nomial.
"""
function poly_eval_of_TM(b, f, I::Interval{T}, x0, n, bounds) where {T<:AbstractFloat}

    #M = TaylorModel{T}(n, x0, I, Taylor1(zeros(n+1)), Interval{T}(0))
    M = zero(TaylorModel{T}, n, x0, I, bounds)

    for i in n:-1:0
        M *= f
        M.p[0] += b[i]  # add constant
    end

    return M
end

doc"""
Calculate the TaylorModel of `(g∘f)` given a function `g` and a TaylorModel `f`.
"""
function TMcomposition(g, f::TaylorModel)
    x0, I, n = f.x0, f.I, f.n

    Bf = bound(f.p, x0, I)
    a, Δf = f.p, f.Δ

    Mg = make_Taylor_model(g, n, a[0], Bf + Δf, f.bounds)

    b, Δg = Mg.p, Mg.Δ

    a[0] = 0  # zero out first element

    M1 = TaylorModel(n, x0, I, a, Δf, f.bounds)
    M = poly_eval_of_TM(b, M1, I, x0, n, f.bounds)

    c, Δ = M.p, M.Δ

    return TaylorModel(n, x0, I, c, Δ+Δg, f.bounds)
end


# allow exp(t) etc. for a TaylorModel t:
for g in (:exp, :sin, :inv)
    @eval $g(f::TaylorModel) = TMcomposition($g, f)
end


function /(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    M = f * inv(g)
end

doc"""
Integrate a TaylorModel.
`x0` is optional constant to add.
"""
function integrate(f::TaylorModel, x0=0)

    p2 = integrate(f.p)

    high_order_term = f.p[end]
    Δ = integral_bound(f)

    t = TaylorModel(f.n, f.x0, f.I, p2, Δ, f.bounds)
    t.p[0] = x0  # constant term

    return t

end

function integral_bound(f::TaylorModel)
    n = degree(f.p)
    high_order_term = f.p[n]

    coeff = bound(high_order_term, f.bounds)
    power = (f.I - f.x0)^n

    (coeff * power + f.Δ) * diam(f.I)
end

# evaluate a TaylorModel at a point:
# (f::TaylorModel{T})(x) where {T<:AbstractFloat} = (f.p)(Interval{T}(x)) + f.Δ
function (f::TaylorModel)(x)
    if x in f.I
        return (f.p)(Interval(x)) + f.Δ
    else
        throw(ArgumentError("Cannot evaluate TaylorModel at point $x outside interval of definition $(f.I)"))
    end
end


# plot recipe for plotting 1D TaylorModels
@recipe function g(f::TaylorModel)

    x0, I, n, p, Δ = f.x0, f.I, f.n, f.p, f.Δ

    alpha --> 0.5
    seriestype := :shape

    xs = linspace(I.lo, I.hi, 100)
    evals = f.(xs)

    ylos = [y.lo for y in evals]
    yhis = [y.hi for y in evals]

    xs = [xs; reverse(xs); xs[1]]
    ys = [ylos; reverse(yhis); ylos[1]]

    xs, ys
end



end
