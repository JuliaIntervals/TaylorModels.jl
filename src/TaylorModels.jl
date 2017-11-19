module TaylorModels

using IntervalArithmetic, TaylorSeries
using RecipesBase

const Interval = IntervalArithmetic.Interval
import TaylorSeries.integrate

import Base: exp, sin, inv, cos, identity, +, *, /, ^, -

export TaylorModel, bound, make_Taylor_model, TMcomposition,
        taylor_var, integrate, degree,
        calculate_set, Taylor_step


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

include("arithmetic.jl")
include("functions.jl")

include("integrate.jl")
include("draw.jl")



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
function (f::TaylorModel)(t)
    if t in f.I
        return (f.p)(t - f.x0) + f.Δ
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
