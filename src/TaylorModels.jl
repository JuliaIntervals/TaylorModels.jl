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

doc"""
A `TaylorModel` represents a (single-variable) polynomial approximation to a function $f(t)$. The coefficients of the polynomial
may be `TaylorN`.

The fields are:
- `n`: degree of the polynomial
- `x0`: expansion point
- `I`: interval over which the Taylor model is defined / valid
- `p`: the polynomial, represented as `TaylorSeries.Taylor1`
- `Δ`: the interval bound
- `bounds`: an array of intervals representing the bounds of the variables that occur in the coefficients.
"""
struct TaylorModel{T,S}
    n::Int      # degree
    x0::Interval{T}  # expansion point
    I::Interval{T}   # interval over which the Taylor model is valid
    p::Taylor1{S}  # Taylor Taylor1nomial
    Δ::Interval{T}   # interval remainder bound
    bounds::Vector{Interval{T}}
end

include("arithmetic.jl")
include("functions.jl")
include("bound.jl")

include("integrate.jl")
include("draw.jl")



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

"""
Evaluate a TaylorModel at a point
"""
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
