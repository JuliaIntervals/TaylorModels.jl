

# missing methods in TaylorSeries.jl:

degree(f::TaylorN) = f.order
setindex!(f::TaylorN, x, i) = f.coeffs[i+1] = x

doc"""
A `Taylor1Model` represents a polynomial approximation to a function $f(t)$ of a single variable.

The fields are:
- `n`: degree of the polynomial
- `x0`: expansion point as IntervalBox
- `I`: IntervalBox over which the Taylor model is defined / valid
- `p`: the polynomial, of type `TaylorSeries.TaylorN`
- `Δ`: the interval bound
"""
struct TaylorNModel{N,T,S}
    n::Int      # degree
    x0::IntervalBox{N,T}  # expansion point
    I::IntervalBox{N,T}   # interval over which the Taylor model is valid
    p::TaylorN{S}  # Taylor Taylor1nomial
    Δ::Interval{T}   # interval remainder bound
    order_bounds::Vector{Interval{T}}
end

# constructors that calculate order bounds automatically:
function TaylorNModel(n, x0, I::IntervalBox{N,T}, p, Δ) where {N,T}
    order_bounds = [homog([I...]) for homog in p.coeffs]  # bound homogeneous polynomials

    return TaylorNModel(n, x0, I, p, Δ, order_bounds)
end

TaylorNModel(n, x0, I::IntervalBox{N,T}, p) where {N,T} =  TaylorNModel(n, x0, I, p, Interval{T}(0))

include("arithmetic.jl")
# include("functions.jl")
# include("bound.jl")
#
include("integrate.jl")
# include("draw.jl")

import Base.copy
copy(f::TaylorNModel) = TaylorNModel(f.n, f.x0, f.I, copy(f.p), f.Δ, copy(f.order_bounds))


# Taylor1Model for a constant:
# Taylor1Model(n::Int, x0, I, c::T) where {T<:AbstractFloat} = Taylor1Model{T}(n, x0, I, Taylor1{Interval{T}}(c), Interval{T}(0), [])
#
# Taylor1Model(n::Int, x0, I) = Taylor1Model{Float64}(n, x0, I, Taylor1{Interval{Float64}}([0.0, 1.0]), Interval{Float64}(0.0))
#

#Taylor1Model(n::Int, x0::Interval{T}, I::Interval{T}, p::Taylor1{S}, Δ::Interval{T}) where {T, S} = Taylor1Model{T, S}(n, x0, I, p, Δ)

# Taylor1Model for a variable:
