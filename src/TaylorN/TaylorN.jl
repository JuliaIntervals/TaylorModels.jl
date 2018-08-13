

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

# convenience constructors
TaylorNModel(f, p, Δ) = TaylorNModel(f.n, f.x0, f.I, p, Δ)
TaylorNModel(f, Δ) = TaylorNModel(f.n, f.x0, f.I, f.p, Δ)

include("arithmetic.jl")
# include("functions.jl")
# include("bound.jl")
#
include("integrate.jl")
# include("draw.jl")

import Base.copy
copy(f::TaylorNModel) = TaylorNModel(f.n, f.x0, f.I, copy(f.p), f.Δ, copy(f.order_bounds))

import Base.show
function show(io::IO, f::TaylorNModel)
    print(io,
    """Taylor model of degree $(f.n):
       - x0: $(f.x0)
       - I: $(f.I)
       - p: $(f.p)
       - Δ: $(f.Δ)
    """
    )
end
