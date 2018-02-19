
# missing methods in TaylorSeries.jl:

degree(f::Taylor1) = f.order
setindex!(f::Taylor1, x, i) = f.coeffs[i+1] = x

doc"""
A `Taylor1Model` is a Taylor model in 1 variable, providing a bounding box for a function as a polynomial $p$ plus an interval bound $\Delta$ for a function $f(x)$ of a single variable.

Fields:
- `n`: degree of the polynomial
- `x0`: expansion point
- `I`: interval over which the Taylor model is defined / valid
- `p`: the polynomial, represented as `TaylorSeries.Taylor1`
- `Δ`: the interval bound

The approximation is $f(x) \simeq \sum_i p_i (x - x_0)^i$.
"""
struct Taylor1Model{T,S}
    n::Int      # degree
    x0::Interval{T}  # expansion point
    I::Interval{T}   # interval over which the Taylor model is valid
    p::Taylor1{S}  # Taylor Taylor1nomial
    Δ::Interval{T}   # interval remainder bound
end

Taylor1Model(n, x0, I::Interval{T}, p, Δ) where{T} = Taylor1Model(n, interval(x0), I, Taylor1{Interval{T}}(p), interval(Δ))

# convenience constructors with same n, x0, I:
Taylor1Model(f, p, Δ) = Taylor1Model(f.n, f.x0, f.I, p, Δ)
Taylor1Model(f, Δ) = Taylor1Model(f, f.p, Δ)


include("arithmetic.jl")
include("functions.jl")
include("bound.jl")

include("integrate.jl")
include("draw.jl")


import Base.show
function show(io::IO, f::Taylor1Model)
    print(io,
    """Taylor1 model of degree $(f.n):
     - x0:  $(f.x0)
     -  I:  $(f.I)
     -  p: $(f.p)
     -  Δ:  $(f.Δ)
    """
    )
end


zero(f::Taylor1Model{T}) where {T} = Taylor1Model(f, zeros(length(f.p)), 0)

degree(f::Taylor1Model) = degree(f.p)


# Taylor1Model for a variable:
taylor1_var(n::Int, x0, I::T) where {T} = Taylor1Model(n, Interval(x0), I, Taylor1(T[x0, 1], n), T(0))

taylor1_var(n::Int, I::Interval) = taylor1_var(n, mid(I), I)


import Base.copy
copy(f::Taylor1Model) = Taylor1Model(f, copy(f.p), f.Δ)

"""
Evaluate a Taylor1Model at a point
"""
function (f::Taylor1Model)(x)

    x0 = f.x0

    if x ∈ f.I
        return (f.p)(x - x0) + f.Δ
    else
        throw(ArgumentError("Cannot evaluate Taylor1Model at point $x outside interval of definition $(f.I)"))
    end
end

"Note that TaylorSeries indices start at 0."
Base.getindex(f::Taylor1Model, i::Integer) = f.p[i]


# plot recipe for plotting 1D Taylor1Models
@recipe function g(f::Taylor1Model)

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
