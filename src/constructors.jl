# constructors.jl

const tupleTMs = (:TaylorModel1, :RTaylorModel1)
const NumberNotSeries = TaylorSeries.NumberNotSeries

#=
Structs `TaylorModel1` and `RTaylorModel1` are essentially identical, except
the way the remainder is handled; the remainder for `TaylorModel1`
must contain 0.
=#
for TM in tupleTMs
    @eval begin
        struct $(TM){T,S} <: AbstractSeries{T}
            pol  :: Taylor1{T}    # polynomial approx (of order `ord`)
            rem  :: Interval{S}   # remainder
            x0   :: Interval{S}   # expansion point
            dom  :: Interval{S}   # interval of interest

            # Inner constructor
            function $(TM){T,S}(pol::Taylor1{T}, rem::Interval{S},
                    x0::Interval{S}, dom::Interval{S}) where {T,S}
                if $(TM) == TaylorModel1
                    @assert zero(S) ∈ rem && x0 ⊆ dom
                else
                    @assert x0 ⊆ dom
                end
                return new{T,S}(pol, rem, x0, dom)
            end
        end

        # Outer constructors
        $(TM)(pol::Taylor1{T}, rem::Interval{S},
            x0::Interval{S}, dom::Interval{S}) where {T,S} = $(TM){T,S}(pol, rem, x0, dom)

        # Short-cut for independent variable
        $(TM)(ord::Integer, x0, dom::Interval{T}) where {T} =
            $(TM)(x0 + Taylor1(eltype(x0), ord), zero(dom), interval(x0), dom)

        # Short-cut for a constructor expanding around midpoint by default
        $(TM)(ord::Integer, dom::Interval{T}) where {T} =
            $(TM)(ord, Interval(mid(dom)), dom)

        # Short-cut for a constant TM
        $(TM)(a::Interval{T}, ord::Integer, x0::Interval{T},
            dom::Interval{T}) where {T} = $(TM)(Taylor1([a], ord), zero(dom), x0, dom)

        # convenience constructors with same n, x0, I:
        # TaylorModel1(f, p, Δ) = TaylorModel1(f.n, f.x0, f.dom, p, Δ)
        # TaylorModel1(f, Δ) = TaylorModel1(f, f.p, Δ)

        # Functions to retrieve the order and remainder
        get_order(tm::$TM) = get_order(tm.pol)
        remainder(tm::$TM) = tm.rem
        polynomial(tm::$TM) = tm.pol
        domain(tm::$TM) = tm.dom
    end
end


@doc doc"""
    TaylorModel1{T,S}

Taylor model in 1 variable, providing a rigurous polynomial approximation
(around `x_0`) and an absolute remainder `\Delta` for a function `f(x)` in one variable,
valid in the interval `dom`. Corresponds to definition 2.1.3 of
Mioara Joldes' thesis.

Fields:
- `pol`: polynomial approximation, represented as `TaylorSeries.Taylor1`
- `rem`: the interval bound
- `x0` : expansion point
- `dom`: domain, interval over which the Taylor model is defined / valid

The approximation `f(x) = \sum_{i=0}^n p_i (x - x_0)^i + \Delta` is
satisfied for all `x\in dom` (`0\in\Delta`); `n` is the order (degree)
of the polynomial `p(x)`.

""" TaylorModel1

@doc doc"""
    RTaylorModel1{T,S}

Taylor model in 1 variable, providing a rigurous polynomial approximation
(around `x_0`) and a relative remainder `\delta` for a function `f(x)` in one variable,
valid in the interval `dom`. Corresponds to definition 2.3.2 of
Mioara Joldes' thesis.

Fields:
- `pol`: polynomial approximation, represented as `TaylorSeries.Taylor1`
- `rem`: the interval bound
- `x0` : expansion point
- `dom`: domain, interval over which the Taylor model is defined / valid

The approximation `f(x) = \sum_i p_i (x - x_0)^i + \delta (x - x_0)^{n+1}` is
satisfied for all `x\in dom`; `n` is the order (degree) of the polynomial `p(x)`.

""" RTaylorModel1


# TaylorModelN's struct
"""
    TaylorModelN{N,T,S}

Taylor Models with absolute remainder for `N` independent variables.

"""
struct TaylorModelN{N,T,S} <: AbstractSeries{T}
    pol  :: TaylorN{T}        # polynomial approx (of order `ord`)
    rem  :: Interval{S}       # remainder
    x0   :: IntervalBox{N,S}  # expansion point
    dom  :: IntervalBox{N,S}  # interval of interest

    # Inner constructor
    function TaylorModelN{N,T,S}(pol::TaylorN{T}, rem::Interval{S},
            x0::IntervalBox{N,S}, dom::IntervalBox{N,S}) where {N,T<:NumberNotSeries,S<:Real}

        @assert N == get_numvars()
        @assert zero(S) ∈ rem && x0 ⊆ dom

        return new{N,T,S}(pol, rem, x0, dom)
    end
end

# Outer constructors
TaylorModelN(pol::TaylorN{T}, rem::Interval{S}, x0::IntervalBox{N,S}, dom::IntervalBox{N,S}) where {N,T,S} =
    TaylorModelN{N,T,S}(pol, rem, x0, dom)

# Short-cut for independent variable
TaylorModelN(nv::Integer, ord::Integer, x0::IntervalBox{N,T}, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(x0[nv] + TaylorN(Interval{T}, nv, order=ord), zero(dom[1]), x0, dom)

# Short-cut for a constant
TaylorModelN(a::Interval{T}, ord::Integer, x0::IntervalBox{N,T}, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(TaylorN(a, ord), zero(dom[1]), x0, dom)
TaylorModelN(a::T, ord::Integer, x0::IntervalBox{N,T}, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(TaylorN(a, ord), zero(dom[1]), x0, dom)

# Functions to retrieve the order and remainder
get_order(tm::TaylorModelN) = tm.pol.order
remainder(tm::TaylorModelN) = tm.rem
polynomial(tm::TaylorModelN) = tm.pol
domain(tm::TaylorModelN) = tm.dom
