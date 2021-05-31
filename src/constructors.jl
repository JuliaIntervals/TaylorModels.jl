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
            x0, dom::Interval{S}) where {T,S} = $(TM){T,S}(pol, rem, interval(x0), dom)

        # Constructor just chainging the remainder
        $(TM)(u::$(TM){T,S}, Δ::Interval{S}) where {T,S} = $(TM){T,S}(u.pol, Δ, u.x0, u.dom)

        # Short-cut for independent variable
        $(TM)(ord::Integer, x0, dom::Interval{T}) where {T} =
            $(TM)(x0 + Taylor1(eltype(x0), ord), zero(dom), interval(x0), dom)

        # Short-cut for a constructor expanding around midpoint by default
        $(TM)(ord::Integer, dom::Interval{T}) where {T} =
            $(TM)(ord, Interval(mid(dom)), dom)

        # Short-cut for a constant TM
        $(TM)(a::T, ord::Integer, x0::Interval{S}, dom::Interval{S}) where
            {T,S} = $(TM)(Taylor1([a], ord), zero(dom), x0, dom)

        # Functions to retrieve the order and remainder
        @inline get_order(tm::$TM) = get_order(tm.pol)
        @inline remainder(tm::$TM) = tm.rem
        @inline polynomial(tm::$TM) = tm.pol
        @inline domain(tm::$TM) = tm.dom
        @inline expansion_point(tm::$TM) = tm.x0
    end
end


@doc doc"""
    TaylorModel1{T,S}

Absolute Taylor model in 1 variable, providing a rigurous polynomial approximation
given by a Taylor polynomial `pol` (around `x0`) and an absolute remainder
`rem` for a function `f(x)` in one variable, valid in the interval `dom`.
This corresponds to definition 2.1.3 of Mioara Joldes' thesis.

Fields:
- `pol`: polynomial approximation, represented as `TaylorSeries.Taylor1`
- `rem`: the interval bound
- `x0` : expansion point
- `dom`: domain, interval over which the Taylor model is defined / valid

The approximation ``f(x) = p(x) + \Delta`` is satisfied for all
``x\in \mathcal{D}`` (``0\in \Delta``); `n` is the order (degree)
of the polynomial ``p(x)=\sum_{i=0}^n p_i (x - x_0)^i``.

""" TaylorModel1

@doc doc"""
    RTaylorModel1{T,S}

Relative Taylor model in 1 variable, providing a rigurous polynomial approximation
given by a Taylor polynomial `pol` (around `x0`) and a relative remainder
`rem` for a function `f(x)` in one variable, valid in the interval `dom`.
This corresponds to definition 2.3.2 of Mioara Joldes' thesis.

Fields:
- `pol`: polynomial approximation, represented as `TaylorSeries.Taylor1`
- `rem`: the interval bound
- `x0` : expansion point
- `dom`: domain, interval over which the Taylor model is defined / valid

The approximation ``f(x) = p(x) + \delta (x - x_0)^{n+1}`` is satisfied for all
``x\in \mathcal{D}``; `n` is the order (degree) of the polynomial
``p(x)=\sum_{i=0}^n p_i (x - x_0)^i``.

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
TaylorModelN(pol::TaylorN{T}, rem::Interval{S}, x0::IntervalBox{N,S},
    dom::IntervalBox{N,S}) where {N,T,S} = TaylorModelN{N,T,S}(pol, rem, x0, dom)

# Constructor for just changing the remainder
TaylorModelN(u::TaylorModelN{N,T,S}, Δ::Interval{S}) where {N,T,S} =
    TaylorModelN{N,T,S}(u.pol, Δ, u.x0, u.dom)

# Short-cut for independent variable
TaylorModelN(nv::Integer, ord::Integer, x0::IntervalBox{N,T}, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(x0[nv] + TaylorN(Interval{T}, nv, order=ord), zero(dom[1]), x0, dom)

# Short-cut for a constant
TaylorModelN(a::Interval{T}, ord::Integer, x0::IntervalBox{N,T}, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(TaylorN(a, ord), zero(dom[1]), x0, dom)
TaylorModelN(a::T, ord::Integer, x0::IntervalBox{N,T}, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(TaylorN(a, ord), zero(dom[1]), x0, dom)

# Functions to retrieve the order and remainder
@inline get_order(tm::TaylorModelN) = tm.pol.order
@inline remainder(tm::TaylorModelN) = tm.rem
@inline polynomial(tm::TaylorModelN) = tm.pol
@inline domain(tm::TaylorModelN) = tm.dom
@inline expansion_point(tm::TaylorModelN) = tm.x0


"""
    TMSol{N,T,V1,V2,M}

Structure containing the solution of a validated integration.

# Fields
    `time :: AbstractVector{T}`  Vector containing the expansion time of the `TaylorModel1` solutions

    `fp   :: AbstractVector{IntervalBox{N,T}}`  IntervalBox vector representing the flowpipe

    `xTMv :: AbstractMatrix{TaylorModel1{TaylorN{T},T}}`  Matrix whose entry `xTMv[i,t]` represents 
    the `TaylorModel1` of the i-th dependent variable, obtained at time time[t].
"""
struct TMSol{N,T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{IntervalBox{N,T}},
        M<:AbstractMatrix{TaylorModel1{TaylorN{T},T}}}
    time :: V1
    fp   :: V2
    xTM  :: M

    function TMSol(time::V1, fp::V2, xTM::M) where 
            {N,T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{IntervalBox{N,T}},
            M<:AbstractMatrix{TaylorModel1{TaylorN{T},T}}}
        @assert length(time) == length(fp) == size(xTM,2) && N == size(xTM,1)
        return new{N,T,V1,V2,M}(time, fp, xTM)
    end
end

get_time(a::TMSol) = getfield(a,:time)
get_fp(a::TMSol)   = getfield(a,:fp)
get_xTM(a::TMSol)  = getfield(a,:xTM)
get_time(a::TMSol, n::Int) = getindex(getfield(a,:time),n)
get_fp(a::TMSol, n::Int)   = getindex(getfield(a,:fp),n)
get_xTM(a::TMSol, n::Int)  = getindex(getfield(a,:xTM),:,n)
