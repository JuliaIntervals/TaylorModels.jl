# constructors.jl

const tupleTMs = (:TM1AbsRem, :TM1RelRem)
const NumberNotSeries = TaylorSeries.NumberNotSeries

#=
Structs `TM1AbsRem{T}` and `TM1RelRem{T}` are essentially identical, except
the way the remainder is computed and that the remainder for `TM1AbsRem{T}`
must contain 0.
=#
for TM in tupleTMs
    @eval begin
        struct $(TM){T,S} <: AbstractSeries{T}
            pol  :: Taylor1{T}    # polynomial approx (of order `ord`)
            rem  :: Interval{S}   # remainder
            x0   :: Interval{S}   # expansion point
            iI   :: Interval{S}   # interval of interest

            # Inner constructor
            function $(TM){T,S}(pol::Taylor1{T}, rem::Interval{S},
                    x0::Interval{S}, iI::Interval{S}) where {T,S}
                if $(TM) == TM1AbsRem
                    @assert zero(S) ∈ rem && x0 ⊆ iI
                else
                    @assert x0 ⊆ iI
                end
                return new{T,S}(pol, rem, x0, iI)
            end
        end

        # Outer constructors
        $(TM)(pol::Taylor1{T}, rem::Interval{S},
            x0::Interval{S}, iI::Interval{S}) where {T,S} = $(TM){T,S}(pol, rem, x0, iI)

        # Short-cut for independent variable
        $(TM)(ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
            $(TM)(x0 + Taylor1(Interval{T}, ord), zero(iI), x0, iI)

        # Short-cut for a constant
        $(TM)(a::Interval{T}, ord::Int, x0::Interval{T},
            iI::Interval{T}) where {T} = $(TM)(Taylor1([a], ord), zero(iI), x0, iI)

        # Functions to retrieve the order and remainder
        get_order(tm::$TM) = get_order(tm.pol)
        remainder(tm::$TM) = tm.rem
        polynomial(tm::$TM) = tm.pol
    end
end


@doc """
    TM1AbsRem{T,S}

Taylor Models with Absolute Remainder. Corresponds to definition 2.1.3
(Mioara Joldes thesis). `T` is the type of the coefficients of the polynomial,
ans `S` is the type of the `Interval`s.

""" TM1AbsRem

@doc """
    TM1RelRem{T,S}

Taylor Models with Relative Remainder. Corresponds to definition 2.3.2
(Mioara Joldes thesis). `T` is the type of the coefficients of the polynomial,
ans `S` is the type of the `Interval`s.

""" TM1RelRem


# TMNAbsRem's struct
"""
    TMNAbsRem{N,T,S}

Taylor Models with Absolute Remainder for `N` independent variables.
`T` is the type of the `TaylorN` polynomial, and `S` the type of
the `IntervalBox`es.

"""
struct TMNAbsRem{N,T,S} <: AbstractSeries{T}
    pol  :: TaylorN{T}        # polynomial approx (of order `ord`)
    rem  :: Interval{S}       # remainder
    x0   :: IntervalBox{N,S}  # expansion point
    iI   :: IntervalBox{N,S}  # interval of interest

    # Inner constructor
    function TMNAbsRem{N,T,S}(pol::TaylorN{T}, rem::Interval{S},
            x0::IntervalBox{N,S}, iI::IntervalBox{N,S}) where {N,T<:NumberNotSeries,S<:Real}

        @assert N == get_numvars()
        @assert zero(S) ∈ rem && x0 ⊆ iI

        return new{N,T,S}(pol, rem, x0, iI)
    end
end

# Outer constructors
TMNAbsRem(pol::TaylorN{T}, rem::Interval{S}, x0::IntervalBox{N,S}, iI::IntervalBox{N,S}) where {N,T,S} =
    TMNAbsRem{N,T,S}(pol, rem, x0, iI)

# Short-cut for independent variable
TMNAbsRem(nv::Int, ord::Int, x0::IntervalBox{N,T}, iI::IntervalBox{N,T}) where {N,T} =
    TMNAbsRem(x0[nv] + TaylorN(Interval{T}, nv, order=ord), zero(iI[1]), x0, iI)

# Short-cut for a constant
TMNAbsRem(a::Interval{T}, ord::Int, x0::IntervalBox{N,T}, iI::IntervalBox{N,T}) where {N,T} =
    TMNAbsRem(TaylorN(a, ord), zero(iI[1]), x0, iI)

# Functions to retrieve the order and remainder
get_order(tm::TMNAbsRem) = tm.pol.order
remainder(tm::TMNAbsRem) = tm.rem
polynomial(tm::TMNAbsRem) = tm.pol
