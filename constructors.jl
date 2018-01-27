# constructors.jl

"""
    TMAbsRem{T,S}

Taylor Models with Absolute Remainder. Corresponds to definition 2.1.3
(Mioara Joldes thesis).

"""
struct TMAbsRem{T}
    pol  :: Taylor1{Interval{T}}    # polynomial approx (of order `ord`)
    arem :: Interval{T}             # absolute remainder
    x0   :: Interval{T}             # expansion point
    iI   :: Interval{T}             # interval of interest

    # Inner constructor
    function TMAbsRem{T}(pol::Taylor1{Interval{T}}, arem::Interval{T},
            x0::Interval{T}, iI::Interval{T}) where {T}
        @assert zero(T) ∈ arem && x0 ⊆ iI
        return new{T}(pol, arem, x0, iI)
    end
end

# Outer constructors
TMAbsRem(pol::Taylor1{Interval{T}}, arem::Interval{T},
    x0::Interval{T}, iI::Interval{T}) where {T} = TMAbsRem{T}(pol, arem, x0, iI)

# short-cuts for independent variable
TMAbsRem(ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
    TMAbsRem(Taylor1(Interval{T}, ord), zero(iI), x0, iI)

# TMAbsRem(T::Type, ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
#     TMAbsRem(Taylor1(Interval{T}, ord), zero(iI), x0, iI)

# short-cut for a constant
TMAbsRem(a::Interval{T}, ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
    TMAbsRem(Taylor1([a], ord), zero(iI), x0, iI)

# Functions to inspect the entities of `TMAbsRem`
get_order(tm::TMAbsRem) = tm.pol.order

remainder(tm::TMAbsRem) = tm.arem
