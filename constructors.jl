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

# Short-cut for independent variable
TMAbsRem(ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
    TMAbsRem(x0 + Taylor1(Interval{T}, ord), zero(iI), x0, iI)

# Short-cut for a constant
TMAbsRem(a::Interval{T}, ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
    TMAbsRem(Taylor1([a], ord), zero(iI), x0, iI)



"""
    TMRelRem{T,S}

Taylor Models with Relative Remainder. Corresponds to definition 2.1.3
(Mioara Joldes thesis).

"""
struct TMRelRem{T}
    pol  :: Taylor1{Interval{T}}    # polynomial approx (of order `ord+1`)
    # rrem :: Interval{T}             # relative remainder
    x0   :: Interval{T}             # expansion point
    iI   :: Interval{T}             # interval of interest

    # Inner constructor
    function TMRelRem{T}(pol::Taylor1{Interval{T}}, # rrem::Interval{T},
            x0::Interval{T}, iI::Interval{T}) where {T}
        @assert zero(T) ∈ pol[end] && x0 ⊆ iI
        return new{T}(pol, x0, iI)
    end
end

# Outer constructors
TMRelRem(pol::Taylor1{Interval{T}}, # rrem::Interval{T},
    x0::Interval{T}, iI::Interval{T}) where {T} = TMRelRem{T}(pol, x0, iI)

# Short-cut for independent variable
function TMRelRem(ord::Int, x0::Interval{T}, iI::Interval{T}) where {T}
    TMRelRem(x0 + Taylor1(Interval{T}, ord+1), x0, iI)
end

# Short-cut for a constant
TMRelRem(a::Interval{T}, ord::Int, x0::Interval{T}, iI::Interval{T}) where {T} =
    TMRelRem(Taylor1([a], ord+1), x0, iI)



# Functions to inspect the entities of `TMAbsRem` and `TMRelRem`
get_order(tm::TMAbsRem) = tm.pol.order
get_order(tm::TMRelRem) = tm.pol.order-1

remainder(tm::TMAbsRem) = tm.arem
remainder(tm::TMRelRem) = tm.pol[end]
