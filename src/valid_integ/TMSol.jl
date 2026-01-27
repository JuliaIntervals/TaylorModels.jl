"""
    TMSol{T,V1,V2,M}

Structure containing the solution of a validated integration (methods
`validated_integ`, `validated_integ2`).

# Fields
    `time :: AbstractVector{T}`: Vector containing the expansion time of the
    `TaylorModel1` solutions

    `fp   :: AbstractVector{Vector{N,Interval{T}}}`: Vector of interval vectors
    representing the flowpipe

    `xTMv :: AbstractMatrix{TaylorModel1{TaylorN{T},T}}`: Matrix whose entry `xTMv[i,t]`
    represents the `TaylorModel1` of the i-th dependent variable, obtained at time time[t].
"""
struct TMSol{T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{Vector{Interval{T}}},
        M<:AbstractMatrix{TaylorModel1{TaylorN{T},T}}}
    time :: V1
    fp   :: V2
    xTM  :: M
    function TMSol(time::V1, fp::V2, xTM::M) where
            {T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{Vector{Interval{T}}},
            M<:AbstractMatrix{TaylorModel1{TaylorN{T},T}}}
        @assert length(time) == length(fp) == size(xTM,2)
        return new{T,V1,V2,M}(time, fp, xTM)
    end
end


"""
    TMSol3{N,T,V1,V2,M}

Structure containing the solution of a validated integration (method
`validated_integ3`).

# Fields
    `time :: AbstractVector{T}`: Vector containing the expansion time of the
    `TaylorModel1` solutions

    `fp   :: AbstractVector{Vector{Interval{T}}}`: Vector of interval vectors
    representing the flowpipe

    `xTMv :: AbstractMatrix{TaylorModel1{TaylorModelN{N,T,T},T}}`: Matrix whose
    entry `xTMv[i,t]` represents the `TaylorModel1` of the i-th dependent variable,
    obtained at time time[t].
"""
struct TMSol3{N,T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{Vector{Interval{T}}},
        M<:AbstractMatrix{TaylorModel1{TaylorModelN{N,T,T},T}}}
    time :: V1
    fp   :: V2
    xTM  :: M
    function TMSol3(time::V1, fp::V2, xTM::M) where
            {N,T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{Vector{Interval{T}}},
            M<:AbstractMatrix{TaylorModel1{TaylorModelN{N,T,T},T}}}
        @assert length(time) == length(fp) == size(xTM,2)
        return new{N,T,V1,V2,M}(time, fp, xTM)
    end
end

for TMS in (:TMSol, :TMSol3)
    @eval begin
        @inline TM.expansion_point(a::$TMS) = getfield(a,:time)
        @inline TM.expansion_point(a::$TMS, n::Int) = getindex(getfield(a,:time),n)
        @inline flowpipe(a::$TMS) = getfield(a,:fp)
        @inline flowpipe(a::$TMS, n::Int) = getindex(getfield(a,:fp),n)
        @inline get_xTM(a::$TMS) = getfield(a,:xTM)
        @inline get_xTM(a::$TMS, n::Int) = getindex(getfield(a,:xTM),:,n)
        @inline TM.domain(a::$TMS) = domain.(getindex(getfield(a, :xTM), 1, :)) # vector!
        @inline TM.domain(a::$TMS, n::Int) = domain(getindex(getfield(a, :xTM), 1, n))
        @inline TS.get_numvars(a::$TMS) = get_numvars(a.xTM[1][0])

        # TMSol utilities
        @inline Base.firstindex(a::$TMS) = firstindex(a.time)
        @inline Base.lastindex(a::$TMS)  = lastindex(a.time)
        @inline Base.length(a::$TMS) = length(a.time)
        @inline Base.iterate(a::$TMS, state=0) = state ≥ lastindex(a) ? nothing : (a[state+1], state+1)
        @inline Base.eachindex(a::$TMS) = firstindex(a):lastindex(a)

        Base.getindex(a::$TMS, n::Integer) = getindex(get_xTM(a),:,n)
        Base.getindex(a::$TMS, u::UnitRange) = getindex(get_xTM(a),:,u)
        Base.getindex(a::$TMS, c::Colon) = getindex(get_xTM(a),:,c)
        Base.getindex(a::$TMS, n::Integer, m::Integer) = getindex(get_xTM(a),m,n)
        Base.getindex(a::$TMS, u::UnitRange, m::Integer) = getindex(get_xTM(a),m,u)
        Base.getindex(a::$TMS, c::Colon, m::Integer) = getindex(get_xTM(a),m,c)
    end
end
