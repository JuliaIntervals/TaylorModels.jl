# Methods ported from TaylorIntegration, specialized for validated integrations

# VectorCacheVI
# @unpack tv, xv, xaux, t, x, dx, rv, xauxI, tI, xI, dxI, rvI,
#         xTMN, xTM1v, rem, parse_eqs = cacheVI
struct VectorCacheVI{
        TV,XV,
        XAUX,TT,X,DX,RV,
        # TVI,
        XAUXI,TTI,XI,DXI,RVI,
        XTMN,XTM1V,REM,PARSE_EQS} <: TI.AbstractVectorCache
    tv::TV
    xv::XV
    xaux::XAUX
    t::TT
    x::X
    dx::DX
    rv::RV
    # tvI::TVI
    # psolI::PSOLI
    xauxI::XAUXI
    tI::TTI
    xI::XI
    dxI::DXI
    rvI::RVI
    xTMN::XTMN
    xTM1v::XTM1V
    rem::REM
    parse_eqs::PARSE_EQS
end


# init_cache_VI
"""
    init_cache_VI(t0::T, x0::Array{Interval{U},1},
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true)
    init_cache_VI(t0::T, x0::Array{TaylorModel1{TaylorN{U},U},1},
        maxsteps::Int, orderT::Int, ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true)

Initialize the internal integration variables and normalize the given interval
box to the domain `[-1, 1]^n`. If `x0` corresponds to a vector of TaylorModel1,
it is assumed that the domain of the TaylorN variables is normalized to the domain
`[-1, 1]^n`.

"""
function init_cache_VI(t0::T, x0::Array{Interval{U},1},
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true) where {U,T,F}

    dof = length(x0)
    zI = zero(Interval{U})
    S  = symmetric_box(dof, U)
    zB = zero(S)
    qaux = mid.(x0)
    δq0 = x0 .- qaux
    # Similar to `normalize_taylor`, without returning a TaylorN{Interval{T}}, but TaylorN{T}
    q0 = Array{TaylorN{typeof(qaux[1])}}(undef, dof)
    @inbounds for ind in eachindex(q0)
        q0[ind] = qaux[ind] + TaylorN(ind, order=orderQ) * radius(δq0[ind])
    end

    # Initialize the vector of Taylor1{TaylorN{U}} expansions
    t, x, dx = TI.init_expansions(t0, q0, orderT)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsX, rv = TI._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Initialize variables for Taylor integration with intervals
    tI, xI, dxI = TI.init_expansions(t0, x0, orderT+1)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsI, rvI = TI._determine_parsing!(parse_eqs, f!, tI, xI, dxI, params)

    if parse_eqs && parse_eqsX == parse_eqsI == true
        t, x, dx = TI.init_expansions(t0, q0, orderT)
        tI, xI, dxI = TI.init_expansions(t0, x0, orderT+1)
    end

    # More initializations
    xTMN  = Array{TaylorModelN{dof,T,T}}(undef, dof)
    @. xTMN = TaylorModelN(constant_term(x), zI, (zB,), (S,))
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    rem   = Array{Interval{T}}(undef, dof)
    for ind in eachindex(x)
        rem[ind] = zI
        xTM1v[ind, 1] = TaylorModel1(deepcopy(x[ind]), zI, zI, zI)
    end

    # Initialize cache
    # @unpack tv, xv, xaux, t, x, dx, rv, xauxI, tI, xI, dxI, rvI,
    #         xTMN, xTM1v, rem, parse_eqs = cacheVI
    cacheVI = VectorCacheVI(
            Array{T}(undef, maxsteps + 1),
            Vector{Vector{Interval{U}}}(undef, maxsteps + 1),
            Array{Taylor1{TaylorN{U}}}(undef, dof),
            t, x, dx, rv,
            Array{Taylor1{Interval{U}}}(undef, dof),
            tI, xI, dxI, rvI,
            xTMN, xTM1v, rem,
            parse_eqsX)

    return cacheVI
end
function init_cache_VI(t0::T, xTM::Array{TaylorModel1{TaylorN{U},U},1},
        maxsteps::Int, orderT::Int, ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true) where {U,T,F}

    dof = length(xTM)
    zI = zero(Interval{U})
    S  = symmetric_box(dof, U)
    zB = zero(S)

    # Initialize the vector of Taylor1{TaylorN} expansions
    q0 = constant_term.(polynomial.(xTM))
    x0 = evaluate(evaluate.(polynomial.(xTM)), Vector(S))

    # Initialize the vector of Taylor1{TaylorN{U}} expansions
    t, x, dx = TI.init_expansions(t0, q0, orderT)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsX, rv = TI._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Initialize variables for Taylor integration with intervals
    tI, xI, dxI = TI.init_expansions(t0, x0, orderT+1)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsI, rvI = TI._determine_parsing!(parse_eqs, f!, tI, xI, dxI, params)

    if parse_eqs && parse_eqsX == parse_eqsI == true
        t, x, dx = TI.init_expansions(t0, q0, orderT)
        tI, xI, dxI = TI.init_expansions(t0, x0, orderT+1)
    end

    # More initializations
    xTMN  = Array{TaylorModelN{dof,T,T}}(undef, dof)
    @. xTMN = TaylorModelN(constant_term(x), zI, (zB,), (S,))
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    rem   = Array{Interval{T}}(undef, dof)
    @. begin
        rem = zI
        xTM1v[:, 1] = TaylorModel1(deepcopy(x), zI, zI, zI)
    end

    # Initialize cache
    # @unpack tv, xv, xaux, t, x, dx, rv, xauxI, tI, xI, dxI, rvI,
    #         xTMN, xTM1v, rem, parse_eqs = cacheVI
    cacheVI = VectorCacheVI(
            Array{T}(undef, maxsteps + 1),
            Vector{Vector{Interval{U}}}(undef, maxsteps + 1),
            Array{Taylor1{TaylorN{U}}}(undef, dof),
            t, x, dx, rv,
            Array{Taylor1{Interval{U}}}(undef, dof),
            tI, xI, dxI, rvI,
            xTMN, xTM1v, rem,
            parse_eqsX)

    return cacheVI
end


# stepsize
function TI.stepsize(x::Taylor1{Interval{U}}, epsilon::T) where {T<:Real,U<:Number}
    R = promote_type(typeof(norm(constant_term(x), Inf)), T)
    ord = x.order
    h = typemax(R)
    z = zero(R)
    for k in (ord - 1, ord)
        @inbounds aux = norm(x[k], Inf)
        isequal_interval(aux, z) && continue
        aux1 = TI._stepsize(aux, epsilon, k)
        h = min(h, aux1)
    end
    return h::R
end

function TI.stepsize(q::AbstractArray{Taylor1{Interval{U}},N}, epsilon::T) where
        {T<:Real,U<:IA.NumTypes,N}
    R = promote_type(typeof(norm(constant_term(q[1]), Inf)), T)
    h = typemax(R)
    for i in eachindex(q)
        @inbounds hi = TI.stepsize(q[i], epsilon)
        h = min(h, hi)
    end

    # If `isinf(h)==true`, we use the maximum (finite)
    # step-size obtained from all coefficients as above.
    # Note that the time step is independent from `epsilon`.
    if isequal_interval(h, typemax(R))
        h = zero(R)
        for i in eachindex(q)
            @inbounds hi = TI._second_stepsize(q[i], epsilon)
            h = max(h, hi)
        end
    end
    return h::R
end
