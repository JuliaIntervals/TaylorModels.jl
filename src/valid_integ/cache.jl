# Methods ported from TaylorIntegration, specialized for validated integrations

# VectorCacheVI
struct VectorCacheVI{
        TV,XV,
        XAUX,TT,X,DX,RV,
        XAUXI,TTI,XI,DXI,RVI,
        XTMN,XTM1V,REM,PARSE_EQS} <: TI.AbstractVectorCache
    tv::TV
    xv::XV
    xaux::XAUX
    t::TT
    x::X
    dx::DX
    rv::RV
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
    # Initialize Taylor1{TaylorN} expansions explicitly
    q0 = Array{TaylorN{U}}(undef, dof)
    @inbounds for ind in eachindex(q0)
        q0[ind] = mid(x0[ind]) + TaylorN(ind, order=orderQ) * radius(x0[ind])
    end

    # Initialize the vector of Taylor1{TaylorN{U}} expansions
    t, x, dx = TI.init_expansions(t0, q0, orderT)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsX, rv = TI._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Initialize variables for Taylor integration with intervals
    tI, xI, dxI = TI.init_expansions(t0, x0, orderT+1)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsI, rvI = TI._determine_parsing!(parse_eqs, f!, tI, xI, dxI, params)

    if parse_eqsX && parse_eqsI
        t, x, dx = TI.init_expansions(t0, q0, orderT)
        tI, xI, dxI = TI.init_expansions(t0, x0, orderT+1)
    end

    # More initializations
    xTMN  = Array{TaylorModelN{dof,Interval{T},T}}(undef, dof)
    @. xTMN = TaylorModelN(TaylorN(getcoeff(xI[:], 0), orderQ), zI, (zB,), (S,))
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    rem   = Array{Interval{T}}(undef, dof)
    for ind in eachindex(x)
        rem[ind] = zI
        xTM1v[ind, 1] = TaylorModel1(deepcopy(x[ind]), zI, zI, zI)
    end

    # Initialize cache
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
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true) where {U,T,F}
    dof = length(xTM)
    S  = symmetric_box(dof, U)
    # Initialize the vector of Taylor1{TaylorN} expansions
    x0 = evaluate(evaluate.(polynomial.(xTM)), Vector(S))
    return init_cache_VI(t0, x0, maxsteps, orderT, orderQ,
                f!, params; parse_eqs = parse_eqs)
end

struct VectorCacheVI3{
        TV,XV,
        XAUX,TT,X,RV,
        TM1TMN,VTM1TMN,VTN,VTMN,TN,
        XTM1V,REM,PARSE_EQS} <: TI.AbstractVectorCache
    tv::TV
    xv::XV
    xaux::XAUX
    t::TT
    x::X
    dx::X
    rv::RV
    t1N::TM1TMN
    x1N::VTM1TMN
    dx1N::VTM1TMN
    x2N::VTM1TMN
    z1N::TM1TMN
    vTN::VTN
    vTMN::VTMN
    auxN::TN
    xTM1v::XTM1V
    x0New::REM
    rem1::REM
    rem2::REM
    parse_eqs::PARSE_EQS
end

function init_cache_VI3(t0::T, x0::Array{Interval{U},1},
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true) where {U,T,F}

    dof = length(x0)
    zI = zero(Interval{U})
    symIbox = symmetric_box(dof, U)
    zbox = zero(symIbox)
    vTN = Array{TaylorN{U}}(undef, dof)

    # Initialize the vector of Taylor1{TaylorN{U}} expansions explicitly
    @inbounds for ind in eachindex(vTN)
        vTN[ind] = mid(x0[ind]) + TaylorN(ind, order=orderQ) * radius(x0[ind])
    end
    t, x, dx = TI.init_expansions(t0, vTN, orderT)
    auxN = zero(vTN[1])

    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsX, rv = TI._determine_parsing!(parse_eqs, f!, t, x, dx, params)
    if parse_eqsX
        TI.init_expansions!(x, dx, vTN, orderT)
    end

    # More initializations
    zN = TaylorModelN(zero(x[1][0]), zI, zbox, symIbox)
    uN = TaylorModelN( one(x[1][0]), zI, zbox, symIbox)
    t1N = TaylorModel1(Taylor1([zN, uN], orderT), zI, 0.0, zI)
    z1N = zero(t1N)

    TT = TaylorModel1{TaylorModelN{dof,T,U}, U}
    x1N = Array{TT}(undef, dof)
    dx1N = Array{TT}(undef, dof)
    x2N = Array{TT}(undef, dof)
    xTM1v = Array{TT}(undef, dof, maxsteps+1)
    vTMN = Vector{TaylorModelN{dof,T,U}}(undef, dof)
    x0New = Vector{Interval{T}}(undef, dof)
    rem1   = Array{Interval{T}}(undef, dof)
    rem2   = Array{Interval{T}}(undef, dof)
    for i in eachindex(x1N)
        dx1N[i] = deepcopy(z1N)
        x2N[i]  = deepcopy(z1N)
        x0New[i] = zI
        rem1[i] = zI
        rem2[i] = zI
        # x1N[i] = TaylorModel1(deepcopy(x[i]), zI, 0.0, zI)
        # xTM1v[i, 1] = deepcopy(x1N[i])
        x1N[i] = deepcopy(z1N)
        xTM1v[i, 1] = deepcopy(z1N)
        for it in eachindex(x[i].coeffs)
            for iq in eachindex(x[i].coeffs[it].coeffs)
                for hq in eachindex(x[i].coeffs[it].coeffs[iq].coeffs)
                    x1N[i].pol.coeffs[it].pol.coeffs[iq].coeffs[hq] =
                        x[i].coeffs[it].coeffs[iq].coeffs[hq]
                    xTM1v[i,1].pol.coeffs[it].pol.coeffs[iq].coeffs[hq] =
                        x[i].coeffs[it].coeffs[iq].coeffs[hq]
                end
            end
        end
        vTMN[i] = evaluate(polynomial(x1N[i]), 0.0)
    end

    # Initialize cache
    cacheVI = VectorCacheVI3(
            Array{T}(undef, maxsteps + 1), #tv
            Vector{Vector{Interval{U}}}(undef, maxsteps + 1), #xv
            Array{Taylor1{TaylorN{U}}}(undef, dof), #xaux
            t, x, dx, rv,
            t1N, x1N, dx1N, x2N, z1N, vTN, vTMN, auxN,
            xTM1v, x0New, rem1, rem2,
            parse_eqsX)

    return cacheVI
end
function init_cache_VI3(t0::T, xTM::Array{TaylorModel1{TaylorModelN{N,T,U},U},1},
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F, params = nothing;
        parse_eqs::Bool = true) where {N,U,T,F}
    dof = length(xTM)
    symIbox  = symmetric_box(dof, U)
    # Initialize the vector of Taylor1{TaylorN} expansions
    x0 = evaluate(evaluate.(polynomial.(xTM)), Vector(symIbox))
    return init_cache_VI3(t0, x0, maxsteps, orderT, orderQ, f!, params;
        parse_eqs = parse_eqs)
end
