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
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F,
        localsp::JetSpace, params = nothing;
        parse_eqs::Bool = true) where {U,T,F}

    dof = length(x0)
    zI = zero(Interval{U})
    S  = symmetric_box(dof, U)
    zB = zero(S)
    # Initialize Taylor1{TaylorN} expansions explicitly
    q0 = Array{TaylorN{U}}(undef, dof)
    # @inbounds for ind in eachindex(q0)
    #     q0[ind] = mid(x0[ind]) + TaylorN(localsp, ind, order=orderQ) * radius(x0[ind])
    # end
    q0 .= mid.(x0) .+ variables(localsp; order=orderQ) .* radius.(x0)
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
    xTMN .= TM.unsafe_TaylorModelN.(
        TaylorN.((localsp,), getcoeff.(xI[:], 0), orderQ), zI, (zB,), (S,))
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    rem   = Array{Interval{T}}(undef, dof)
    for ind in eachindex(x)
        rem[ind] = zI
        xTM1v[ind, 1] = TM.unsafe_TaylorModel1(deepcopy(x[ind]), zI, zI, zI)
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
        maxsteps::Int, orderT::Int, orderQ::Int, f!::F,
        localsp::JetSpace, params = nothing;
        parse_eqs::Bool = true) where {U,T,F}
    dof = length(xTM)
    S  = symmetric_box(dof, U)
    # Initialize the vector of Taylor1{TaylorN} expansions
    x0 = evaluate(evaluate.(polynomial.(xTM)), Vector(S))
    return init_cache_VI(t0, x0, maxsteps, orderT, orderQ,
                f!, localsp, params; parse_eqs = parse_eqs)
end

struct VectorCacheVI3{N,T,U} <: TI.AbstractVectorCache
    # Ouput stuff
    tv::Vector{T}
    xv::Vector{Vector{Interval{U}}}
    xTM1v::Matrix{TaylorModel1{TaylorModelN{N,T,U},U}}
    # Internals related to the integration
    xaux::Vector{Taylor1{TaylorN{T}}}
    t::Taylor1{T}
    x::Vector{Taylor1{TaylorN{T}}}
    dx::Vector{Taylor1{TaylorN{T}}}
    rv::TI.RetAlloc{Taylor1{TaylorN{T}}}
    t1N::TaylorModel1{TaylorModelN{N,T,U}}
    x1N::Vector{TaylorModel1{TaylorModelN{N,T,U},U}}
    dx1N::Vector{TaylorModel1{TaylorModelN{N,T,U},U}}
    x2N::Vector{TaylorModel1{TaylorModelN{N,T,U},U}}
    z1N::TaylorModel1{TaylorModelN{N,T,U}}
    vTN::Vector{TaylorN{T}}
    auxI::TaylorN{Interval{U}}
    x0New::Vector{Interval{U}}
    rem1::Vector{Interval{U}}
    rem2::Vector{Interval{U}}
    rem0::Vector{Interval{U}}
    # Preconditioning stuff
    vTMN::Vector{TaylorModelN{N,T,U}}
    leftTMN::Vector{TaylorModelN{N,T,U}}
    rightTMN::Vector{TaylorModelN{N,T,U}}
    remsQR::Vector{Interval{U}}
    linTN::Matrix{T}
    scaleV::Vector{T}
    # Param for the integration
    parse_eqs::Bool
end

init_cache_VI3(f!::F, t0::T, x0::Array{Interval{U},1},
        maxsteps::Int, orderT::Int, orderQ::Int,
        localsp::JetSpace, params = nothing;
        parse_eqs::Bool = true) where {U,T,F} =
    init_cache_VI3(f!, t0, SVector(x0...),
        maxsteps, orderT, orderQ, localsp, params;
        parse_eqs = parse_eqs)

function init_cache_VI3(f!::F, t0::T, x0::SVector{N,Interval{U}},
        maxsteps::Int, orderT::Int, orderQ::Int,
        localsp::JetSpace, params = nothing;
        parse_eqs::Bool = true) where {N,U,T,F}
    # Internal vars for assignements
    @assert N == length(x0) == get_numvars(localsp)
    zI = zero(Interval{U})
    symIbox = symmetric_box(N, U)
    zbox = zero(symIbox)
    vTN = Array{TaylorN{U}}(undef, N)
    # Initialize the vector of Taylor1{TaylorN{U}} expansions explicitly
    vTN .= mid.(x0) .+ variables(localsp; order=orderQ) .* radius.(x0)
    t, x, dx = TI.init_expansions(t0, vTN, orderT)
    # Determine if specialized jetcoeffs! method exists/works
    parse_eqsX, rv = TI._determine_parsing!(parse_eqs, f!, t, x, dx, params)
    if parse_eqsX
        TI.init_expansions!(x, dx, vTN, orderT)
    end

    # More initializations
    zN = TM.unsafe_TaylorModelN(zero(x[1][0]), zI, zbox, symIbox)
    uN = TM.unsafe_TaylorModelN( one(x[1][0]), zI, zbox, symIbox)
    t1N = TM.unsafe_TaylorModel1(Taylor1([zN, uN], orderT), zI, 0.0, zI)
    z1N = zero(t1N)
    # and initializations
    TT = TaylorModel1{TaylorModelN{N,T,U}, U}
    x1N = Array{TT}(undef, N)
    dx1N = Array{TT}(undef, N)
    x2N = Array{TT}(undef, N)
    xTM1v = Array{TT}(undef, N, maxsteps+1)
    vTMN = Vector{TaylorModelN{N,T,U}}(undef, N)
    x0New = Vector{Interval{T}}(undef, N)
    rem1 = Array{Interval{T}}(undef, N)
    rem2 = Array{Interval{T}}(undef, N)
    rem0 = Array{Interval{T}}(undef, N)
    leftTMN  = Vector{TaylorModelN{N,T,U}}(undef, N)
    rightTMN = Vector{TaylorModelN{N,T,U}}(undef, N)
    remsQR   = Array{Interval{U}}(undef, N)
    for i in eachindex(x1N)
        dx1N[i] = deepcopy(z1N)
        x2N[i]  = deepcopy(z1N)
        x0New[i] = zI
        rem1[i] = zI
        rem2[i] = zI
        rem0[i] = zI
        # x1N[i] = TaylorModel1(deepcopy(x[i]), zI, 0.0, zI)
        # xTM1v[i, 1] = deepcopy(x1N[i])
        x1N[i] = deepcopy(z1N)
        for k in 1:maxsteps+1
            xTM1v[i, k] = deepcopy(z1N)
        end
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
        leftTMN[i] = zero(vTMN[i])
        rightTMN[i] = zero(vTMN[i])
        remsQR[i] = zI
    end
    linTN = zeros(T, N, N)
    scaleV = zeros(T, N)
    auxI = TM._aux_horner(vTMN[1])

    # Initialize cache
    return VectorCacheVI3{N,T,U}(
            Array{T}(undef, maxsteps + 1), #tv
            Vector{Vector{Interval{U}}}(undef, maxsteps + 1), #xv
            xTM1v,
            Array{Taylor1{TaylorN{U}}}(undef, N), #xaux
            t, x, dx, rv,
            t1N, x1N, dx1N, x2N, z1N, vTN, auxI,
            x0New, rem1, rem2, rem0,
            vTMN, leftTMN, rightTMN, remsQR, linTN, scaleV,
            parse_eqsX)
end

function init_cache_VI3(f!::F, t0::T,
        xTM::Array{TaylorModel1{TaylorModelN{N,T,U},U},1},
        maxsteps::Int, orderT::Int, orderQ::Int,
        localsp::JetSpace, params = nothing;
        parse_eqs::Bool = true) where {N,U,T,F}

    # Initialize the vector of Taylor1{TaylorN} expansions
    x0 = evaluate(evaluate.(xTM, xTM[1].x0), centered_dom(xTM[1].pol[0]))
    return init_cache_VI3(f!, t0, SVector(x0...), maxsteps, orderT, orderQ,
        localsp, params; parse_eqs = parse_eqs)
end
