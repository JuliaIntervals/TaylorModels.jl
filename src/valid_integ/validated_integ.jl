# Some methods for validated integration of ODEs

const _DEF_MINABSTOL = 1.0e-50


function validated_integ(f!, X0, t0::T, tmax::T, orderQ::Int, orderT::Int,
        abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        check_property::F=(t, x)->true) where {T<:Real, F}

    # Initialize the Taylor1 expansions
    # N = get_numvars()
    # dof = N
    dof = length(X0)
    t  = t0 + Taylor1(T, orderT )
    tI = t0 + Taylor1(T, orderT+1)

    # Internals: jet transport integration
    x     = Array{Taylor1{TaylorN{T}}}(undef, dof)
    dx    = Array{Taylor1{TaylorN{T}}}(undef, dof)

    # Internals: Taylor1{Interval{T}} integration
    xI    = Array{Taylor1{Interval{T}}}(undef, dof)
    dxI   = Array{Taylor1{Interval{T}}}(undef, dof)

    # Aux vars (parse_eqs == false)
    xaux  = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xauxI = Array{Taylor1{Interval{T}}}(undef, dof)

    # Set initial conditions
    initialize!(X0, orderQ, orderT, x, xI)
    f!(dx, x, params, t)
    f!(dxI, xI, params, tI)

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
    parse_eqsI, rvI = TI._determine_parsing!( parse_eqs, f!, tI, xI, dxI, params)
    parse_eqs,  rv  = TI._determine_parsing!( parse_eqs, f!, t, x, dx, params)
    @assert parse_eqs == parse_eqsI

    # Re-initialize the Taylor1 expansions
    t  = t0 + Taylor1(T, orderT )
    tI = t0 + Taylor1(T, orderT+1)
    initialize!(X0, orderQ, orderT, x, xI)
    return _validated_integ!(f!, t0, tmax, x, dx, xI, dxI,
        xaux, xauxI, rv, rvI, orderT, abstol, params,
        parse_eqs, maxsteps, adaptive, minabstol, absorb, check_property)
end



"""
    initialize!(X0::AbstractVector{<:Interval}, orderQ, orderT, x, xI)
    initialize!(X0::AbstractVector{<:Interval}, orderQ, orderT, x)

Initialize the internal integration variables and normalize the given interval
box to the domain `[-1, 1]^n`.
"""
function initialize!(X0::SVector{N,Interval{T}}, orderQ, orderT, x, xI) where {N,T}
    @assert N == get_numvars()

    # center of the box and vector of widths
    q0 = mid.(X0)
    δq0 = X0 .- q0

    qaux = normalize_taylor.(q0 .+ TaylorN.(1:N, order=orderQ), (δq0,), true)
    @. begin
        x = Taylor1(qaux, orderT)
        # dx = Taylor1(zero(qaux), orderT)
        xI = Taylor1(X0, orderT+1)
        # dxI = Taylor1(zero(X0), orderT+1)
    end

    return nothing
end
function initialize!(X0::SVector{N,Interval{T}}, orderQ, orderT, x) where {N,T}
    @assert N == get_numvars()
    q0 = mid.(X0)
    δq0 = X0 .- q0

    qaux = normalize_taylor.(q0 .+ TaylorN.(1:N, order=orderQ), (δq0,), true)
    @. begin
        x = Taylor1(qaux, orderT)
        # dx = x
    end
    return nothing
end

"""
    initialize!(X0::Vector{TaylorModel1{TaylorN{T}, T}}, orderQ, orderT, x, xI)
    initialize!(X0::Vector{TaylorModel1{TaylorN{T}, T}}, orderQ, orderT, x)

Initialize the auxiliary integration variables assuming that the given vector
of taylor models `X0` is normalized to the domain `[-1, 1]^n` in space.
"""
function initialize!(X0::Vector{TaylorModel1{TaylorN{T},T}}, orderQ, orderT, x, xI) where {T}
    # nomalized domain
    N = get_numvars()
    S = symmetric_box(N, T)

    qaux = constant_term.(polynomial.(X0))
    @. begin
        x = Taylor1(qaux, orderT)
        # dx = x
        # we assume that qaux is normalized to S=[-1, 1]^N
        xI = Taylor1(evaluate(qaux, (S,)), orderT+1)
        # dxI = xI
    end
    return nothing
end
function initialize!(X0::Vector{TaylorModel1{TaylorN{T},T}}, orderQ, orderT, x) where {T}
    # nomalized domain
    N = get_numvars()

    qaux = constant_term.(polynomial.(X0))
    @. begin
        x = Taylor1(qaux, orderT)
        # dx = x
    end
    return nothing
end



function _validated_integ!(f!, t0::T, tmax::T, x, dx, xI, dxI,
                xaux, xauxI, rv, rvI, orderT::Int, abstol::T, params,
                parse_eqs::Bool, maxsteps::Int, adaptive::Bool, minabstol::T,
                absorb::Bool, check_property::F) where {T<:Real,F}

    # Set proper parameters for jet transport
    N = get_numvars()
    dof = N

    # Some variables
    zt = zero(t0)
    zI = zero(Interval{T})
    zB = fill(zero(Interval{T}), SVector{N})
    S  = symmetric_box(N, T)
    t  = t0 + Taylor1(T, orderT)
    tI = t0 + Taylor1(T, orderT+1)
    xTMN  = Array{TaylorModelN{N,T,T}}(undef, dof)
    @. xTMN = TaylorModelN(constant_term(x), zI, (zB,), (S,))

    # Allocation of vectors
    # Output
    tv    = Array{T}(undef, maxsteps+1)
    xv    = Array{SVector{N,Interval{T}}}(undef, maxsteps+1)
    @inbounds tv[1] = t0
    @inbounds xv[1] = evaluate(xTMN, S)
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    rem   = Array{Interval{T}}(undef, dof)
    @. begin
        rem = zI
        xTM1v[:, 1] = TaylorModel1(deepcopy(x), zI, zI, zI)
    end

    # Direction of the integration
    sign_tstep = copysign(1, tmax-t0)

    local _success # if true, the validation step succeeded
    red_abstol = abstol

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax

        # Validated step of the integration
        (_success, δt, red_abstol) = validated_step!(Val(parse_eqs), f!,
                    t, x, dx, tI, xI, dxI, xaux, xauxI, rv, rvI,
                    t0, tmax, sign_tstep, xTMN, rem, zB, S,
                    orderT, red_abstol, params,
                    adaptive, minabstol, absorb, check_property)
        δtI = sign_tstep * interval(zt, sign_tstep*δt)

        # New initial conditions and time, and output vectors
        nsteps += 1
        @inbounds tv[nsteps] = t0
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tI[0] = t0
        @. begin
            xTM1v[:, nsteps] = TaylorModel1(deepcopy(x), rem, zI, δtI) # deepcopy is needed!
            x = Taylor1(evaluate(x, δt), orderT)
            # dx = Taylor1(zero(constant_term(x)), orderT)
            xI = Taylor1(evaluate(xTMN, (S,)), orderT+1)
            # dxI = xI
        end
        xv[nsteps] = evaluate(xTMN, S) # interval box

        # Try to increase `red_abstol` if `adaptive` is true
        if adaptive
            red_abstol = min(abstol, 10*red_abstol)
        end

        # If the integration step is unsuccessfull, break with a warning; note that the
        # last integration step (which was not successfull) is returned
        if !_success
            @warn("""
            Exiting due to unsuccessfull step
            """, _success, t0)
            break
        end

        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end

    end

    return TMSol(view(tv,1:nsteps), view(xv,1:nsteps), view(xTM1v,:,1:nsteps))
end



"""
    validated-step!
"""
function validated_step!(vB::Val{B}, f!,
        t::Taylor1{T}, x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        tI::Taylor1{T}, xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}},
        xaux::Vector{Taylor1{TaylorN{T}}}, xauxI::Vector{Taylor1{Interval{T}}},
        rv::TI.RetAlloc{Taylor1{TaylorN{T}}}, rvI::TI.RetAlloc{Taylor1{Interval{T}}},
        t0::T, tmax::T, sign_tstep::Int,
        xTMN::Vector{TaylorModelN{N,T,T}}, rem::Vector{Interval{T}},
        zbox::SVector{N,Interval{T}}, symIbox::SVector{N,Interval{T}},
        orderT::Int, abstol::T, params,
        adaptive::Bool, minabstol::T, absorb::Bool,
        check_property::F) where {B,N,T,F}

    # One step integration (non-validated)
    δt = TI.taylorstep!(vB, f!, t, x, dx, xaux, abstol, params, rv)
    f!(dx, x, params, t)  # Update last t coeff `dx[:][orderT]`

    # One step integration for the *initial box*
    TI.taylorstep!(vB, f!, tI, xI, dxI, xauxI, abstol, params, rvI)
    f!(dxI, xI, params, tI)  # Update last t coeff `dxI[:][orderT+1]`

    # Step size
    δt = min(δt, sign_tstep*(tmax-t0))
    δt = sign_tstep * δt

    return _validation(f!, t, x, dx, xI, dxI, δt, sign_tstep,
        xTMN, rem, zbox, symIbox, orderT, abstol, params,
        adaptive, minabstol, absorb, check_property)
end

function _validation(f!, t::Taylor1{T}, x::Vector{Taylor1{TaylorN{T}}},
        dx::Vector{Taylor1{TaylorN{T}}}, xI::Vector{Taylor1{Interval{T}}},
        dxI::Vector{Taylor1{Interval{T}}},
        δt, sign_tstep::Int,
        xTMN::Vector{TaylorModelN{N,T,T}}, rem::Vector{Interval{T}},
        zbox::SVector{N,Interval{T}}, symIbox::SVector{N,Interval{T}},
        orderT::Int, abstol::T, params,
        adaptive::Bool, minabstol::T, absorb::Bool,
        check_property::Function=(t, x)->true) where {N,T}

    # Test if `check_property` is satisfied; if not, half the integration time.
    # If after 25 checks `check_property` is not satisfied, throw an error.
    issatisfied = false
    rem_old = copy(rem)

    local _success
    reduced_abstol = abstol
    bool_red = true
    while bool_red
        # Validate the solution: remainder consistent with Schauder thm
        δtI = sign_tstep * interval(0, sign_tstep*δt)
        (_success, Δ) = remainder_taylorstep!(f!, t, x, dx, xI, dxI, symIbox, δtI, params)
        # (_success, Δ, δtI) = remainder_taylorstep2!(f!, t, x, dx, xI, dxI, symIbox, δtI, params)
        # δt = sup(δtI)

        # Shrink stepsize δt if adaptive is true and _success is false
        if !_success
            if adaptive
                bool_red = reduced_abstol > minabstol
                if bool_red
                    reduced_abstol = reduced_abstol/10
                    δt = δt * 0.1^(1/orderT)
                    continue
                else
                    @warn("Minimum absolute tolerance reached: ", t[0], δt, reduced_abstol, Δ)
                end
            else
                @warn("It cannot prove existence and unicity of the solution: ", t[0], δt, Δ, _success)
            end
        end

        # Remainder
        rem .= rem_old .+ Δ

        # Create TaylorModelN to store remainders and evaluation
        @inbounds begin
            for i in eachindex(x)
                xTMN[i] = fp_rpa( TaylorModelN(x[i](δtI), rem[i], zbox, symIbox) )

                # If remainder is still too big, do it again
                j = 0
                while absorb && (j < 10) && (mag(rem[i]) > 1.0e-10)
                    t[0] == 0 && println("absorb_remainder ")
                    j += 1
                    xTMN[i] = absorb_remainder(xTMN[i])
                    rem[i] = remainder(xTMN[i])
                end
            end
            xvv = evaluate(xTMN, symIbox) # interval box

            issatisfied = check_property(t[0]+δt, xvv)
            if !issatisfied
                # δt = δt/2
                bool_red = reduced_abstol > minabstol
                @info("issatisfied: ", bool_red, δt)
                if bool_red
                    reduced_abstol = reduced_abstol / 10
                    δt = δt * 0.1^(1/orderT)
                    continue
                else
                    @warn("Minimum absolute tolerance reached (issatisfied): ", t[0], δt, reduced_abstol, Δ)
                end
            end
        end # @inbounds
        break
    end

    if !issatisfied && !adaptive
        @warn("""
            Something went wrong:
            """, issatisfied, adaptive, bool_red, reduced_abstol,
            t[0], δt, Δ, xvv, check_property(t0+δt, xvv)
        )
    end

    return (_success, δt, reduced_abstol)
end


"""
    remainder_taylorstep!(f!, t, x, dx, xI, dxI, δI, δt, params)

Returns a remainder for the integration step for the dependent variables (`x`)
checking that the solution satisfies the criteria for existence and uniqueness.
"""
function remainder_taylorstep!(f!::Function, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}},
        δI::SVector{N,Interval{T}}, δtI::Interval{T}, params) where {N,T}

    orderT = get_order(dx[1])
    aux = δtI^(orderT+1)
    Δx  = SVector{length(xI)}([xI[i][orderT+1] for i in eachindex(xI)]) * aux
    Δdx = SVector{length(xI)}([dxI[i][orderT+1] for i in eachindex(xI)]) * aux
    Δ0  = SVector{length(x)}([dx[i][orderT](δI) for i in eachindex(x)]) * aux / (orderT+1)
    Δ = Δ0 + Δdx * δtI
    Δxold = Δx

    # Checking existence and uniqueness
    iscontractive(Δ, Δx) && return (true, Δx)

    # If the check didn't work, compute new remainders. A new Δx is proposed,
    # and the corresponding Δdx is computed
    xxI  = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
    dxxI = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
    vv = Array{Interval{T}}(undef, N)
    for _ = 1:50
        # Remainder of Picard iteration
        Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δtI, Δx, Δ0, params)

        # Checking existence and uniqueness
        iscontractive(Δ, Δx) && return (true, Δx)
        # iscontractive(Δ, Δx) && return (true, _contract_iteration!(f!, t, x, dx, xxI, dxxI, δI, δtI, Δx, Δdx, Δ0, params))

        # Expand Δx in the directions needed
        Δxold = Δx
        if issubset_interval(Δ, Δx)
            @inbounds for ind in 1:N
                # Widen the directions where ⊂ does not hold
                vv[ind] = Δx[ind]
                if Δ[ind] == Δx[ind]
                    vv[ind] = widen.(Δ[ind])
                end
            end
            Δx = vv
            continue
        end
        Δx = Δ
    end

    return (false, Δx)
end

# function remainder_taylorstep2!(f!::Function, t::Taylor1{T},
#         x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
#         xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}},
#         δI::IntervalBox{N,T}, δtI::Interval{T}, params) where {N,T}

#     orderT = get_order(dx[1])
#     aux = δtI^(orderT+1)
#     # Estimate of remainder for x(t) based on interval integration xI
#     ΔxC = IntervalBox([xI[i][orderT+1] for i in eachindex(xI)])
#     Δx  = ΔxC * aux
#     # Estimate of remainder for \dot{x}(t) based on interval integration dxI
#     ΔdxC = IntervalBox([dxI[i][orderT+1] for i in eachindex(xI)])
#     Δdx  = ΔdxC * aux
#     # Integrated last term of \dot{x}(t) based on interval integration dxI
#     Δ0C = IntervalBox([dx[i][orderT](δI) for i in eachindex(x)])
#     Δ0  = Δ0C * aux / (orderT+1)
#     Δ = Δ0 + Δdx * δtI

#     # @show(ΔxC, ΔdxC, Δ0C)
#     # Checking existence and uniqueness
#     iscontractive(Δ, Δx) && return (true, Δx, δtI)

#     xxI  = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
#     dxxI = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)

#     δt = δtI
#     # @show(t[0], δtI)
#     for i = 1:10
#         δt /= 2
#         aux = δt^(orderT+1)
#         ΔΔx  = ΔxC * aux
#         ΔΔdx  = ΔdxC * aux
#         ΔΔ0  = Δ0C * aux / (orderT+1)
#         ΔΔ = ΔΔ0 + ΔΔdx * δt
#         # iscontractive(ΔΔ, ΔΔx) && (@show(δt, ΔΔx, Δx); return (true, ΔΔx, δt))
#         iscontractive(ΔΔ, ΔΔx) && return (true, _contract_iteration!(f!, t, x, dx, xxI, dxxI, δI, δt, ΔΔx, ΔΔdx, ΔΔ0, params), δt)
#         # if iscontractive(ΔΔ, ΔΔx)
#         #     ΔΔ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, ΔΔx, ΔΔ0, params)
#         #     @show(iscontractive(ΔΔ, ΔΔx))
#         #     break
#         # end
#     end

#     # # Checking existence and uniqueness
#     # iscontractive(Δ, Δx) && return (true, Δx)

#     # If the check didn't work, compute new remainder Δ using Δx, Δ0 and
#     # the evaluation of f! over intervals.
#     # A new Δx is proposed, and the corresponding Δdx is computed
#     # xxI  = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
#     # dxxI = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
#     vv = Array{Interval{T}}(undef, N)
#     Δxold = Δx
#     for _ = 1:50
#         Δx = Δ
#         # Remainder of Picard iteration
#         Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δtI, Δx, Δ0, params)

#         # Checking existence and uniqueness
#         iscontractive(Δ, Δx) && return (true, Δx, δtI)
#         # iscontractive(Δ, Δx) && return _contract_iteration!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δdx, Δ0, params)

#         # Expand Δx in the directions needed
#         Δxold = Δx
#         if Δ ⊆ Δx
#             @inbounds for ind in 1:N
#                 # Widen the directions where ⊂ does not hold
#                 vv[ind] = Δx[ind]
#                 if Δx[ind] ⊆ Δ[ind]
#                     vv[ind] = widen.(Δ[ind])
#                 end
#             end
#             Δ = IntervalBox(vv)
#             continue
#         end
#         # Δx = Δ
#     end

#     return (false, Δx, δtI)
# end



"""
    iscontractive(Δ, Δx)

Checks if `Δ .⊂ Δx` is satisfied. If ``Δ ⊆ Δx` is satisfied, it returns
`true` if all cases where `==` holds corresponds to the zero `Interval`.
"""
function iscontractive(Δ::Interval{T}, Δx::Interval{T}) where{T}
    (Δ ⊂ Δx || Δ == Δx == zero(Δ)) && return true
    return false
end
iscontractive(Δ::SVector{N,Interval{T}}, Δx::SVector{N,Interval{T}}) where{N,T} =
    all(iscontractive.(Δ[:], Δx[:]))


"""
    picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0, params)

Return the remainder of Picard operator
"""
function picard_remainder!(f!::Function, t::Taylor1{T},
    x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
    xxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    δI::SVector{N,Interval{T}}, δt::Interval{T},
    Δx::SVector{N,Interval{T}}, Δ0::SVector{N,Interval{T}}, params) where {N,T}

    # Extend `x` and `dx` to have interval coefficients
    zI = zero(δt)
    @. begin
        xxI = x + Δx
        # dxxI = dx + zI
    end

    # Compute `dxxI` from the equations of motion
    f!(dxxI, xxI, params, t)

    # Picard iteration, considering only the bound of `f` and the last coeff of f
    Δdx = Interval.(evaluate.( (dxxI - dx)(δt), Ref(δI) ))
    Δ = Δ0 + Δdx * δt
    return Δ
end


# # Picard iterations to contract further Δx, once Δ ⊂ Δx holds
# # **Currently not used**
# function _contract_iteration!(f!::Function, t::Taylor1{T},
#         x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
#         xxI::Vector{Taylor1{TaylorN{Interval{T}}}}, dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
#         δI::IntervalBox{N,T}, δt::Interval{T},
#         Δx::IntervalBox{N,T}, Δdx::IntervalBox{N,T}, Δ0::IntervalBox{N,T}, params) where {N,T}

#     # Some abbreviations
#     Δ = Δ0 + Δdx * δt
#     Δxold = Δx

#     # Picard contractions
#     for its = 1:10
#         # Remainder of Picard iteration
#         Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0, params)

#         # If contraction doesn't hold, return old bound
#         iscontractive(Δ, Δx) || return Δxold

#         # Contract estimate
#         Δxold = Δx
#         Δx = Δ
#     end

#     return Δxold
# end


"""
    absorb_remainder(a::TaylorModelN)

Returns a TaylorModelN, equivalent to `a`, such that the remainder
is mostly absorbed in the constant and linear coefficients. The linear shift assumes
that `a` is normalized to the interval box `(-1..1)^N`.

Ref: Xin Chen, Erika Abraham, and Sriram Sankaranarayanan,
"Taylor Model Flowpipe Construction for Non-linear Hybrid
Systems", in Real Time Systems Symposium (RTSS), pp. 183-192 (2012),
IEEE Press.
"""
function absorb_remainder(a::TaylorModelN{N,T,T}) where {N,T}
    Δ = remainder(a)
    orderQ = get_order(a)
    δ = symmetric_box(N, T)
    aux = diam(Δ)/(2N)
    rem = zero(Δ)

    # Linear shift
    lin_shift = mid(Δ) + sum((aux*TaylorN(i, order=orderQ) for i in 1:N))
    bpol = a.pol + lin_shift

    # Compute the new remainder
    aI = a(δ)
    bI = bpol(δ)

    if issubset_interval(bI, aI)
        rem = interval(inf(aI)-inf(bI), sup(aI)-sup(bI))
    elseif issubset_interval(aI, bI)
        rem = interval(inf(bI)-inf(aI), sup(bI)-sup(aI))
    else
        r_lo = inf(aI)-inf(bI)
        r_hi = sup(aI)-sup(bI)
        if r_lo > 0
            rem = interval(-r_lo, r_hi)
        else
            rem = interval( r_lo, -r_hi)
        end
    end

    return TaylorModelN(bpol, rem, expansion_point(a), domain(a))
end


# Postverify and define Taylor models to be returned
function scalepostverify_sw!(xTMN::Vector{TaylorModelN{N,T,T}},
        X::Vector{TaylorN{T}}) where {N,T}
    postverify = true
    x0 = expansion_point(xTMN[1])
    B = domain(xTMN[1])
    zI = zero(Interval{T})
    @inbounds for i in eachindex(xTMN)
        pol = polynomial(xTMN[i])
        ppol = fp_rpa(TaylorModelN(pol(X), zI, x0, B ))
        postverify = postverify && (xTMN[i](B) ⊆ ppol(B))
        xTMN[i] = copy(ppol)
    end
    @assert postverify """
        Failed to post-verify shrink-wrapping:
        X = $(linear_polynomial(X))
        xTMN = $(xTMN)
        """
    return postverify
end


"""
    shrink_wrapping!(xTMN::TaylorModelN)

Returns a modified inplace `xTMN`, which has absorbed the remainder
by the modified shrink-wrapping method of Florian Bünger.
The domain of `xTMN` is the normalized interval box `[-1,1]^N`.

Ref: Florian B\"unger, Shrink wrapping for Taylor models revisited,
Numer Algor 78:1001–1017 (2018), https://doi.org/10.1007/s11075-017-0410-1
"""
function shrink_wrapping!(xTMN::Vector{TaylorModelN{N,T,T}}) where {N,T}
    # Original domain of TaylorModelN should be the symmetric normalized box
    B = symmetric_box(N, T)
    @assert all(domain.(xTMN) .== (B,))
    zI = zero(Interval{T})
    x0 = zero(B)
    @assert all(expansion_point.(xTMN) .== (x0,))

    # Vector of independent TaylorN variables
    order = get_order(xTMN[1])
    X = [TaylorN(T, i, order=order) for i in 1:N]

    # Remainder of original TaylorModelN and componentwise mag
    rem = remainder.(xTMN)
    r = mag.(rem)
    qB = r .* B
    one_r = ones(eltype(r), N)

    # Shift to remove constant term
    xTN0 = constant_term.(xTMN)
    xTNcent = polynomial.(xTMN) .- xTN0
    xTNcent_lin = linear_polynomial(xTNcent)

    # Step 4 of Bünger algorithm: Jacobian (at zero) and its inverse
    jac = TaylorSeries.jacobian(xTNcent_lin)
    # If the conditional number is too large (inverse of jac is ill defined),
    # don't change xTMN
    cond(jac) > 1.0e4 && return one_r
    # Inverse of the Jacobian
    invjac = inv(jac)

    # Componentwise bound
    r̃ = mag.(invjac * qB) # qB <-- r .* B
    qB´ = r̃ .* B
    @assert issubset_interval(invjac * qB, qB´)

    # Step 6 of Bünger algorithm: compute g
    g = invjac*xTNcent .- X
    # g = invjac*(xTNcent .- xTNcent_lin)
    # ... and its jacobian matrix (full dependence!)
    jacmatrix_g = TaylorSeries.jacobian(g, X)

    # Alternative to Step 7: Check the validity of Eq 16 (or 17) for Lemma 2
    # of Bünger's paper, for s=0, and s very small. If it satisfies it,
    # postverify and return. Otherwise, use Bünger's step 7.
    q = 1.0 .+ r̃
    s = zero(q)
    @. q = 1.0 + r̃ + s
    jaq_q1 = jacmatrix_g * (q .- 1.0)
    eq16 = all(mag.(evaluate.(jaq_q1, Ref(q .* B))) .≤ s)
    if eq16
        postverify = scalepostverify_sw!(xTMN, q .* X)
        postverify && return q
    end
    s .= eps.(q)
    @. q = 1.0 + r̃ + s
    jaq_q1 = jacmatrix_g * (q .- 1.0)
    eq16 = all(mag.(evaluate.(jaq_q1, Ref(q .* B))) .≤ s)
    if eq16
        postverify = scalepostverify_sw!(xTMN, q .* X)
        postverify && return q
    end

    # Step 7 of Bünger algorithm: estimate of `q`
    # Some constants/parameters
    q_tol = 1.0e-12
    q = 1.0 .+ r̃
    ff = 65/64
    q_max = ff .* q
    s = zero(q)
    q_old = similar(q)
    q_1 = similar(q)
    jaq_q1 .= jacmatrix_g * r̃
    iter_max = 100
    improve = true
    iter = 0
    while improve && iter < iter_max
        qB .= q .* B
        q_1 .= q .- 1.0
        q_old .= q
        mul!(jaq_q1, jacmatrix_g, q_1)
        eq16 = all(evaluate.(jaq_q1, Ref(qB)) .≤ s)
        eq16 && break
        @inbounds for i in eachindex(xTMN)
            s[i] = mag( jaq_q1[i](qB) )
            q[i] = 1.0 + r̃[i] + s[i]
            # If q is too large, return xTMN unchanged
            q[i] > q_max[i] && return -one_r
        end
        improve = any( ((q .- q_old)./q) .> q_tol )
        iter += 1
    end
    # (improve || q == one_r) && return one_r
    # Compute final q and rescale X
    @. q = 1.0 + r̃ + ff * s
    @. X = q * X

    # Postverify
    postverify = scalepostverify_sw!(xTMN, X)

    return q
end
