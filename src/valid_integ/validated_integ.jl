# Some methods for validated integration of ODEs

const _DEF_MINABSTOL = 1.0e-50


function validated_integ(f!, X0::AbstractVector{Interval{U}},
        t0::T, tmax::T, orderQ::Int, orderT::Int,
        abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        check_property::F=(t, x)->true) where {T<:Real, U, F}

    # Initialize cache
    vX0 = Vector(X0)
    cacheVI = init_cache_VI(t0, vX0, maxsteps, orderT, orderQ, f!, params; parse_eqs)

    return _validated_integ!(f!, vX0, t0, tmax, abstol, cacheVI, params,
        maxsteps, adaptive, minabstol, absorb, check_property)
end

function validated_integ(f!, X0::Vector{TaylorModel1{TaylorN{T}, U}},
        t0::T, tmax::T, orderQ::Int, orderT::Int,
        abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        check_property::F=(t, x)->true) where {T<:Real, U, F}

    # Initialize cache
    cacheVI = init_cache_VI(t0, X0, maxsteps, orderT, orderQ, f!, params; parse_eqs)
    q0 = evaluate(constant_term.(polynomial.(X0)), Vector(symmetric_box(length(X0),U)))

    return _validated_integ!(f!, q0, t0, tmax, abstol, cacheVI, params,
        maxsteps, adaptive, minabstol, absorb, check_property)
end


function _validated_integ!(f!, q0, t0::T, tmax::T, abstol::T,
        cacheVI::VectorCacheVI, params,
        maxsteps::Int, adaptive::Bool, minabstol::T,
        absorb::Bool, check_property::F) where {T<:Real,F}

    # Unpack caches
    @unpack tv, xv, xaux, t, x, dx, rv, xauxI, tI, xI, dxI, rvI,
            xTMN, xTM1v, rem, parse_eqs = cacheVI

    # Initial conditions
    sign_tstep = copysign(1, tmax - t0)
    dof = length(q0)
    orderT = get_order(t)
    zt = zero(t0)
    zI = zero(Interval{T})
    symIbox = symmetric_box(dof, T)
    zbox = zero(symIbox)
    @inbounds xv[1] = q0
    @inbounds tv[1] = t0

    # Integration
    nsteps = 1
    local _success # if true, the validation step succeeded
    red_abstol = abstol
    VV = Val(parse_eqs)
    while sign_tstep*t0 < sign_tstep*tmax
        # Validated step of the integration
        (_success, δt, red_abstol) = validated_step!(VV, f!,
                    t, x, dx, tI, xI, dxI, xaux, xauxI, rv, rvI,
                    t0, tmax, sign_tstep,
                    xTMN, rem,
                    zbox, symIbox, orderT, red_abstol, params,
                    adaptive, minabstol, absorb, check_property)
        δtI = interval(zt, sign_tstep*δt)

        # Save different objects
        nsteps += 1
        for ind in eachindex(x)
            xTM1v[ind, nsteps] = TaylorModel1(deepcopy(x[ind]), rem[ind], zI, δtI) # deepcopy is needed!
            x[ind] = Taylor1(evaluate(x[ind], δt), orderT)
            # dx = Taylor1(zero(constant_term(x)), orderT)
            xI[ind] = Taylor1(evaluate(xTMN[ind], symIbox), orderT+1)
            # dxI = xI
        end

        # New initial conditions and time, and output vectors
        @inbounds tv[nsteps] = t0
        xv[nsteps] = copy.(constant_term.(xI[:]))       # Vector{Interval}
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tI[0] = t0

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

    return TMSol(view(tv,1:nsteps), view(xv,1:nsteps,), view(xTM1v,:,1:nsteps))
end



"""
    validated_step!
"""
function validated_step!(vB::Val{B}, f!,
        t::Taylor1{T}, x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        tI::Taylor1{T}, xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}},
        xaux::Vector{Taylor1{TaylorN{T}}}, xauxI::Vector{Taylor1{Interval{T}}},
        rv::TI.RetAlloc{Taylor1{TaylorN{T}}}, rvI::TI.RetAlloc{Taylor1{Interval{T}}},
        t0::T, tmax::T, sign_tstep::Int,
        xTMN::Vector{TaylorModelN{N,Interval{T},T}}, rem::Vector{Interval{T}},
        zbox::AbstractVector{Interval{T}}, symIbox::AbstractVector{Interval{T}},
        orderT::Int, abstol::T, params,
        adaptive::Bool, minabstol::T, absorb::Bool,
        check_property::F) where {N,B,T,F}

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
        xTMN::Vector{TaylorModelN{N,Interval{T},T}}, rem::Vector{Interval{T}},
        zbox::AbstractVector{Interval{T}}, symIbox::AbstractVector{Interval{T}},
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
                xTMN[i] = TaylorModelN(x[i](δtI), rem[i], SVector{N}(zbox), SVector{N}(symIbox))

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
        δI::AbstractVector{Interval{T}}, δtI::Interval{T}, params) where {T}

    orderT = get_order(dx[1])
    aux = δtI^interval(orderT+1)
    N = length(x)
    Δx  = [xI[i][orderT+1] for i in eachindex(xI)] * aux
    Δdx = [dxI[i][orderT+1] for i in eachindex(xI)] * aux
    Δ0  = [dx[i][orderT](δI) for i in eachindex(x)] * aux / (orderT+1)
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
                if isequal_interval(Δ[ind], Δx[ind])
                    # vv[ind] = widen.(Δ[ind])
                    vv[ind] = interval(prevfloat(inf(Δ[ind])), nextfloat(sup(Δ[ind])))
                end
            end
            Δx = copy(vv)
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
#         if issubset_interval(Δ, Δx)
#             @inbounds for ind in 1:N
#                 # Widen the directions where ⊂ does not hold
#                 vv[ind] = Δx[ind]
#                 if issubset_interval(Δx[ind], Δ[ind])
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
    picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0, params)

Return the remainder of Picard operator
"""
function picard_remainder!(f!::Function, t::Taylor1{T},
    x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
    xxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    δI::AbstractVector{Interval{T}}, δt::Interval{T},
    Δx::Vector{Interval{T}}, Δ0::Vector{Interval{T}}, params) where {T}

    # Extend `x` and `dx` to have interval coefficients
    # zI = zero(δt)
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
