# Methods related to validated integration

const _DEF_MINABSTOL = 1.0e-50


"""
    picard_lindelof(f!, dxTM1TMN::Vector{TaylorModel1{T,S}},
        xTM1TMN::Vector{TaylorModel1{T,S}}, params, t)

Compute the Picard-Lindelof operator to validate a solution
of a differential equation.
"""
function picard_lindelof(f!,
        dx1N::Vector{TaylorModel1{TaylorModelN{N,T,S},S}},
        x1N ::Vector{TaylorModel1{TaylorModelN{N,T,S},S}},
        t1N ::TaylorModel1{TaylorModelN{N,T,S},S},
        params) where {N,T,S}
    x2N = Vector{TaylorModel1{TaylorModelN{N,T,S},S}}
    z = zero(x1N[1])
    x2N .= z
    picard_lindelof!(f!, x2N, dx1N, x1N, params, t1N)
    return x2N
end

function picard_lindelof!(f!,
        x2N ::Vector{TaylorModel1{TaylorModelN{N,T,S},S}},
        dx1N::Vector{TaylorModel1{TaylorModelN{N,T,S},S}},
        x1N ::Vector{TaylorModel1{TaylorModelN{N,T,S},S}},
        t1N ::TaylorModel1{TaylorModelN{N,T,S},S},
        params) where {N,T,S}
    dom = domain(x1N[1][0])
    f!(dx1N, x1N, params, t1N)
    @inbounds for ind in eachindex(x1N)
        x2N[ind] = integrate(dx1N[ind], x1N[ind][0], dom)
    end
    return nothing
end


function _abs_rems!(vTMN, x0New, x1N, δt, symIbox, ind)
    tmn = evaluate(x1N[ind], δt)
    x0New[ind] = evaluate(tmn, symIbox)
    vTMN[ind] = fp_rpa(tmn)
    Δ = remainder(vTMN[ind])
    # Absorb remainder shifts
    vTMN[ind].pol.coeffs[1].coeffs[1] += mid(Δ)
    vTMN[ind].pol.coeffs[2].coeffs[ind] +=
        copysign(radius(Δ), vTMN[ind].pol[1].coeffs[ind])
    return nothing
end


function validated_integ3(f!,
        X0::AbstractVector{Interval{U}}, t0::T, tmax::T,
        orderQ::Int, orderT::Int, abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        check_property::F=(t, x)->true) where {T<:Real, U, F}

    # Initialize cache
    vX0 = Vector(X0)
    cacheVI = init_cache_VI3(t0, vX0, maxsteps, orderT, orderQ, f!, params;
        parse_eqs) :: VectorCacheVI3

    return _validated_integ3!(f!, vX0, t0, tmax, abstol, cacheVI, params,
        maxsteps, adaptive, minabstol, absorb, check_property)
end

function validated_integ3(f!,
        X0::Vector{TaylorModel1{TaylorModelN{N,T,U}, U}}, t0::T, tmax::T,
        orderQ::Int, orderT::Int, abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        check_property::F=(t, x)->true) where {N, T<:Real, U, F}

    # Initialize cache
    cacheVI = init_cache_VI3(t0, X0, maxsteps, orderT, orderQ, f!, params;
        parse_eqs) :: VectorCacheVI3
    q0 = evaluate(constant_term.(polynomial.(X0)), Vector(symmetric_box(length(X0),U)))

    return _validated_integ3!(f!, q0, t0, tmax, abstol, cacheVI, params,
        maxsteps, adaptive, minabstol, absorb, check_property)
end


function _validated_integ3!(f!, q0, t0::T, tmax::T, abstol::T,
        cacheVI::VectorCacheVI3, params,
        maxsteps::Int, adaptive::Bool, minabstol::T,
        absorb::Bool, check_property::F) where {T<:Real,F}

    # Unpack caches
    @unpack tv, xv, xaux, t, x, dx, rv,
            t1N, x1N, dx1N, x2N, z1N, vTN, vTMN,
            xTM1v, x0New, rem1, rem2, parse_eqs = cacheVI

    # Initial conditions
    sign_tstep = copysign(1, tmax - t0)
    dof = length(q0)
    orderT = get_order(t)
    orderQ = get_order(x[1][0])
    symIbox = symmetric_box(dof, T)
    zbox = zero(q0)
    @inbounds xv[1] = q0
    @inbounds tv[1] = t0
    zz = zero(x[1][0][0][1])

    # Integration
    nsteps = 1
    local _success # if true, the validation step succeeded
    red_abstol = abstol
    VV = Val(parse_eqs)
    while sign_tstep*t0 < sign_tstep*tmax
        # Polynomial approx and time step
        δt = validated_step3!(VV, f!,
            t, x, dx, xaux, rv,
            t0, tmax, sign_tstep,
            red_abstol, params)

        # Validated step of the integration
        (_success, δt, reduced_abstol) = _validation3(
            f!, t, x,
            t1N, x1N, dx1N, x2N, z1N,
            δt, sign_tstep,
            rem1, rem2, zbox, symIbox,
            orderT, abstol, params,
            adaptive, minabstol, absorb, check_property)

        # Output
        nsteps += 1
        cdom = centered_dom(x1N[1])
        tv[nsteps] = t0
        xv[nsteps] = copy(zbox)
        for ind in eachindex(x)
            _abs_rems!(vTMN, x0New, x1N, δt, symIbox, ind)
            xTM1v[ind, nsteps] = deepcopy(x1N[ind])
            xv[nsteps][ind] = evaluate(evaluate(x1N[ind], cdom), symIbox)
            #
            # Update initial state
            if dof == 1
                # No issue with the wrapping effect in 1-d
                normalize_taylorNs!(vTN, x0New, orderQ)
                TI.init_expansions!(x, dx, vTN, orderT)
            else
                # x[ind] = Taylor1(polynomial(vTMN[ind]), orderT)
                x1N[ind].rem = vTMN[ind].rem
                for ordT in eachindex(x1N[ind].pol.coeffs)
                    for ordQ in eachindex(x1N[ind].pol.coeffs[ordT].pol.coeffs)
                        for h in eachindex(x1N[ind].pol.coeffs[ordT].pol.coeffs[ordQ].coeffs)
                            x[ind].coeffs[ordT].coeffs[ordQ].coeffs[h] = zz
                            # dx[ind].coeffs[ordT].coeffs[ordQ].coeffs[h] = zz
                        end
                    end
                end
                # Constant coeff
                for ordQ in eachindex(x1N[ind].pol.coeffs[1].pol.coeffs)
                    for h in eachindex(x1N[ind].pol.coeffs[1].pol.coeffs[ordQ].coeffs)
                        x[ind].coeffs[1].coeffs[ordQ].coeffs[h] =
                            vTMN[ind].pol.coeffs[ordQ].coeffs[h]
                        dx[ind].coeffs[1].coeffs[ordQ].coeffs[h] = zz
                    end
                end
            end
        end
        # Update time
        t0 += δt
        t[0] = t0
        t1N.pol[0].pol[0][1] = t0

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
function validated_step3!(vB::Val{B}, f!,
        t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}},
        dx::Vector{Taylor1{TaylorN{T}}},
        xaux::Vector{Taylor1{TaylorN{T}}},
        rv::TI.RetAlloc{Taylor1{TaylorN{T}}},
        t0::T, tmax::T, sign_tstep::Int,
        abstol::T, params,) where {B,T}
    # One step integration (non-validated)
    δt = TI.taylorstep!(vB, f!, t, x, dx, xaux, abstol, params, rv)
    # f!(dx, x, params, t)  # Update last t coeff `dx[:][orderT]`
    # Step size
    δt = min(δt, sign_tstep*(tmax-t0))
    δt = sign_tstep * δt
    return δt
end


function _validation3(f!, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}},
        t1N::TaylorModel1{TaylorModelN{N,T,T},T},
        x1N::Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        dx1N::Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        x2N::Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        z1N::TaylorModel1{TaylorModelN{N,T,T},T},
        δt, sign_tstep::Int,
        rem1::Vector{Interval{T}},
        rem2::Vector{Interval{T}},
        zbox::AbstractVector{Interval{T}},
        symIbox::AbstractVector{Interval{T}},
        orderT::Int, abstol::T, params,
        adaptive::Bool, minabstol::T, absorb::Bool,
        check_property::Function=(t, x)->true) where {N,T}

    #
    local _success = false
    local bool_red = true
    local δtI
    reduced_abstol = abstol
    local issatisfied = false
    rem_old = copy(remainder.(x1N[:]))

    while bool_red
        # Verify Picard contraction
        δtI = sign_tstep * interval(0, sign_tstep*δt)
        z = zbox[1]
        t1N.dom = δtI
        z1N.dom = δtI
        # Reuse memory
        for i in eachindex(x1N)
            for ordT in eachindex(x1N[i].pol.coeffs)
                for ordQ in eachindex(x1N[i].pol.coeffs[ordT].pol.coeffs)
                    for h in eachindex(x1N[i].pol.coeffs[ordT].pol.coeffs[ordQ].coeffs)
                        x1N[i].pol.coeffs[ordT].pol.coeffs[ordQ].coeffs[h] =
                            x[i].coeffs[ordT].coeffs[ordQ].coeffs[h]
                    end
                end
                # x1N[i].pol.coeffs[ordT].rem = z
                # x1N[i].pol.coeffs[ordT].x0  = zbox
                # x1N[i].pol.coeffs[ordT].dom = symIbox
            end
            x1N[i].rem = rem_old[i]
            x1N[i].x0 = 0.0
            x1N[i].dom = δtI
            dx1N[i].pol = z1N.pol
            dx1N[i].rem = z1N.rem
            dx1N[i].x0 = x1N[i].x0
            dx1N[i].dom = δtI
            x2N[i].pol = z1N.pol
            x2N[i].rem = z1N.rem
            x2N[i].x0 = x1N[i].x0
            x2N[i].dom = δtI
            rem1[i] = z
            rem2[i] = z
        end
        for _ in 1:50
            rem1 .= total_remainder.(x1N)
            picard_lindelof!(f!, x2N, dx1N, x1N, t1N, params)
            rem2 .= total_remainder.(x2N)
            iscontractive(rem2, rem1) && break
            # x1N .= TaylorModel1.(x1N, 1.1 .* rem2)
            # # x1N .= TaylorModel1.(x2N, 1.1 .* remainder.(x2N))
            for i in eachindex(x1N)
                # x1N[i] = TaylorModel1(x1N[i], 1.1 * rem2[i])
                x1N[i].rem = 1.1 * rem2[i]
                # x1N[i].pol = copy(x1N[i].pol)
                # x1N[i].rem = 1.1 * remainder(x2N[i])
            end
        end
        _success = iscontractive(rem2, rem1)
        # @assert iscontractive(remainder.(x2N), remainder.(x1N))
        rem1 .= remainder.(x1N)

        # Shrink stepsize δt if adaptive is true and _success is false
        if !_success
            if @show(adaptive)
                bool_red = reduced_abstol > minabstol
                if bool_red
                    reduced_abstol = reduced_abstol/10
                    δt = δt / 2 #* 0.1^(1/orderT)
                    continue
                else
                    @warn("Minimum absolute tolerance reached: ",
                        t[0], δt, reduced_abstol, rem1)
                end
            else
                @warn("It cannot prove existence and unicity of the solution: ",
                    t[0], δt, rem1, _success)
            end
        end

        # xTMN .= evaluate.(x1N, δtI)

        # TO BE UPDATED
        # # Test if `check_property` is satisfied; if not, half the integration time.
        # # If after 25 checks `check_property` is not satisfied, throw an error.

        # # Create TaylorModelN to store remainders and evaluation
        # @inbounds begin
        #     for i in eachindex(x)
        #         # xTMN[i] = fp_rpa( TaylorModelN(x[i](δtI), rem[i], SVector{N}(zbox), SVector{N}(symIbox)) )
        #         xTMN[i] = TaylorModelN(x[i](δtI), rem[i], SVector{N}(zbox), SVector{N}(symIbox))

        #         # If remainder is still too big, do it again
        #         j = 0
        #         while absorb && (j < 10) && (mag(rem[i]) > 1.0e-10)
        #             t[0] == 0 && println("absorb_remainder ")
        #             j += 1
        #             xTMN[i] = absorb_remainder(xTMN[i])
        #             rem[i] = remainder(xTMN[i])
        #         end
        #     end
        #     xvv = evaluate(xTMN, symIbox) # interval box

        #     issatisfied = check_property(t[0]+δt, xvv)
        #     if !issatisfied
        #         # δt = δt/2
        #         bool_red = reduced_abstol > minabstol
        #         @info("issatisfied: ", bool_red, δt)
        #         if bool_red
        #             reduced_abstol = reduced_abstol / 10
        #             δt = δt * 0.1^(1/orderT)
        #             continue
        #         else
        #             @warn("Minimum absolute tolerance reached (issatisfied): ", t[0], δt, reduced_abstol, Δ)
        #         end
        #     end
        # end # @inbounds
            break
        end # while

    # if !issatisfied && !adaptive
    #     @warn("""
    #         Something went wrong:
    #         """, issatisfied, adaptive, bool_red, reduced_abstol,
    #         t[0], δt, Δ, xvv, check_property(t0+δt, xvv)
    #     )
    # end

    return (_success, δt, reduced_abstol)
end
