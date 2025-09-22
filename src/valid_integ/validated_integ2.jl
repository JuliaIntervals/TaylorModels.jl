# Some methods for validated integration of ODEs (second approach)

function validated_integ2(f!, X0::AbstractVector{Interval{U}}, t0::T, tmax::T, orderQ::Int, orderT::Int,
        abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        validatesteps::Int=30, ε::T=1e-10, δ::T=1e-6,
        absorb_steps::Int=3) where {T <: Real, U}

    # Initialize cache
    vX0 = Vector(X0)
    cacheVI = init_cache_VI(t0, vX0, maxsteps, orderT, orderQ, f!, params; parse_eqs)

    return _validated_integ2!(f!, vX0, t0, tmax, abstol, cacheVI, params,
        maxsteps, adaptive, minabstol, absorb,
        validatesteps, ε, δ, absorb_steps)
end

function validated_integ2(f!, X0::Vector{TaylorModel1{TaylorN{T}, U}}, t0::T, tmax::T, orderQ::Int, orderT::Int,
        abstol::T, params=nothing;
        maxsteps::Int=2000, parse_eqs::Bool=true,
        adaptive::Bool=true, minabstol::T=T(_DEF_MINABSTOL), absorb::Bool=false,
        validatesteps::Int=30, ε::T=1e-10, δ::T=1e-6,
        absorb_steps::Int=3) where {T <: Real, U}

    # Initialize cache
    cacheVI = init_cache_VI(t0, X0, maxsteps, orderT, orderQ, f!, params; parse_eqs)
    q0 = evaluate(constant_term.(polynomial.(X0)), Vector(symmetric_box(length(X0),U)))

    return _validated_integ2!(f!, q0, t0, tmax, abstol, cacheVI, params,
        maxsteps, adaptive, minabstol, absorb,
        validatesteps, ε, δ, absorb_steps)
end


function _validated_integ2!(f!, q0, t0::T, tf::T, abstol::T, cacheVI::VectorCacheVI, params,
        maxsteps::Int, adaptive::Bool, minabstol::T, absorb::Bool,
        validatesteps::Int, ε::T, δ::T, absorb_steps::Int) where {T <: Real}

    # Unpack caches
    @unpack tv, xv, xaux, t, x, dx, rv, xauxI, tI, xI, dxI, rvI,
            xTMN, xTM1v, rem, parse_eqs = cacheVI

    # Set proper parameters for jet transport
    sign_tstep = copysign(1, tf - t0)
    dof = length(q0)
    orderT = get_order(t)
    zt = zero(t0)
    zI = interval(zero(T))
    S  = symmetric_box(dof, T)
    zB = zero(S)
    @inbounds xv[1] = q0
    @inbounds tv[1] = t0

    xTM1 = Array{TaylorModel1{TaylorN{T},T}}(undef, dof)
    dxTM1 = Array{TaylorModel1{TaylorN{T},T}}(undef, dof)
    low_ratiov = Array{T}(undef, dof)
    hi_ratiov = Array{T}(undef, dof)
    # rem = Array{Interval{T}}(undef, dof)
    E = Array{Interval{T}}(undef, dof)
    E′ = Array{Interval{T}}(undef, dof)

    # Initializations
    @. begin
    #     xTMN = TaylorModelN(constant_term(x), zI, (zB,), (S,))
        xTM1 = TaylorModel1(deepcopy(x), zI, zI, zI)
    #     rem = zI
    #     xTM1v[:, 1] = TaylorModel1(deepcopy(x), zI, zI, zI)
    end
    polv = polynomial.(xTM1)
    fill!(E, zI)
    fill!(E′, zI)
    # @inbounds tv[1] = t0
    # @inbounds xv[1] = evaluate(xTMN, S)

    # Direction of the integration
    # sign_tstep = copysign(1, tf - t0)

    red_abstol = abstol

    # Integration
    nsteps = 1
    while t0 * sign_tstep < tf * sign_tstep
        δt = TI.taylorstep!(Val(parse_eqs), f!, t, x, dx, xaux, abstol, params, rv)
        f!(dx, x, params, t)

        δt = min(δt, sign_tstep*(tf-t0))
        δt = sign_tstep * δt

        # Reuse previous TaylorModel1 to save some allocations
        (_success, δt, red_abstol) = _validate_step!(xTM1, f!, dxTM1, xTMN, params, x, t,
                                    S, dof, rem, red_abstol, δt, sign_tstep, E, E′,
                                    polv, low_ratiov, hi_ratiov,
                                    adaptive, minabstol,
                                    ε=ε, δ=δ,
                                    validatesteps=validatesteps)
        domt = sign_tstep * interval(zt, sign_tstep*δt)

        # δtI = intersect_interval(interval(δt, δt), domt) # assure it is inside the domain in t
        nsteps += 1
        @inbounds tv[nsteps] = t0
        t0 += δt
        @inbounds t[0] = t0

        # Flowpipe
        @. begin
            rem = remainder(xTM1)
            xTMN = fp_rpa(TaylorModelN(copy(evaluate(xTM1, domt)), rem, (zB,), (S,)))
        end
        xv[nsteps] = evaluate(xTMN, S)

        # New initial condition
        @inbounds for i in eachindex(x)
            aux_pol = evaluate(xTM1[i], δt) #δtI
            # rem[i] = remainder(xTM1[i])
            xTMN[i] = fp_rpa(TaylorModelN(deepcopy(aux_pol), rem[i], zB, S))
            # xTMN[i] = fp_rpa(TaylorModelN(deepcopy(aux_pol), 0 .. 0, zB, S))

            # Absorb remainder
            j = 0
            while absorb && (j < absorb_steps) && (mag(rem[i]) > 1.0e-10)
                t[0] == 0 && println("absorb_remainder ")
                j += 1
                xTMN[i] = absorb_remainder(xTMN[i])
                rem[i] = remainder(xTMN[i])
            end

            x[i] = Taylor1(polynomial(xTMN[i]), orderT)
            xTM1v[i, nsteps] = xTM1[i]
        end

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
    _validate_step(xTM1K, f!, dx, x0, params, t, box, dof; ε=1e-10, maxsteps=20, extrasteps=50)

Validate the Taylor Model solution for the current integration time step.
This function implements the epsilon inflation algorithm proposed by Bünger
with some custom adaptations.

Ref: Florian B\"unger, A Taylor model toolbox for solving ODEs implemented in MATLAB/INTLAB,
J. Comput. Appl. Math. 368, 112511, https://doi.org/10.1016/j.cam.2019.112511
"""
function _validate_step!(xTM1K, f!, dx, x0, params, x, t, box, dof, rem, abstol,
                        δt, sign_tstep, E, E′, polv, low_ratiov, hi_ratiov,
                        adaptive::Bool, minabstol;
                        ε=1e-10, δ=1e-6, validatesteps=20, extrasteps=50)
    #
    T = IA.numtype(box[1])
    zI = interval(zero(T))
    zt = zero(t[0])
    domT = sign_tstep * interval(zt, sign_tstep*δt)
    orderT = get_order(t)
    @. begin
        polv = deepcopy.(x)
        xTM1K = TaylorModel1(polv, zI, zI, domT)
        # xTM1K = TaylorModel1(polv, rem, zI, domT)
        E = remainder(xTM1K)
        # E = remainder(x0)
    end
    εi = interval(1 - ε, 1 + ε)
    δi = interval(-δ, δ)

    _success = false
    reduced_abstol = abstol
    # Try to prove existence and uniqueness up to numchecks, including reducing the
    # time step
    bool_red = true
    # for nchecks = 1:numchecks
    while bool_red

        # Try to prove existence and uniqueness up to validatesteps
        nsteps = 0
        E′ .= picard_iteration(f!, dx, xTM1K, params, t, x0, box, Val(true)) # 0-th iteration
        @. xTM1K = TaylorModel1(polv, E′, zI, domT)
        while nsteps < validatesteps
            E′ .= picard_iteration(f!, dx, xTM1K, params, t, x0, box)
            all(iscontractive.(E′, E)) && break

            # Only inflates the required component
            @inbounds for i in eachindex(dx)
                if !iscontractive(E′[i], E[i])
                    E[i] = E′[i] * εi + δi
                end
                xTM1K[i] = TaylorModel1(polv[i], E[i], zI, domT)
            end
            nsteps += 1
        end
        _success = all(iscontractive.(E′, E))
        _success && break

        # Shrink stepsize `δt` if `adaptive` is `true` and `_success` is false,
        # up to some minimum
        if !_success
            if adaptive
                bool_red = reduced_abstol > minabstol
                if bool_red
                    reduced_abstol = reduced_abstol/10
                    δt = δt * 0.1^(1/orderT)
                    domT = sign_tstep * interval(0, sign_tstep*δt)
                    @. begin
                        xTM1K = TaylorModel1(polv, zI, zI, domT)
                        # xTM1K = TaylorModel1(polv, rem, zI, domT)
                        E = remainder(xTM1K)
                        # E = remainder(x0)
                    end
                else
                    @warn("Minimum absolute tolerance reached: ", t[0], E′, E,
                    _success, all(iscontractive.(E′, E)), reduced_abstol)
                end
            else
                @warn("""It cannot prove existence and unicity of the solution:
                    t0 = $(t[0])
                    δt = $(δt)
                    Δx = $(E)
                """)
                break
            end
        end
    end

    if !all(iscontractive.(E′, E))
        @warn("Maximum number of validate steps reached.", t[0], E′, E,
            _success, all(iscontractive.(E′, E)))
        return (_success, δt, reduced_abstol)
    end

    # @. begin
    #     low_ratiov = inf(E′) / inf(E)
    #     hi_ratiov  = sup(E′) / sup(E)
    # end

    # # Contract further the remainders if the last contraction improves more than 5%
    # for ind = 1:extrasteps
    #     minimum(low_ratiov) > 0.90 && minimum(hi_ratiov) > 0.90 && break
    #     E .= remainder.(xTM1K)
    #     E′ .= picard_iteration(f!, dx, xTM1K, params, t, x0, box)
    #     @. begin
    #         xTM1K = TaylorModel1(polv, E′, zI, dom)
    #         low_ratiov = inf(E′) / inf(E)
    #         hi_ratiov  = sup(E′) / sup(E)
    #     end
    #     # @show(ind, E, E′, all(iscontractive.(E′, E)), low_ratiov)
    # end

    return (_success, δt, reduced_abstol)
end



"""
    picard(dx, x0, box)

Computes the picard (integral) operator for the initial condition `x0`.
`dx` must be the rhs of the differential equation.
"""
function _picard(dx, x0, box)
    ∫f = integrate(dx, 0., box)
    pol = ∫f.pol + x0.pol
    Δk = remainder(∫f)
    return pol, Δk
end

_picard_rem(dx, box) = remainder(integrate(dx, box))


"""
    picard_iteration(f!, dx, xTM1K, params, t, x0, box, ::Val{true})

Computes the picard (integral) operator for the set of equations `f!` and the initial condition`x0`.
The Val parameter enables the selection of the desired method. This one returns the remainder of
the resulting Taylor Model with the remainder of the initial condition included.
"""
function picard_iteration(f!, dx, xTM1K, params, t, x0, box, ::Val{true})
    f!(dx, xTM1K, params, t)
    return @. _picard_rem(dx, (box,)) + remainder(x0)
end

"""
    picard_iteration(f!, dx, xTM1K, params, t, x0, box)

Computes the picard (integral) operator for the set of equations `f!` and the initial condition`x0`.
This method returns the remainder of the resulting Taylor Model without the remainder of the initial condition.
"""
function picard_iteration(f!, dx, xTM1K, params, t, x0, box)
    f!(dx, xTM1K, params, t)
    return _picard_rem.(dx, (box,))
end
