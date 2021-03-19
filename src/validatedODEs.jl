# Some methods for validated integration of ODEs

"""
    remainder_taylorstep!(f!, t, x, dx, xI, dxI, δI, δt, params)

Returns a remainder for the integration step for the dependent variables (`x`)
checking that the solution satisfies the criteria for existence and uniqueness.
"""
function remainder_taylorstep!(f!::Function, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}},
        δI::IntervalBox{N,T}, δt::Interval{T}, params) where {N,T}

    orderT = get_order(dx[1])
    aux = δt^(orderT+1)
    Δx  = IntervalBox([xI[i][orderT+1] for i in eachindex(xI)]) * aux
    Δdx = IntervalBox([dxI[i][orderT+1] for i in eachindex(xI)]) * aux
    Δ0  = IntervalBox([dx[i][orderT](δI) for i in eachindex(x)]) * aux / (orderT+1)
    Δ = Δ0 + Δdx * δt
    Δxold = Δx

    # Checking existence and uniqueness
    iscontractive(Δ, Δx) && return Δx

    # If the check didn't work, compute new remainders. A new Δx is proposed,
    # and the corresponding Δdx is computed
    xxI  = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
    dxxI = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
    vv = Array{Interval{T}}(undef, N)
    for its = 1:50
        # Remainder of Picard iteration
        Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0, params)

        # Checking existence and uniqueness
        iscontractive(Δ, Δx) && return Δx
        # iscontractive(Δ, Δx) && return _contract_iteration!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δdx, Δ0, params)

        # Expand Δx in the directions needed
        Δxold = Δx
        if Δ ⊆ Δx
            @inbounds for ind in 1:N
                # Widen the directions where ⊂ does not hold
                vv[ind] = Δx[ind]
                if Δ[ind] == Δx[ind]
                    vv[ind] = widen.(Δ[ind])
                end
            end
            Δx = IntervalBox(vv)
            continue
        end
        Δx = Δ
    end

    # If it didn't work, throw an error
    @format full
    error("""
    Error: it cannot prove existence and unicity of the solution:
        t0 = $(t[0])
        δt = $(δt)
        Δ  = $(Δ)
        Δxo = $(Δxold)
        Δx = $(Δx)
        $(Δ .⊆ Δxold)
    """)
end


"""
    iscontractive(Δ, Δx)

Checks if `Δ .⊂ Δx` is satisfied. If ``Δ ⊆ Δx` is satisfied, it returns
`true` if all cases where `==` holds corresponds to the zero `Interval`.
"""
function iscontractive(Δ::IntervalBox{N,T}, Δx::IntervalBox{N,T}) where{N,T}
    zI = zero_interval(T)
    @inbounds for ind in 1:N
        Δ[ind] ⊂ Δx[ind] && continue
        Δ[ind] == Δx[ind] == zI && continue
        return false
    end
    return true
end


"""
    picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0, params)

Return the remainder of Picard operator
"""
function picard_remainder!(f!::Function, t::Taylor1{T},
    x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
    xxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    δI::IntervalBox{N,T}, δt::Interval{T},
    Δx::IntervalBox{N,T}, Δ0::IntervalBox{N,T}, params) where {N,T}

    # Extend `x` and `dx` to have interval coefficients
    zI = zero_interval(T)
    @inbounds for ind in eachindex(x)
        xxI[ind]  = x[ind] + Δx[ind]
        dxxI[ind] = dx[ind] + zI
    end

    # Compute `dxxI` from the equations of motion
    f!(dxxI, xxI, params, t)

    # Picard iteration, considering only the bound of `f` and the last coeff of f
    Δdx = IntervalBox( evaluate.( (dxxI - dx)(δt), δI ) )
    Δ = Δ0 + Δdx * δt
    return Δ
end


# Picard iterations to contract further Δx, once Δ ⊂ Δx holds
# **Currently not used**
function _contract_iteration!(f!::Function, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        xxI::Vector{Taylor1{TaylorN{Interval{T}}}}, dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
        δI::IntervalBox{N,T}, δt::Interval{T},
        Δx::IntervalBox{N,T}, Δdx::IntervalBox{N,T}, Δ0::IntervalBox{N,T}, params) where {N,T}

    # Some abbreviations
    Δ = Δ0 + Δdx * δt
    Δxold = Δx

    # Picard contractions
    for its = 1:10
        # Remainder of Picard iteration
        Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0, params)

        # If contraction doesn't hold, return old bound
        iscontractive(Δ, Δx) || return Δxold

        # Contract estimate
        Δxold = Δx
        Δx = Δ
    end

    return Δxold
end


"""
    absorb_remainder(a::TaylorModelN)

Returns a TaylorModelN, equivalent to `a`, such that the remainder
is mostly absorbed in the constant and linear coefficients. The linear shift assumes
that `a` is normalized to the `IntervalBox(-1..1, Val(N))`.

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
    rem = Interval{T}(0, 0)

    # Linear shift
    lin_shift = mid(Δ) + sum((aux*TaylorN(i, order=orderQ) for i in 1:N))
    bpol = a.pol + lin_shift

    # Compute the new remainder
    aI = a(δ)
    bI = bpol(δ)

    if bI ⊆ aI
        rem = Interval(aI.lo-bI.lo, aI.hi-bI.hi)
    elseif aI ⊆ bI
        rem = Interval(bI.lo-aI.lo, bI.hi-aI.hi)
    else
        r_lo = aI.lo-bI.lo
        r_hi = aI.hi-bI.hi
        if r_lo > 0
            rem = Interval(-r_lo, r_hi)
        else
            rem = Interval( r_lo, -r_hi)
        end
    end

    return TaylorModelN(bpol, rem, a.x0, a.dom)
end


# Postverify and define Taylor models to be returned
function scalepostverify_sw!(xTMN::Vector{TaylorModelN{N,T,T}},
        X::Vector{TaylorN{T}}) where {N,T}
    postverify = true
    x0 = xTMN[1].x0
    B = domain(xTMN[1])
    zI = zero(B[1])
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
    zI = zero_interval(T)
    x0 = zero_box(N, T)
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
    @assert invjac * qB ⊆ qB´

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
    eq16 = all(mag.(evaluate.(jaq_q1, (q .* B,))) .≤ s)
    if eq16
        postverify = scalepostverify_sw!(xTMN, q .* X)
        postverify && return q
    end
    s .= eps.(q)
    @. q = 1.0 + r̃ + s
    jaq_q1 = jacmatrix_g * (q .- 1.0)
    eq16 = all(mag.(evaluate.(jaq_q1, (q .* B,))) .≤ s)
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
        eq16 = all(evaluate.(jaq_q1, (qB,)) .≤ s)
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


"""
    validated-step!
"""
function validated_step!(f!, t::Taylor1{T}, x::Vector{Taylor1{TaylorN{T}}},
        dx::Vector{Taylor1{TaylorN{T}}}, xaux::Vector{Taylor1{TaylorN{T}}},
        tI::Taylor1{T}, xI::Vector{Taylor1{Interval{T}}},
        dxI::Vector{Taylor1{Interval{T}}}, xauxI::Vector{Taylor1{Interval{T}}},
        t0::T, tmax::T, sign_tstep::Int,
        xTMN::Vector{TaylorModelN{N,T,T}}, xv::Vector{IntervalBox{N,T}},
        rem::Vector{Interval{T}}, zbox::IntervalBox{N,T}, symIbox::IntervalBox{N,T},
        nsteps::Int, orderT::Int, abstol::T, params, parse_eqs::Bool,
        check_property::Function=(t, x)->true) where {N,T}

    # One step integration (non-validated)
    # TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux, params)
    # δt = TaylorIntegration.stepsize(x, abstol)
    δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux, abstol, params, parse_eqs)
    f!(dx, x, params, t)  # Update `dx[:][orderT]`

    # One step integration for the initial box
    # TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, tI, xI, dxI, xauxI, params)
    # δtI = TaylorIntegration.stepsize(xI, abstol)
    δtI = TaylorIntegration.taylorstep!(f!, tI, xI, dxI, xauxI, abstol, params, parse_eqs)
    f!(dxI, xI, params, tI)  # Update `dxI[:][orderT+1]`

    # Step size
    δt = min(δt, sign_tstep*(tmax-t0))
    δt = sign_tstep * δt

    # Test if `check_property` is satisfied; if not, half the integration time.
    # If after 25 checks `check_property` is not satisfied, throw an error.
    nsteps += 1
    issatisfied = false
    rem_old = copy(rem)
    for nchecks = 1:25
        # Validate the solution: remainder consistent with Schauder thm
        δtI = sign_tstep * Interval(0, sign_tstep*δt)
        Δ = remainder_taylorstep!(f!, t, x, dx, xI, dxI, symIbox, δtI, params)
        rem .= rem_old .+ Δ

        # Create TaylorModelN to store remainders and evaluation
        @inbounds begin
            for i in eachindex(x)
                xTMN[i] = fp_rpa( TaylorModelN(x[i](δtI), rem[i], zbox, symIbox) )

                # If remainder is still too big, do it again
                j = 0
                while (j < 10) && (mag(rem[i]) > 1.0e-10)
                    j += 1
                    xTMN[i] = absorb_remainder(xTMN[i])
                    rem[i] = remainder(xTMN[i])
                end
            end
            xv[nsteps] = evaluate(xTMN, symIbox) # IntervalBox

            if !check_property(t0+δt, xv[nsteps])
                δt = δt/2
                continue
            end
        end # @inbounds

        issatisfied = true
        break
    end

    if !issatisfied
        error("""
            `check_property` is not satisfied:
            $t0 $nsteps $δt
            $(xv[nsteps])
            $(check_property(t0+δt, xv[nsteps]))""")
    end

    return δt
end

"""
    initialize!(X0::IntervalBox{N, T}, orderQ, orderT, x, dx, xTMN, xI, dxI, rem, xTM1v) where {N, T}

Initialize the auxiliary integration variables and normalize the given interval
box to the domain `[-1, 1]^n`.
"""
function initialize!(X0::IntervalBox{N, T}, orderQ, orderT, x, dx, xTMN, xI, dxI, rem, xTM1v) where {N, T}
    @assert N == get_numvars()

    # center of the box and vector of widths
    q0 = mid.(X0)
    δq0 = X0 .- q0

    # normalized domain
    zI = zero_interval(T)
    zB = zero_box(N, T)
    S = symmetric_box(N, T)

    @inbounds for i in eachindex(x)
        qaux = normalize_taylor(q0[i] + TaylorN(i, order=orderQ), δq0, true)
        x[i] = Taylor1(qaux, orderT)
        dx[i] = x[i]
        xTMN[i] = TaylorModelN(qaux, zI, zB, S)

        xI[i] = Taylor1(X0[i], orderT+1)
        dxI[i] = xI[i]
        rem[i] = zI

        xTM1v[i, 1] = TaylorModel1(deepcopy(x[i]), zI, zI, zI)
    end
end

"""
    initialize!(X0::Vector{TaylorModel1{TaylorN{T}, T}}, orderQ, orderT, x, dx, xTMN, xI, dxI, rem, xTM1v) where {T}

Initialize the auxiliary integration variables assuming that the given vector
of taylor models `X0` is normalized to the domain `[-1, 1]^n` in space.
"""
function initialize!(X0::Vector{TaylorModel1{TaylorN{T}, T}}, orderQ, orderT, x, dx, xTMN, xI, dxI, rem, xTM1v) where {T}
    # nomalized domain
    N = get_numvars()
    zI = zero_interval(T)
    zB = zero_box(N, T)
    S = symmetric_box(N, T)

    @inbounds for i in eachindex(x)
        yi = X0[i]
        pi = polynomial(yi)

        # we only keep the t^0 coefficient
        qaux = pi.coeffs[1]
        @assert all(iszero, pi.coeffs[2:end])

        x[i] = Taylor1(qaux, orderT)
        dx[i] = x[i]
        xTMN[i] = TaylorModelN(qaux, zI, zB, S)

        # we assume that qaux is normalized to [-1, 1]^N
        pi_int = evaluate(qaux, S)
        xI[i] = Taylor1(pi_int, orderT+1)
        dxI[i] = xI[i]

        # remainder
        rem[i] = remainder(yi)

        # expansion point in time assumed zero
        x0t = expansion_point(yi)
        @assert x0t == zI

        # domain in time assumed zero
        domt = domain(yi)
        @assert domt == zI
        xTM1v[i, 1] = TaylorModel1(deepcopy(x[i]), rem[i], x0t, domt)
    end
end

function validated_integ(f!, X0, t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T,
                         params=nothing; maxsteps::Int=500, parse_eqs::Bool=true,
                         check_property::Function=(t, x)->true) where {T<:Real}

    # Set proper parameters for jet transport
    N = get_numvars()
    dof = N

    # Some variables
    zI = zero_interval(T)
    zB = zero_box(N, T)
    S  = symmetric_box(N, T)
    t  = t0 + Taylor1(orderT)
    tI = t0 + Taylor1(orderT+1)

    # Allocation of vectors

    # Output
    tv    = Array{T}(undef, maxsteps+1)
    xv    = Array{IntervalBox{N,T}}(undef, maxsteps+1)
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    rem = Array{Interval{T}}(undef, dof)

    # Internals: jet transport integration
    x     = Array{Taylor1{TaylorN{T}}}(undef, dof)
    dx    = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xaux  = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xTMN  = Array{TaylorModelN{N,T,T}}(undef, dof)

    # Internals: Taylor1{Interval{T}} integration
    xI    = Array{Taylor1{Interval{T}}}(undef, dof)
    dxI   = Array{Taylor1{Interval{T}}}(undef, dof)
    xauxI = Array{Taylor1{Interval{T}}}(undef, dof)

    # Set initial conditions
    initialize!(X0, orderQ, orderT, x, dx, xTMN, xI, dxI, rem, xTM1v)
    sign_tstep = copysign(1, tmax-t0)

    # Output vectors
    @inbounds tv[1] = t0
    @inbounds xv[1] = evaluate(xTMN, S)

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx, params)
        catch
            parse_eqs = false
        end
    end

    # Integration
    nsteps = 1
    while sign_tstep*t0 < sign_tstep*tmax

        # Validated step of the integration
        δt = validated_step!(f!, t, x, dx, xaux, tI, xI, dxI, xauxI,
            t0, tmax, sign_tstep, xTMN, xv, rem, zB, S,
            nsteps, orderT, abstol, params, parse_eqs, check_property)

        # New initial conditions and time
        nsteps += 1
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tI[0] = t0
        @inbounds tv[nsteps] = t0
        @inbounds for i in eachindex(x)
            δtI = sign_tstep * Interval{T}(0, sign_tstep*δt)
            xTM1v[i, nsteps] = TaylorModel1(deepcopy(x[i]), rem[i], zI, δtI)
            aux = x[i](δt)
            x[i]  = Taylor1(aux, orderT)
            dx[i] = Taylor1(zero(aux), orderT)
            auxI = xTMN[i](S)
            xI[i] = Taylor1(auxI, orderT+1)
            dxI[i] = xI[i]
        end

        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end

    end

    return view(tv,1:nsteps), view(xv,1:nsteps), view(xTM1v, :, 1:nsteps)
end

"""
    picard(dx, x0, box)

Computes the picard (integral) operator for the initial condition `x0`.
`dx` must be the rhs of the differential equation.
"""
function picard(dx, x0, box)
    ∫f = integrate(dx, 0., box)
    pol = ∫f.pol + x0.pol # picard operator
    return TaylorModel1(deepcopy(pol), ∫f.rem + x0.rem, ∫f.x0, ∫f.dom)
end

function _picard(dx, x0, box)
    ∫f = integrate(dx, 0., box)
    pol = ∫f.pol + x0.pol
    Δk = ∫f.rem
    return pol, Δk
end

"""
    verify_contraction(f!, dx, xTM1K, params, t, x0, box)

Verifies contraction property for a given set of `f!`, `dx`, `xTM1K`,
`params`, `t`, `x0`, `box`. It asserts the contention remainder(∫f!(dx, xTM1K, params, t)) ⊆ remainder(xTM1K).
"""
function verify_contraction(f!, dx, xTM1K, params, t, x0, box)
    f!(dx, xTM1K, params, t)
    @inbounds for i in eachindex(dx)
        pii, Δ = _picard(dx[i], x0[i], box)
        @assert (Δ + x0[i].rem) ⊆ xTM1K[i].rem
    end
end

function _bound_integration(a::Taylor1, Δr, δ, δI)
    order = get_order(a)
    aux = δ^order / (order+1)
    Δ = δ * (Δr + getcoeff(a, order)(δI) * aux)
    return Δ
end

function _validate_step!(xTM1K, f!, dx, x0, params, t, box, dof; ε=1e-10, δ=1e-5,
                      maxsteps=20, extrasteps=50)
    E = [interval(0) for i in eachindex(x0)]
    E′ = [interval(0) for _ in 1:dof]
    polv = polynomial.(xTM1K)
    x0v = expansion_point.(xTM1K)
    domv = domain.(xTM1K)
    low_ratiov = Array{Float64, 1}(undef, dof)
    hi_ratiov = Array{Float64, 1}(undef, dof)
    nsteps = 0
   
    while nsteps < maxsteps
        nsteps += 1
        f!(dx, xTM1K, params, t)
        @inbounds for i in eachindex(dx)
            pᵢ₁, Δ = _picard(dx[i], x0[i], box)
            E′[i] = Δ + x0[i].rem
        end

        # Only inflates the required component
        @inbounds for i in eachindex(dx)
            if E′[i] ⊆ E[i] || E′[i] == E[i] == Interval(0)
                xTM1K[i] = TaylorModel1(polv[i], E′[i], x0v[i], domv[i])
            else
                εi = (1 - ε) .. (1 + ε)
                δi = -δ .. δ
                E[i] = E′[i] * εi + δi
                xTM1K[i] = TaylorModel1(polv[i], E[i], x0v[i], domv[i])
            end
        end

        all(E′ .⊆ E) && break

    end
    
    @inbounds for i in eachindex(dx)
        low_ratiov[i] = E′[i].lo / E[i].lo
        hi_ratiov[i] = E′[i].hi / E[i].hi
    end
    
    # Contract further the remainders if the last contraction improves more than 5%
    for ind = 1:extrasteps
        all(low_ratiov .> 0.95) && all(hi_ratiov .> 0.95) && break
        f!(dx, xTM1K, params, t)
        @inbounds for i in eachindex(dx)
            E[i] = E′[i]
            _, Δ = _picard(dx[i], x0[i], box)
            E′[i] = Δ + x0[i].rem
            xTM1K[i] = TaylorModel1(polv[i], E′[i], x0v[i], domv[i])
            low_ratiov[i] = E′[i].lo / E[i].lo
            hi_ratiov[i] = E′[i].hi / E[i].hi
        end
    end


    verify_contraction(f!, dx, xTM1K, params, t, x0, box)

    if nsteps == maxsteps
        @warn ("Maximum number of validate steps reached.")
    end

    return nothing
end

function validated_integ2(f!, qq0, δq0::IntervalBox{N, T}, t0, tf, orderQ, orderT,
                         abstol, params=nothing; parse_eqs=true, maxsteps=500,
                         validatesteps=30, ε=1e-10, δ=1e-3,
                         absorb_steps=3) where {N, T}
    dof = N
    @assert N == get_numvars()
    zI = zero(Interval{T})
    zbox = IntervalBox(zI, Val(N))
    symIbox = IntervalBox(Interval{T}(-1 .. 1), Val(N))
    q0 = IntervalBox(qq0)
    t = t0 + Taylor1(orderT)
    
    tv = Array{T}(undef, maxsteps+1)
    xv = Array{IntervalBox{N, T}}(undef, maxsteps+1)
    xTM1v = Array{TaylorModel1{TaylorN{T}, T}}(undef, dof, maxsteps+1)
    x = Array{Taylor1{TaylorN{T}}}(undef, dof)
    dx = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xaux = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xTMN = Array{TaylorModelN{N, T, T}}(undef, dof)
    dxTM1 = Array{TaylorModel1{TaylorN{T}, T}}(undef, dof)
    xTM1 = Array{TaylorModel1{TaylorN{T}, T}}(undef, dof)
    xTM1K = Array{TaylorModel1{TaylorN{T}, T}}(undef, dof)

    rem = Array{Interval{T}}(undef, dof)

    @inbounds for i in eachindex(x)
        qaux = normalize_taylor(qq0[i] + TaylorN(i, order=orderQ), δq0, true)
        x[i] = Taylor1(qaux, orderT)
        dx[i] = x[i]
        xTMN[i] = TaylorModelN(qaux, zI, zbox, symIbox)
        rem[i] = zI
        xTM1v[i, 1] = TaylorModel1(deepcopy(x[i]), zI, zI, zI)
    end

    sign_tstep = copysign(1, tf - t0)
    nsteps = 1
    @inbounds tv[1] = t0
    @inbounds xv[1] = evaluate(xTMN, symIbox)

    parse_eqs = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    while t0 * sign_tstep < tf * sign_tstep
        δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux, abstol, params, parse_eqs)
        f!(dx, x, params, t)

        δt = min(δt, sign_tstep*(tf-t0))
        δt = sign_tstep * δt

        @inbounds for i in eachindex(x)
            dom = sign_tstep > 0 ? 0 .. δt : δt .. 0
            x0 = sign_tstep > 0 ? dom.lo : dom.hi
            Δ = zero(Interval{Float64})
            xTM1[i] = TaylorModel1(deepcopy(x[i]), Δ, x0, dom)
        end
        
        # to reuse the previous TaylorModel and save some allocations
        _validate_step!(xTM1, f!, dxTM1, xTMN, params, t, symIbox, dof,
                        maxsteps=validatesteps, ε=ε, δ=δ)
        t0 += δt
        nsteps += 1
        
        @inbounds t[0] = t0
        @inbounds tv[nsteps] = t0
        xv[nsteps] = evaluate(xTMN, symIbox)

        @inbounds for i in eachindex(x)
            aux_pol = evaluate(xTM1[i], Interval(δt))
            rem[i] = remainder(xTM1[i])
            xTMN[i] = fp_rpa(TaylorModelN(deepcopy(aux_pol), 0 .. 0, zbox, symIbox))
            # temporal solution
            j = 0
            while (j < absorb_steps) && (mag(rem[i]) > 1.0e-10)
                j += 1
                xTMN[i] = absorb_remainder(xTMN[i])
                rem[i] = remainder(xTMN[i])
            end
            x[i] = Taylor1(xTMN[i].pol, orderT)
            xTM1v[i, nsteps] = xTM1[i]
        end

        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    return view(tv, 1:nsteps), view(xv, 1:nsteps), view(xTM1v, :, 1:nsteps)
end
