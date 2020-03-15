# Some methods for validated integration of ODEs

"""
    iscontractive(Δ, Δx)

Checks if `Δ .⊂ Δx` is satisfied. If ``Δ ⊆ Δx` is satisfied, it returns
`true` if all cases where `==` holds corresponds to the zero `Interval`.
"""
@inline function iscontractive(Δ::IntervalBox{N,T}, Δx::IntervalBox{N,T}) where{N,T}
    zi = Interval{T}(0, 0)
    @inbounds for ind in 1:N
        Δ[ind] ⊂ Δx[ind] && continue
        Δ[ind] == Δx[ind] == zi && continue
        return false
    end
    return true
end


# Shift remainders of each term of the TM1TMN and define a T1 polynomial with them
@inline function shift_reset_remainder!(
        x::TaylorModel1{TaylorModelN{N,T,T},T}) where{N,T}
    x0 = x.x0
    δtI = domain(x)
    zrem = zero(remainder(x))
    t_rem = Taylor1( remainder.(x[:]) )
    # x = TaylorModel1(Taylor1(polynomial(x)), zrem, x0, δtI)
    return t_rem
end


"""
    picard_remainder!(f!, dxTM1TMN, xTM1TMN, t, δtI, symIbox, params)

Returns the remainder of a Picard contraction, without performing the
integration, which is equivalent to assume that the polynomial part
is not changed.
"""
function picard_remainder!(f!,
        dx  :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        x   :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        t   :: Taylor1{T},
        δtI :: Interval{T},
        symIbox :: IntervalBox{N,T}, params) where {N,T}
    f!(dx, x, params, t)
    Δ_integration = bound_integration.(dx, δtI)
    Δx_picard = IntervalBox(evaluate.(Δ_integration, (symIbox,)))
    return Δx_picard
end


function update_remainder!(
        x  :: Vector{TaylorModel1{TaylorModelN{N,T,S},S}},
        Δx :: IntervalBox{N,S}) where {N,T,S}
    @inbounds x0 = x[1].x0
    @inbounds δtI = domain(x[1])
    @inbounds for ind in eachindex(x)
        x[ind] = TaylorModel1( polynomial(x[ind]), Δx[ind], x0, δtI)
    end
    return nothing
end

function validate_solution!(f!,
        xTM1TMN  :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        dxTM1TMN :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        xTMN     :: Vector{TaylorModelN{N,T,T}},
        t        :: Taylor1{T},
        symIbox  :: IntervalBox{N,T},
        sign_δt  :: Int,
        params = nothing;
        check_property :: Function = (t, x)->true,
        maxiter  :: Int = 50) where {N,T}

    # t_rem = Vector{Taylor1{Interval{T}}}(undef, N)
    # t_rem .= shift_reset_remainder!.(xTM1TMN)

    # Some abbreviations
    @inbounds x0 = xTM1TMN[1].x0
    @inbounds δtI = domain(xTM1TMN[1])
    δt = sign_δt * mag(δtI)
    Δx = IntervalBox(remainder.(xTM1TMN))
    Δx_picard = Δx
    solution_exists = false
    checkprop_satisfied = false

    # Adjust Δx to verify existence and uniqueness, including verifying
    # that check_property is fulfilled; reduce time if necessary.
    time_reduct = 1
    while time_reduct ≤ 40
        iter_picard = 1
        while iter_picard ≤ maxiter
            # Picard iteration; update remainder
            Δx_picard = picard_remainder!(f!, dxTM1TMN, xTM1TMN, t, δtI, symIbox, params)

            # If Δx_picard is too large, reduce the time interval
            any(isentire.(Δx_picard)) && break

            # If existence is validated, shrink remainder until a tolerance
            # condition is satisfied.
            solution_exists = iscontractive(Δx_picard, Δx)
            if solution_exists
                if maximum(mag.(Δx)) > eps(maximum(constant_term.(constant_term.(xTM1TMN))))
                    # ts0 = time()
                    iterate_picard_contraction!(f!, xTM1TMN, dxTM1TMN, t, symIbox, params)
                    # ts1 = time()
                    # println("iterate_picard_contraction! ", ts1-ts0)
                    Δx = IntervalBox(remainder.(xTM1TMN))
                end
                # @show(iter_picard)

                xTMN .= fp_rpa.(evaluate(xTM1TMN, (δtI,)))
                # @show(Δx, remainder.(xTMN))
                xv = evaluate(xTMN, symIbox)
                δt  = sign_δt * mag(δtI)
                checkprop_satisfied = check_property(t[0]+δt, xv)
                if !(checkprop_satisfied)
                    time_reduct += 1
                    δtI = δtI / 2
                    δt  = sign_δt * mag(δtI)
                    # Use half Δx
                    Δx = Δx / 2
                    update_remainder!(xTM1TMN, Δx)
                    if iter_picard ≤ maxiter && time_reduct ≤ 40
                        continue
                    else
                        break
                    end
                end

                # @show(time_reduct)
                # @show(t_rem(δtI))
                return δtI
            end

            # Adjust Δx and update `xTM1TMN`
            iter_picard += 1
            ff = 2 .- (Δx_picard .⊆ Δx)
            Δx = ff .* Δx_picard
            update_remainder!(xTM1TMN, Δx)
        end

        # If the Picard iterations did not converged, try again with half `δtI`
        time_reduct += 1
        δtI = δtI / 2
        δt  = sign_δt * mag(δtI)
        # Use half Δx
        Δx = Δx / 2
        update_remainder!(xTM1TMN, Δx)
    end

    if !solution_exists
        @format full
        error("""
        Error: it cannot prove existence and unicity of the solution:
            t0  = $(t[0])
            δtI = $(domain(xTM1TMN[1]))
            Δx  = $(IntervalBox(remainder.(xTM1TMN)))
            Δx_picard = $(Δx_picard)
            Δx_picard ⊆ Δx = $(Δx_picard ⊆ Δx)
        """)
    end

    if !checkprop_satisfied
        @format full
        error("""
        Error: `check_property` is not satisfied:
            $t0
            $δt
            xv
            check_property(t,xv)
            """)
    end

    # @show(t_rem(δtI))
    return δtI
end


function iterate_picard_contraction!(f!,
        xTM1TMN :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        dxTM1TMN:: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        t :: Taylor1{T}, symIbox :: IntervalBox{N,T}, params=nothing) where {N,T}

    # Some abbreviations
    @inbounds x0 = xTM1TMN[1].x0
    @inbounds δtI = domain(xTM1TMN[1])
    Δx = IntervalBox(remainder.(xTM1TMN))
    Δx_picard = IntervalBox(remainder.(xTM1TMN))

    j = 0
    while true
        if maximum(mag.(Δx)) < eps(maximum(constant_term.(constant_term.(xTM1TMN))))
            # @show(j)
            return nothing
        end

        j += 1
        Δx = Δx_picard
        update_remainder!(xTM1TMN, Δx)

        # Picard iteration; update remainder
        Δx_picard = picard_remainder!(f!, dxTM1TMN, xTM1TMN, t, δtI, symIbox, params)
        if Δx_picard == Δx
            Δx = widen.(Δx)
            update_remainder!(xTM1TMN, Δx)
            # @show(j)
            return nothing
        end
    end
    return nothing
end


"""
    absorb_remainder!(a::Vector{TaylorModelN})  :: Vector{T}

Updates in place `a::TaylorModelN` in such a way that the remainder
is absorbed by the linear coefficients, without altering the other
terms of the polynomial. This is achived by the transformation
`X[i] -> X[i] + r δ_{ij}`, with `δ_{ij}` the Kronecker delta
and `X[i]` represents the TaylorN independent ith-variable.

The domain of all components of `a` is the
normalized inteval box `-1 .. 1`. It returns the vector `r` used
to achieve this

The idea of the proposed linear transformation
is somewhat similar to the shrink-wrapping method of Florian Bünger (see
[`shrink_wrapping!`](@ref)).
"""
function absorb_remainder!(a::Vector{TaylorModelN{N,T,T}}) where {N,T}
    # Original domain of TaylorModelN should be the symmetric normalized box
    B = IntervalBox(Interval{T}(-1,1), Val(N))
    @assert all(domain.(a) .== (B,))
    zI = Interval{T}(0, 0)
    @inbounds x0 = a[1].x0
    r = Array{T}(undef, N)
    X = Array{TaylorN{Interval{T}}}(undef, N)

    # Vector of independent TaylorN variables
    order = get_order(a[1])

    # Shift linear part of the polynomial and postverify
    @inbounds for i in eachindex(a)
        # TaylorN independent variable
        X[i] = TaylorN(Interval{T}, i, order=order)

        # Polynomial part of a[i], and sign of the i-th linear coefficient
        pol = polynomial(a[i])
        sign_linear_coeff = copysign(1, pol[1][i])

        # Remainder of original TaylorModelN and componentwise mag with
        # appropriate sign to achieve overapproximation
        rem = remainder(a[i])
        r[i] = sign_linear_coeff * mag(rem)

        # Scaled independent variable to perform the shift
        X[i] = r[i] * X[i]

        # New TaylorModelN
        ppol = fp_rpa( TaylorModelN(pol+X[i], zI, x0, B) )

        # Postverification
        bound_lin_a = (linear_polynomial(pol))(B) + remainder(a[i])
        bound_lin_b = (linear_polynomial(ppol))(B) + remainder(ppol)
        postverify = bound_lin_a ⊆ bound_lin_b
        # postverify = a[i](B) ⊆ ppol(B)
        @assert postverify """
        Failed to post-verify `absorb_remainder`
            i = $i
            X[i] = $(X[i])
            a[i](B) = $(a[i](B))
            ppol(B) = $(ppol(B))
            bound_lin_a ⊆ bound_lin_b = $(bound_lin_a ⊆ bound_lin_b)
            a[i](B) ⊆ ppol(B) = $(a[i](B) ⊆ ppol(B))
        """
        a[i] = copy(ppol)
    end

    return r
end


"""
    scale_postverify_sw!(xTMN, X)

Scales `xTMN::Vector{TaylorModelN{N,T,T}}` by `X::Vector{TaylorN{Interval{T}}}`
and sets its remainders to zero (shrink-wrapping method by Florian Bünger)
inplace. It returns a boolean that checks that the range of the original
`xTMN` is inside the range of the new `xTMN`, when evaluating them
in the full domains.
"""
function scale_postverify_sw!(xTMN::Vector{TaylorModelN{N,T,T}},
        X::Vector{TaylorN{Interval{T}}}) where {N,T}
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
    B = IntervalBox(Interval{T}(-1,1), Val(N))
    @assert domain(xTMN[1]) == B# all(domain.(xTMN) .== (B,))
    zI = zero(Interval{T})
    @inbounds x0 = xTMN[1].x0
    @inbounds order = get_order(xTMN[1])

    # Vector of independent TaylorN variables
    X = Array{TaylorN{Interval{T}}}(undef, N)
    qB = Array{Interval{T}}(undef, N)
    xTNcent = Array{TaylorN{T}}(undef, N)
    # X = [TaylorN(Interval{T}, i, order=order) for i in 1:N]

    # Initialize vectors
    @inbounds for i in eachindex(X)
        # TaylorN indep variables
        X[i] = TaylorN(Interval{T}, i, order=order)
        # Remainder of original TaylorModelN and componentwise mag
        rem = remainder(xTMN[i])
        r = mag(rem)
        qB[i] = r * B[i]
        # Shift to remove constant term
        xTNcent[i] = polynomial.(xTMN[i]) - constant_term(xTMN[i])
    end

    # Jacobian (at zero) and its inverse
    jac = TaylorSeries.jacobian(xTNcent)

    # If the conditional number is too large (jac is non-invertible),
    # don't change xTMN
    if cond(jac) > 1.0e4
        # @show(cond(jac))
        return ones(T, N)
    end

    # Inverse of the Jacobian
    invjac = inv(jac)

    # Componentwise bound
    r̃ = mag.(invjac * qB) # qB <-- r .* B
    qB´ = r̃ .* B
    @assert invjac * qB ⊆ qB´

    # Nonlinear part (linear part is close to identity)
    g = invjac*xTNcent .- X
    # ... and its jacobian matrix (full dependence!)
    jacmatrix_g = TaylorSeries.jacobian(g, X)

    # Step 7 of Bünger algorithm: obtain an estimate of `q`
    # Some constants/parameters
    q_tol = 1.0e-12
    q = 1.0 .+ r̃
    ff = 65/64
    q_max = 2 .* q
    zs = zero(q)
    s = similar(q)
    q_old = similar(q)
    q_1 = similar(q)
    jaq_q1 = jacmatrix_g * r̃

    iter_max = 100
    improve = true
    iter = 0
    while improve && iter < iter_max
        qB .= q .* B
        s .= zs
        q_old .= q
        q_1 .= q_old .- 1.0
        mul!(jaq_q1, jacmatrix_g, q_1)
        @inbounds for i in eachindex(xTMN)
            s[i] = mag( jaq_q1[i](qB) )
            q[i] = 1.0 + r̃[i] + s[i]
        end
        # If a component of q is too large, return the "trivial" q
        # any( q .> q_max ) && (@show(iter, s, q ./ q_max); return ones(T, N))
        improve = any( ((q .- q_old)./q) .> q_tol )
        iter += 1
    end
    # (improve || q == one.(q)) && return xTMN

    # Compute final q and rescale X
    @. q = 1.0 + r̃ + ff * s
    any( q .> q_max ) && return ones(T, N)
    # any( q .> q_max ) && (@show(iter, s, q ./ q_max); return ones(T, N))
    @. X = q * X

    # Scale TMNs and postverify
    postverify = scale_postverify_sw!(xTMN, X)
    @assert postverify "Failed to post-verify shrink-wrapping"

    return q
end


"""
    validated_step!
"""
function validated_step!(f!,
        t::Taylor1{T}, tp1::Taylor1{T},
        xT1TMN    :: Vector{Taylor1{TaylorModelN{N,T,T}}},
        dxT1TMN   :: Vector{Taylor1{TaylorModelN{N,T,T}}},
        xauxT1TMN :: Vector{Taylor1{TaylorModelN{N,T,T}}},
        xTM1TMN   :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        dxTM1TMN  :: Vector{TaylorModel1{TaylorModelN{N,T,T},T}},
        xT1TN     :: Vector{Taylor1{TaylorN{T}}},
        vv        :: Vector{TaylorN{T}},
        xTMN      :: Vector{TaylorModelN{N,T,T}},
        t0::T, tmax::T, sign_δt::Int,
        symIbox::IntervalBox{N,T},
        nsteps::Int, orderT::Int, abstol::T, params, parse_eqs::Bool,
        check_property::Function = (t, x)->true;
        max_checkprop::Int = 25) where {N,T}

    # One step integration (non-validated)
    # ts0 = time()
    TaylorIntegration.__jetcoeffs!(
        Val(parse_eqs), f!, tp1, xT1TMN, dxT1TMN, xauxT1TMN, params)

    # Reset auxiliary arrays
    for ind in eachindex(xT1TMN)
        for ord = 0:orderT
            @inbounds vv[ord+1] = polynomial(xT1TMN[ind][ord])
        end
        @inbounds xT1TN[ind] = Taylor1(copy(vv))
    end

    # Proposed time step
    zbox = zero(symIbox[1])
    δt = TaylorIntegration.stepsize(xT1TN, abstol)
    # ts1 = time()
    # println("TaylorIntegration.taylorstep! ", ts1-ts0)
    # ts0 = time()
    # TaylorIntegration.taylorstep!(f!, tp1, xT1TMN, dxT1TMN, xauxT1TMN, abstol, params, parse_eqs)
    # ts1 = time()
    # println("TaylorIntegration.taylorstep! ", ts1-ts0)
    δt  = min(δt, sign_δt*(tmax-t0))
    δtI = sign_δt * Interval(0, δt)

    # Construct xTM1TMN's of order `orderT`
    # ts0 = time()
    Δx = IntervalBox(remainder.(xTM1TMN))# Δx = IntervalBox(Interval{T}(0,0), N)
    # @show(Δx)
    @inbounds for ind in eachindex(xTM1TMN)
        xTM1TMN[ind] = TaylorModel1(
            Taylor1(copy(xT1TMN[ind].coeffs), orderT), Δx[ind], zbox, δtI)
    end
    # ts1 = time()
    # println("Construct xTM1TMN ", ts1-ts0)

    # Validate (existence and unicity) of the solution; the function updates
    # `xTM1TMN` which contains the valid domain δtI where the solution is proved
    # to exist and satisfy unicity, and the remainder.
    # ts0 = time()
    validate_solution!(f!, xTM1TMN, dxTM1TMN, xTMN, t, symIbox,
        sign_δt, params, check_property=check_property)
    # ts1 = time()
    # println("validate_solution! ", ts1-ts0)
    δtI = domain(xTM1TMN[1])

    return sign_δt * mag(δtI)
end


function validated_integ(f!, qq0::AbstractArray{T,1}, δq0::IntervalBox{N,T},
        t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T, params=nothing;
        maxsteps::Int=500, parse_eqs::Bool=true,
        check_property::Function=(t, x)->true) where {N, T<:Real}

    # Set proper parameters for jet transport
    @assert N == get_numvars()
    dof = N

    # Some variables
    R = Interval{T}
    zI = zero(R)
    zbox = IntervalBox(zI, Val(N))
    symIbox = IntervalBox(R(-1, 1), Val(N))
    sign_δt = copysign(1, tmax-t0)
    orderTp1 = orderT+1
    t   = t0 + Taylor1(orderT)
    tp1 = t0 + Taylor1(orderTp1)

    # Allocation of vectors
    # Output
    tout = Array{T}(undef, maxsteps+1)
    xout = Array{TaylorModel1{TaylorModelN{N,T,T},T}}(undef, dof, maxsteps+1)
    qout = Array{R}(undef, dof, maxsteps+1)
    #
    # Internals
    xT1TMN    = Array{Taylor1{TaylorModelN{N,T,T}}}(undef, dof)
    dxT1TMN   = Array{Taylor1{TaylorModelN{N,T,T}}}(undef, dof)
    xauxT1TMN = Array{Taylor1{TaylorModelN{N,T,T}}}(undef, dof)
    xTM1TMN   = Array{TaylorModel1{TaylorModelN{N,T,T},T}}(undef, dof)
    dxTM1TMN  = Array{TaylorModel1{TaylorModelN{N,T,T},T}}(undef, dof)
    xTMN      = Array{TaylorModelN{N,T,T}}(undef, dof)
    # xv        = Vector{R}(undef, dof)
    xT1TN    = Array{Taylor1{TaylorN{T}}}(undef, dof)
    vv       = Array{TaylorN{T}}(undef, orderTp1)
    rem      = Array{R}(undef, dof)
    q_absrem = Array{T}(undef, dof)
    q_shrink = Array{T}(undef, dof)
    # s_method = :shrink_wrapping
    s_method = :absorb_remainder

    # Set initial conditions and output vectors
    # ts0 = time()
    @inbounds tout[1] = t0
    @inbounds for i in eachindex(xTM1TMN)
        qaux = normalize_taylor(qq0[i] + TaylorN(i, order=orderQ), δq0, true)
        xT1TMN[i] = Taylor1(TaylorModelN(qaux, zI, zbox, symIbox), orderTp1)
        xTM1TMN[i] = TaylorModel1(Taylor1(copy(xT1TMN[i].coeffs), orderT), zI, zI, zI)
        rem[i] = zI
        q_absrem[i] = zero(T)
        q_shrink[i] = one(T)
        xout[i,1] = TaylorModel1(Taylor1(copy(xT1TMN[i].coeffs), orderT), zI, zI, zI)
        if s_method == :shrink_wrapping
            qout[i,1] = copy(q_shrink[i])
        else
            qout[i,1] = copy(q_absrem[i])
        end
    end
    # ts1 = time()
    # println("initial conditions ", ts1-ts0)

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx, params)
        catch
            parse_eqs = false
        end
    end
    xBout[1] = evaluate(evaluate(xTM1TMN, zI), symIbox)
    # ts1 = time()
    # println("initial conditions ", ts1-ts0)

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = TaylorIntegration._determine_parsing!(parse_eqs, f!, t, x, dx, params)

    # Integration
    nsteps = 1
    while sign_δt*t0 < sign_δt*tmax
        # @show(t[0])
        # Validated step of the integration
        # ts0 = time()
        δt = validated_step!(f!, t, tp1, xT1TMN, dxT1TMN, xauxT1TMN,
            xTM1TMN, dxTM1TMN, xT1TN, vv, xTMN,
            t0, tmax, sign_δt, symIbox,
            nsteps, orderT, abstol, params, parse_eqs, check_property)
        # ts1 = time()
        # println("validated_step! ", ts1-ts0)
        # @show(IntervalBox(remainder.(xTM1TMN)))

        # Update initial conditions and time, and populate output vectors
        nsteps += 1
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tout[nsteps] = t0
        δtI = domain(xTM1TMN[1])
        @inbounds for i in eachindex(xTM1TMN)
            qaux = fp_rpa( xTM1TMN[i](δt) )
            xTMN[i] = qaux
            rem[i] = remainder(qaux)
            xT1TMN[i] = Taylor1(qaux, orderTp1)
            xout[i, nsteps] = deepcopy(xTM1TMN[i])
            if s_method == :shrink_wrapping
                qout[i, nsteps] = copy(q_shrink[i])
            else
                qout[i, nsteps] = copy(q_absrem[i])
            end
        end

        # Use `absorb_remainder!` if necessary
        if maximum(mag.(rem)) > 1.0e-6
            # @show(rem)
            # ts0 = time()
            if s_method == :shrink_wrapping
                q_shrink .= shrink_wrapping!(xTMN)
            else
                q_absrem .= absorb_remainder!(xTMN)
            end
            # ts1 = time()
            @. rem = remainder(xTMN)
            # @show(rem)
            if s_method == :shrink_wrapping
                # println("shrink_wrapping! ", ts1-ts0)
                # @show(q_shrink)
            else
                # println("absorb_remainder! ", ts1-ts0)
                # @show(q_absrem)
            end
            @inbounds for i in eachindex(xTM1TMN)
                qaux = xTMN[i]
                rem[i] = remainder(qaux)
                xT1TMN[i] = Taylor1(qaux, orderTp1)
                if s_method == :shrink_wrapping
                    qout[i, nsteps] = copy(q_shrink[i])
                    q_shrink[i] = one(T)
                else
                    qout[i, nsteps] = copy(q_absrem[i])
                    q_absrem[i] = zero(T)
                end
            end
        end
        # println()

        # # Use `absorb_remainder!` if necessary
        # if maximum(mag.(rem)) > 1.0e-6
        #     # @show(rem)
        #     # ts0 = time()
        #     if s_method == :shrink_wrapping
        #         q_shrink .= shrink_wrapping!(xTMN)
        #     else
        #         q_absrem .= absorb_remainder!(xTMN)
        #     end
        #     # ts1 = time()
        #     @. rem = remainder(xTMN)
        #     # @show(rem)
        #     if s_method == :shrink_wrapping
        #         # println("shrink_wrapping! ", ts1-ts0)
        #         # @show(q_shrink)
        #     else
        #         # println("absorb_remainder! ", ts1-ts0)
        #         # @show(q_absrem)
        #     end
        #     @inbounds for i in eachindex(xTM1TMN)
        #         qaux = xTMN[i]
        #         rem[i] = remainder(qaux)
        #         xT1TMN[i] = Taylor1(qaux, orderTp1)
        #         if s_method == :shrink_wrapping
        #             qout[i, nsteps] = copy(q_shrink[i])
        #             q_shrink[i] = one(T)
        #         else
        #             qout[i, nsteps] = copy(q_absrem[i])
        #             q_absrem[i] = zero(T)
        #         end
        #     end
        # end
        # # println()

        if nsteps > maxsteps
            @warn("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end

    end

    return view(tout, 1:nsteps), view(xout, :, 1:nsteps), view(qout, :, 1:nsteps)
end
