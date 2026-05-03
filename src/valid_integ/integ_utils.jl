"""
    shrink_wrapping!(xTMN::TaylorModelN)

Inplace modification of `xTMN`, which has absorbed the remainder
by the modified shrink-wrapping method of Florian Bünger.
The domain of `xTMN` is the normalized interval box `[-1,1]^N`.

Ref: Florian B\"unger, Shrink wrapping for Taylor models revisited,
Numer Algor 78:1001–1017 (2018), https://doi.org/10.1007/s11075-017-0410-1
"""
function shrink_wrapping!(xTMN::Vector{TaylorModelN{N,T,S}}) where {N,T,S}
    # Original domain of TaylorModelN should be the symmetric normalized box
    B = symmetric_box(S)
    @assert all(isequal_interval.(domain.(xTMN), (B,)))
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
    jac = TS.jacobian(xTNcent_lin)
    # If the conditional number is too large (inverse of jac is ill defined),
    # don't change xTMN; we use the mid-point jacobian matrix
    cond(mid.(jac)) > 1.0e4 && return one_r
    # Inverse of the Jacobian
    invjac = inv(jac)

    # Componentwise bound
    r̃ = mag.(invjac * qB) # qB <-- r .* B
    qBprime = r̃ .* B
    @assert issubset_interval(invjac * qB, qBprime)

    # Step 6 of Bünger algorithm: compute g
    g = invjac*xTNcent .- X
    # g = invjac*(xTNcent .- xTNcent_lin)
    # ... and its jacobian matrix (full dependence!)
    jacmatrix_g = TS.jacobian(g, X)

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
    jaq_q1 .= jacmatrix_g * (q .- 1.0)
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
        eq16 = all(mag.(evaluate.(jaq_q1, Ref(qB))) .≤ s)
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

# Postverify and define Taylor models to be returned
for TT in (:T, :(Interval{T}))
    @eval function scalepostverify_sw!(xTMN::Vector{TaylorModelN{N,$TT,S}},
            X::Vector{TaylorN{$TT}}) where {N,T<:IANumTypes, S<:IANumTypes}
        postverify = true
        x0 = expansion_point(xTMN[1])
        B = domain(xTMN[1])
        zI = zero(Interval{S})
        oI = one($TT)
        @inbounds for i in eachindex(xTMN)
            pol = polynomial(xTMN[i])
            tmn = TaylorModelN(pol(X), zI, x0, B )
            ppol = fp_rpa(tmn) * oI
            bb = issubset_interval(xTMN[i](B), ppol(B)) ||
                    isequal_interval(xTMN[i](B), ppol(B))
            postverify = postverify && bb
            xTMN[i] = copy(ppol)
        end
        @assert postverify """
            Failed to post-verify shrink-wrapping:
            X = $(linear_polynomial(X))
            xTMN = $(xTMN)
            """
        return postverify
    end
end


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
    δ = symmetric_box(T)
    aux = diam(Δ)/(2N)
    rem = zero(Δ)

    # Linear shift
    lin_shift = mid(Δ) + sum((aux*TaylorN(i, order=orderQ) for i in 1:N))
    bpol = polynomial(a) + lin_shift

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


function _abs_rems!(vTMN::Vector{TaylorModelN{N,T,S}}, x1N, auxN, δt, symIbox) where {N,T,S}
    for ind in eachindex(vTMN)
        # Evaluate x1N at δt (new TMN initial condition with remainder)
        TaylorModels.__evaluate!(vTMN[ind], x1N[ind], δt, auxN)
        Δ = remainder(vTMN[ind])
        radN = radius(Δ)/N
        # Old remainder of constant and linear parts and remainder
        aI = vTMN[ind].pol.coeffs[1].coeffs[1] +
            sum(vTMN[ind].pol.coeffs[2].coeffs) * symIbox[1] + Δ
        # New remainder of constant and linear parts (without Δ, a priori absorbed)
        bI = vTMN[ind].pol.coeffs[1].coeffs[1] + mid(Δ) +
            sum(vTMN[ind].pol.coeffs[2].coeffs .+ radN) * symIbox[1]
        # Compute the new remainder
        r_lo = copysign(inf(aI)-inf(bI), -1)
        r_hi = copysign(sup(aI)-sup(bI), 1)
        # Compare old and proposed new remainders; do nothing if new remainder is wider
        if r_hi-r_lo < diam(Δ)
            # Shifts to absorb remainders
            vTMN[ind].pol.coeffs[1].coeffs[1] += mid(Δ)
            for k in eachindex(vTMN[ind].pol.coeffs[2].coeffs)
                vTMN[ind].pol.coeffs[2].coeffs[k] += radN
            end
            # Store new remainder in TMN init condition
            vTMN[ind].rem = interval(r_lo, r_hi)
        end
    end
    return nothing
end


"""
    qrprecondition(vTMN::Vector{TaylorModelN{N,T,S}})

Returns the left and right preconditioned TaylorModelN's from `vTMN`,
following the explanation of Neher et al (2007).

Ref: M. Neher, K.R. Jackson and N.S. Nedialkov, "On Taylor Model based integration
of ODEs", SIAM J. NUMER. ANAL. 45 (1), pp. 236-262 (2007).
https://doi.org/10.1137/050638448
"""
function qrprecondition(vTMN::Vector{TaylorModelN{N,T,S}}) where {N,T,S}
    leftTMN = zero.(vTMN)
    rightTMN = zero.(vTMN)
    linTN = zero(Array{T}(undef, N, N))
    rems = zero.(remainder.(vTMN))
    scaleV = zero.(mag.(remainder.(vTMN)))
    qrprecondition!(leftTMN, rightTMN, linTN, rems, scaleV, vTMN)
    return leftTMN, rightTMN
end

"""
    qrprecondition!(leftTMN::Vector{TaylorModelN{N,T,S}}, rightTMN::Vector{TaylorModelN{N,T,S}},
        linTN::Matrix{T}, rems::Vector{Interval{S}}, scaleV::Vector{T},
        vTMN::Vector{TaylorModelN{N,T,S}}) where {N,T,S}

In-place implementation of qrprecondition.
"""
function qrprecondition!(
        leftTMN::Vector{TaylorModelN{N,T,S}}, rightTMN::Vector{TaylorModelN{N,T,S}},
        linTN::Matrix{T}, rems::Vector{Interval{S}}, scaleV::Vector{T},
        vTMN::Vector{TaylorModelN{N,T,S}}) where {N,T,S}
    # Initialize TMN polynomials
    for ind in eachindex(vTMN)
        for ordQ in eachindex(vTMN[ind].pol.coeffs)
            for hp in eachindex(vTMN[ind].pol.coeffs[ordQ].coeffs)
                leftTMN[ind].pol.coeffs[ordQ].coeffs[hp] =
                    zero(leftTMN[ind].pol.coeffs[ordQ].coeffs[hp])
                rightTMN[ind].pol.coeffs[ordQ].coeffs[hp] =
                    zero(rightTMN[ind].pol.coeffs[ordQ].coeffs[hp])
            end
        end
    end
    # Linear matrix: TS.jacobian(linear_polynomial(vTMN))
    for jind in eachindex(vTMN)
        for ind in eachindex(vTMN)
            linTN[ind, jind] = vTMN[ind].pol.coeffs[2].coeffs[jind]
        end
    end
    # QR factorization of linTN
    qqr = qr(linTN)
    linTN .= qqr.Q * I # Reuse memory
    # Transform remainders
    mul!(rems, linTN', remainder.(vTMN))
    # Get scaling vector, so range of rightTMN is contained in [-1,1]
    for ind in eachindex(vTMN)
        # Constant term
        leftTMN[ind].pol.coeffs[1].coeffs[1] = vTMN[ind].pol.coeffs[1].coeffs[1]
        # Linear corrections
        for hp in eachindex(vTMN[ind].pol.coeffs[2])
            leftTMN[ind].pol.coeffs[2].coeffs[hp] = linTN[ind, hp]
            rightTMN[ind].pol.coeffs[2].coeffs[hp] = qqr.R[ind, hp]
        end
        # Higher order terms (rightTMN)
        for jind in eachindex(vTMN)
            for ordQ in eachindex(vTMN[ind].pol.coeffs)
                ordQ < 3 && continue
                for hp in eachindex(vTMN[ind].pol.coeffs[ordQ].coeffs)
                    rightTMN[ind].pol.coeffs[ordQ].coeffs[hp] +=
                        linTN[jind, ind] * vTMN[jind].pol.coeffs[ordQ].coeffs[hp]
                end
            end
        end
        # Scaling vector
        scaleV[ind] = mag(evaluate(rightTMN[ind].pol, domain(vTMN[ind])) + rems[ind])
    end
    # Exploit the scaled vars
    for ind in eachindex(vTMN)
        # Constant terms remains unchanged
        # Linear corrections
        for hp in eachindex(vTMN[ind].pol.coeffs[2])
            leftTMN[ind].pol.coeffs[2].coeffs[hp] = linTN[ind, hp] * scaleV[hp]
            rightTMN[ind].pol.coeffs[2].coeffs[hp] *= inv(scaleV[ind])
        end
        # Higher order terms (rightTMN)
        for ordQ in eachindex(vTMN[ind].pol.coeffs)
            ordQ < 3 && continue
            for hp in eachindex(vTMN[ind].pol.coeffs[2])
                rightTMN[ind].pol.coeffs[ordQ].coeffs[hp] *= inv(scaleV[ind])
            end
        end
    end
    # Output
    for ind in eachindex(vTMN)
        for ordQ in eachindex(vTMN[ind].pol.coeffs)
            for hp in eachindex(vTMN[ind].pol.coeffs[ordQ].coeffs)
                leftTMN[ind].pol.coeffs[ordQ].coeffs[hp] =
                    leftTMN[ind].pol.coeffs[ordQ].coeffs[hp]
                rightTMN[ind].pol.coeffs[ordQ].coeffs[hp] =
                    rightTMN[ind].pol.coeffs[ordQ].coeffs[hp]
            end
        end
        leftTMN[ind].rem = zero(rems[ind])
        rightTMN[ind].rem = rems[ind]
        leftTMN[ind].x0 = vTMN[ind].x0
        rightTMN[ind].x0 = vTMN[ind].x0
        leftTMN[ind].dom = vTMN[ind].dom
        rightTMN[ind].dom = vTMN[ind].dom
    end
    return nothing
end
