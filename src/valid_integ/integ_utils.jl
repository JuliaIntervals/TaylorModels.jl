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
    B = symmetric_box(S, space(xTMN[1]))
    @assert all(isequal_interval.(domain.(xTMN), (B,)))
    x0 = zero(B)
    @assert all(expansion_point.(xTMN) .== (x0,))

    # Vector of independent TaylorN variables
    order = TS.order(xTMN[1])
    # X = [TaylorN(T, i, order=order) for i in 1:N]
    X = variables(space(xTMN[1]), order=order)

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
            X::Vector{TaylorN{T}}) where {N,T, S<:IANumTypes}
        postverify = true
        x0 = expansion_point(xTMN[1])
        B = domain(xTMN[1])
        zI = zero(Interval{S})
        oI = one($TT)
        @inbounds for i in eachindex(xTMN)
            pol = polynomial(xTMN[i])
            tmn = TM.unsafe_TaylorModelN(pol(X), zI, x0, B )
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
    orderQ = TS.order(a)
    δ = symmetric_box(T, space(a))
    aux = diam(Δ)/(2N)
    rem = zero(Δ)

    # Linear shift
    lin_shift = mid(Δ) + aux*sum((TaylorN(space(a), i, order=orderQ) for i in 1:N))
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

    return TM.unsafe_TaylorModelN(bpol, rem, expansion_point(a), domain(a))
end


function _abs_rems!(vTMN::Vector{TaylorModelN{N,T,S}}) where {N,T,S}
    for ind in eachindex(vTMN)
        Δ = remainder(vTMN[ind])
        radN = radius(Δ)/N
        # Old remainder of constant and linear parts and remainder
        aI = vTMN[ind].pol.coeffs[1].coeffs[1] +
            sum(vTMN[ind].pol.coeffs[2].coeffs) * vTMN[ind].dom[1] + Δ
        # New remainder of constant and linear parts (without Δ, a priori absorbed)
        bI = vTMN[ind].pol.coeffs[1].coeffs[1] + mid(Δ) +
            sum(vTMN[ind].pol.coeffs[2].coeffs .+ radN) * vTMN[ind].dom[1]
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
    _update_inicond!(x, dx, x1N, vTMN)

In-place update the initial conditions in `x` and `dx` (the latter reset
to zero) from `vTMN`; `x1N` stores the remainder.
"""
function _update_inicond!(x, dx, x1N, vTMN)
    zz = zero(x[1][0][0][1])
    for ind in eachindex(x)
        src1 = vTMN[ind]
        src2 = x1N[ind]
        src2.rem = src1.rem # Store remainder
        # Zero everything
        for ordT in eachindex(src2.pol.coeffs)
            src2_ordT = src2.pol.coeffs[ordT]
            for ordQ in eachindex(src2_ordT.pol.coeffs)
                for h in eachindex(src2_ordT.pol.coeffs[ordQ].coeffs)
                    x[ind].coeffs[ordT].coeffs[ordQ].coeffs[h] = zz
                    # dx[ind].coeffs[ordT].coeffs[ordQ].coeffs[h] = zz
                end
            end
        end
        # Update constant coeff (new initial condition)
        src2_ordT1 = src2.pol.coeffs[1]
        for ordQ in eachindex(src2_ordT1.pol.coeffs)
            for h in eachindex(src2_ordT1.pol.coeffs[ordQ].coeffs)
                x[ind].coeffs[1].coeffs[ordQ].coeffs[h] =
                    src1.pol.coeffs[ordQ].coeffs[h]
                dx[ind].coeffs[1].coeffs[ordQ].coeffs[h] = zz
            end
        end
    end
    return nothing
end


"""
    _update_output!(xTM1v::AbstractVector{TaylorModel1{TaylorModelN{N,T,U}, U}},
        x1N::AbstractVector{TaylorModel1{TaylorModelN{N,T,U}, U}})

In-place update of `xTM1v` with the contents of `x1N`, to avoid `deepcopy`.
Assumes same shape of both entries.
"""
function _update_output!(xTM1v::AbstractVector{<:TaylorModel1}, x1N)
    @inbounds for ind in eachindex(x1N)
        src = x1N[ind]
        tgt = xTM1v[ind]
        tgt.rem = src.rem
        tgt.x0  = src.x0
        tgt.dom = src.dom
        for ordT in eachindex(src.pol.coeffs)
            src_ordT = src.pol.coeffs[ordT]  # a TaylorModelN
            tgt_ordT = tgt.pol.coeffs[ordT]
            tgt_ordT.rem = src_ordT.rem
            tgt_ordT.x0  = src_ordT.x0
            tgt_ordT.dom = src_ordT.dom
            for ordQ in eachindex(src_ordT.pol.coeffs)
                for h in eachindex(src_ordT.pol.coeffs[ordQ].coeffs)
                    tgt_ordT.pol.coeffs[ordQ].coeffs[h] =
                        src_ordT.pol.coeffs[ordQ].coeffs[h]
                end
            end
        end
    end
    return nothing
end


# The following is used to get an exact, floating number representation,
# of the step size δt and of the accumulated time
# Knuth's TwoSum: s + e == a + b exactly
@inline function _two_sum(a::T, b::T) where {T<:AbstractFloat}
    s = a + b
    bb = s - a
    e = (a - (s - bb)) + (b - bb)
    return s, e
end

"""
    _exact_step(t0, δt, tmax, sign_tstep) -> δt′

Step with `t0 + δt′` exact in floating point and `|δt′| ≤ |δt|`, except that
a step reaching `tmax` lands exactly on it. `sign_tstep = ±1` (forward/backward).
"""
function _exact_step(t0::T, δt::T, tmax::T, sign_tstep::Int) where {T<:AbstractFloat}
    t1 = t0 + δt
    sign_tstep*t1 >= sign_tstep*tmax && (t1 = tmax)
    for _ in 1:8
        δt1 = t1 - t0
        s, e = _two_sum(t0, δt1)
        if s == t1 && iszero(e)                          # t0 + δt1 == t1 exactly
            if iszero(δt1)                               # request < ½ ulp(t0): one ulp
                t1 = sign_tstep > 0 ? nextfloat(t0) : prevfloat(t0)
                return t1 - t0
            end
            (t1 == tmax || sign_tstep*δt1 <= sign_tstep*δt) && return δt1
            t1 = sign_tstep > 0 ? prevfloat(t1) : nextfloat(t1)    # ½-ulp overshoot
        else
            # t1 - t0 not representable (t1 not within a factor 2 of t0): take an
            # exact shorter step; later steps reach tmax exactly (Sterbenz).
            if sign_tstep*t0 > 0                         # away from 0: at most doubling
                δt = sign_tstep * min(abs(δt1), abs(t0))
            elseif sign_tstep*tmax >= 0 && abs(t0) <= abs(δt)   # 0 lies on the way
                δt = -t0                                 # stop exactly at 0
            else                                         # toward 0, tmax before 0
                δt = sign_tstep * min(abs(δt1), abs(t0)) / 2
            end
            t1 = t0 + δt
        end
    end
    error("_exact_step: no exact step found (t0 = $t0, δt = $δt)")
end
_exact_step(t0::Real, δt::Real, tmax::Real, ::Int) = δt     # exact types (e.g. Rational)
