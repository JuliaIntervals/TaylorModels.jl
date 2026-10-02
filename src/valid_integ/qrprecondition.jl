# To be included inside `TaylorModels.ValidatedInteg`.
#
# Contents
#   _right_in_box(rightTMN)            : check ρ(B) ⊆ B (remainder included)
#   qrprecondition_rig!(...)           : rigorous QR preconditioning
#   reconstruct_rho(sol)               : rebuild ρ_k from a TMSol3
#
# Types: TaylorModelN{N,T,U}, with coefficient type T ∈ {U, Interval{U}} and
# U <: AbstractFloat the numeric type of remainders, domains and expansion points.
# All scalar arithmetic (scales, linear algebra, eps) is done in U; all
# accumulation in Interval{U}.
#
# Assumptions (asserted): every TMN has expansion point 0 and domain B = [-1,1]^N,
# so that |ξ^α| ≤ 1 for every monomial.


# --- coefficient-type helpers -------------------------------------------------
_midU(x::Interval) = mid(x)
_midU(x::AbstractFloat) = x
_asT(::Type{T}, x::T) where {T<:AbstractFloat} = x
_asT(::Type{Interval{U}}, x::U) where {U<:AbstractFloat} = interval(x)
_coef_equal(x::AbstractFloat, y::AbstractFloat) = x === y
_coef_equal(x::Interval, y::Interval) = isequal_interval(x, y)
function _pol_equal(p, q)
    for k in eachindex(p.coeffs), h in eachindex(p.coeffs[k].coeffs)
        _coef_equal(p.coeffs[k].coeffs[h], q.coeffs[k].coeffs[h]) || return false
    end
    return true
end


"""
    _right_in_box(rightTMN) -> Bool

`true` iff the range enclosure `monotonicity_bounder` of each `rightTMN[i]`
is contained in `[-1,1]`. Must use the same bounder as `qrprecondition_rig!`.
"""
function _right_in_box(rightTMN::Vector{TaylorModelN{N,T,S}}) where {N,T,S}
    B = domain(rightTMN[1])
    @inbounds for i in eachindex(rightTMN)
        issubset_interval(monotonicity_bounder(rightTMN[i], B), B[i]) || return false
    end
    return true
end


# Interval matrix product with explicit loops (outward rounding at every op).
function _imatmul(A::AbstractMatrix{Interval{S}},
        B::AbstractMatrix{Interval{S}}) where {S}
    m, n = size(A)
    n2, p = size(B)
    @assert n == n2
    C = Matrix{Interval{S}}(undef, m, p)
    @inbounds for i in 1:m, j in 1:p
        acc = zero(Interval{S})
        for k in 1:n
            acc += A[i,k] * B[k,j]
        end
        C[i,j] = acc
    end
    return C
end


"""
    _verified_inv(M, P)

Rigorous enclosure of `inv(M)` (M a float or interval matrix; for an interval
matrix, the enclosure holds for every point matrix in `M`) from an approximate
inverse `P`.
With `E = I - P*M` (interval) and `δ = ‖E‖∞ < 1`:
`inv(M) = (I-E)^{-1} P ∈ (I + E + [-r,r]) P`,  `r = δ²/(1-δ)`,
since each entry of `Σ_{m≥2} E^m` is bounded by `‖E‖∞^m`.
Returns `nothing` if `δ ≥ 1`.

Refs: A. Neumaier, Interval Methods for Systems of Equations, Cambridge Univ.
Press (1990); S. M. Rump, Acta Numerica 19, 287-449 (2010).
"""
_verified_inv(M::AbstractMatrix{T}, P::AbstractMatrix{T}) where {T<:AbstractFloat} =
    _verified_inv(interval.(M), P)
function _verified_inv(M::AbstractMatrix{Interval{T}}, P::AbstractMatrix{T}) where {T<:AbstractFloat}
    N = size(M, 1)
    oI, zI = interval(one(T)), zero(Interval{T})
    E = _imatmul(interval.(P), M)
    @inbounds for i in 1:N, j in 1:N
        E[i,j] = (i == j ? oI : zI) - E[i,j]
    end
    # δ: rigorous upper bound of ‖E‖∞
    δ = zero(T)
    @inbounds for i in 1:N
        rowI = zI
        for j in 1:N
            rowI += interval(mag(E[i,j]))
        end
        δ = max(δ, sup(rowI))
    end
    δ < one(T) || return nothing
    δI = interval(δ)
    r = sup(δI * δI / (oI - δI))
    G = similar(E)
    @inbounds for i in 1:N, j in 1:N
        G[i,j] = (i == j ? oI : zI) + E[i,j] + interval(-r, r)
    end
    return _imatmul(G, interval.(P))
end


"""
    _imatvec_tmn!(out, G, vTMN, c)

`out[i] ⊇ Σ_j G[i,j] * (vTMN[j] - c[j])` on `B`, with `c::Vector{U}`.
Coefficients are accumulated in `Interval{U}` and stored with `_split`; any
residual, times a monomial bounded by 1 on B, goes into the remainder.
"""
function _imatvec_tmn!(out::Vector{TaylorModelN{N,T,U}}, G::AbstractMatrix{Interval{U}},
        vTMN::Vector{TaylorModelN{N,T,U}}, c::AbstractVector{U}) where {N,U,T}
    zI = zero(Interval{U})
    uI = interval(-one(U), one(U))
    @inbounds for i in eachindex(out)
        remi = zI
        for ordQ in eachindex(out[i].pol.coeffs)
            cf = out[i].pol.coeffs[ordQ].coeffs
            for h in eachindex(cf)
                acc = zI
                for j in eachindex(vTMN)
                    vjh = interval(vTMN[j].pol.coeffs[ordQ].coeffs[h])
                    ordQ == 1 && (vjh -= interval(c[j]))       # constant of v - c
                    acc += G[i,j] * vjh
                end
                cf[h], res = TM._split(T, acc)
                remi += ordQ == 1 ? res : res * uI
            end
        end
        for j in eachindex(vTMN)
            remi += G[i,j] * remainder(vTMN[j])
        end
        out[i].rem = remi
        out[i].x0  = vTMN[i].x0
        out[i].dom = vTMN[i].dom
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
    qrprecondition!(leftTMN, rightTMN, linTN, rems, scaleV, vTMN;
                        maxiter=4, verbose=false) -> Bool

Rigorous QR preconditioning of `vTMN::Vector{TaylorModelN{N,T,U}}`,
`T ∈ {U, Interval{U}}`.

On success (`true`):
  * `leftTMN[i] = c_i + Σ_j M_L[i,j] σ_j`, zero remainder; `c = mid(v(0))`, or,
    with `recenter = true` (default), `c` shifted so that the computed ranges of
    `Q^{-1}(v - c)` are centred (tighter scales; ρ then has a small constant term);
    `M_L = fl(Q*S)` (stored as thin coefficients of type `T`);
  * `rightTMN = ρ ⊇ M_L^{-1}(v - c)`, with `ρ(B) ⊆ B` verified;
hence `v(ξ) ∈ L(ρ(ξ))` and `v(B) ⊆ L(B)`. On exit, `linTN` holds `M_L`.
`Q` comes from the QR factorization of `mid` of the linear coefficients of `v`;
it only needs to be a good frame, rigor comes from the verified inverses.
Returns `false` if a verified inverse cannot be obtained or `ρ(B) ⊆ B`
cannot be verified in `maxiter` inflations of `scaleV`.

# References
- R. J. Lohner, "Enclosing the solutions of ordinary initial and boundary value
  problems", in E. Kaucher, U. Kulisch, Ch. Ullrich (eds.), Computer Arithmetic:
  Scientific Computation and Programming Languages, Teubner, Stuttgart (1987),
  pp. 255-286. [QR method against the wrapping effect]
- K. Makino, M. Berz, "Suppression of the wrapping effect by Taylor model-based
  verified integrators: long-term stabilization by preconditioning",
  Int. J. Differ. Equ. Appl. 10(4), 353-384 (2005). [left/right TM factorization,
  QR preconditioning]
- M. Neher, K. R. Jackson, N. S. Nedialkov, "On Taylor model based integration
  of ODEs", SIAM J. Numer. Anal. 45(1), 236-262 (2007).
  https://doi.org/10.1137/050638448
- F. Bünger, "Preconditioning of Taylor models, implementation and test cases",
  Nonlinear Theory and Its Applications, IEICE 12(1), 2-40 (2021).
  https://doi.org/10.1587/nolta.12.2 [explicit rigorous implementation, INTLAB]
- S. M. Rump, "Verification methods: Rigorous results using floating-point
  arithmetic", Acta Numerica 19, 287-449 (2010). [verified inverse]
"""
function qrprecondition_rig!(
        leftTMN::Vector{TaylorModelN{N,T,U}}, rightTMN::Vector{TaylorModelN{N,T,U}},
        linTN::Matrix{T}, rems::Vector{Interval{U}}, scaleV::Vector{U},
        vTMN::Vector{TaylorModelN{N,T,U}};
        maxiter::Int = 4, verbose::Bool = false, recenter::Bool = true
        ) where {N, U<:AbstractFloat, T<:Union{U,Interval{U}}}
    B = domain(vTMN[1])
    @assert all(isequal_interval.(B, symmetric_box(N, U)))
    @assert all(iszero.(inf.(vTMN[1].x0))) && all(iszero.(sup.(vTMN[1].x0)))

    # centre c and (mid of) Jacobian A; float QR
    c = Vector{U}(undef, N)
    @inbounds for i in 1:N
        c[i] = _midU(vTMN[i].pol.coeffs[1].coeffs[1])
        for j in 1:N
            linTN[i,j] = vTMN[i].pol.coeffs[2].coeffs[j]     # A, type T
        end
    end
    Amid = _midU.(linTN)                    # QR frame from mid(A) (Q need not be exact)
    if !all(isfinite, Amid) || !all(isfinite, c)
        verbose && @warn "qrprecondition_rig!: non-finite centre or Jacobian" c Amid
        return false
    end
    Q = Matrix(qr(Amid).Q)
    Qt = permutedims(Q)

    # (1) Rigorous scale: s_i ≥ sup_B |(Q^{-1}(v-c))_i|
    Qinv = _verified_inv(Q, Qt)
    if Qinv === nothing
        verbose && @warn "qrprecondition_rig!: stage 1, verified inv(Q) failed" Amid Q
        return false
    end
    _imatvec_tmn!(rightTMN, Qinv, vTMN, c)
    if recenter
        # Shift c ← c + Q m, m_i = mid of the range of (Q^{-1}(v-c))_i, so that the
        # ranges become (nearly) symmetric. Any float c is sound: stage 2 verifies
        # ρ = M_L^{-1}(v - c) for the c actually used.
        @inbounds for i in 1:N
            scaleV[i] = mid(monotonicity_bounder(rightTMN[i], B))     # m, stored temporarily
        end
        if all(isfinite, scaleV)
            @inbounds for i in 1:N
                c[i] += sum(Q[i,j] * scaleV[j] for j in 1:N)
            end
            _imatvec_tmn!(rightTMN, Qinv, vTMN, c)
        end
    end
    @inbounds for i in 1:N
        scaleV[i] = mag(monotonicity_bounder(rightTMN[i], B))
    end
    verbose && println("stage 1: cond(mid A) = ", cond(Amid), "  scaleV = ", repr(scaleV))
    smax = maximum(scaleV)
    if !isfinite(smax)
        verbose && @warn "qrprecondition_rig!: stage 1, non-finite scale" scaleV
        return false
    elseif iszero(smax)                 # v is constant: any s > 0 is sound
        scaleV .= one(U)
    else                                # avoid s_i = 0 (degenerate directions)
        @. scaleV = max(scaleV, eps(U) * smax)
    end

    # (2) Build L, verify ρ = M_L^{-1}(v-c) ⊆ B; inflate s if needed.
    #     M_L = fl(Q*S) = Q̂ S with Q̂[i,j] ∈ fl(Q[i,j] s_j) / s_j (interval division),
    #     so M_L^{-1} ∈ S^{-1} [Q̂^{-1}]; Q̂ ≈ Q is well conditioned for any S.
    Qhat = Matrix{Interval{U}}(undef, N, N)
    G = Matrix{Interval{U}}(undef, N, N)
    for it in 1:maxiter
        @inbounds for i in 1:N
            for ordQ in eachindex(leftTMN[i].pol.coeffs)
                cf = leftTMN[i].pol.coeffs[ordQ].coeffs
                for h in eachindex(cf)
                    cf[h] = zero(T)
                end
            end
            leftTMN[i].pol.coeffs[1].coeffs[1] = _asT(T, c[i])
            for j in 1:N
                mL = Q[i,j] * scaleV[j]               # defines M_L (rounding irrelevant)
                linTN[i,j] = _asT(T, mL)
                leftTMN[i].pol.coeffs[2].coeffs[j] = linTN[i,j]
                Qhat[i,j] = interval(mL) / interval(scaleV[j])
            end
            leftTMN[i].rem = zero(Interval{U})
            leftTMN[i].x0  = vTMN[i].x0
            leftTMN[i].dom = vTMN[i].dom
        end
        Qhinv = _verified_inv(Qhat, Qt)
        if Qhinv === nothing
            verbose && @warn "qrprecondition_rig!: stage 2, verified inv(Q̂) failed" it scaleV
            return false
        end
        @inbounds for i in 1:N, j in 1:N
            G[i,j] = Qhinv[i,j] / interval(scaleV[i])  # S^{-1} Q̂^{-1}
        end
        _imatvec_tmn!(rightTMN, G, vTMN, c)
        ok = true
        verbose && println("stage 2, it = $it  scaleV = ", repr(scaleV))
        @inbounds for i in 1:N
            rems[i] = remainder(rightTMN[i])
            r = monotonicity_bounder(rightTMN[i], B)
            m = mag(r)
            verbose && println("  ρ[$i]: range = ", repr((inf(r), sup(r))),
                               "  rem = ", repr((inf(rems[i]), sup(rems[i]))))
            if !(m <= one(U))                          # also catches NaN
                ok = false
                # Margin beats the O(ε s_j/s_i) rounding noise of ρ_i, which
                # changes between iterations (Q̂ = fl(QS)/S is re-rounded).
                scaleV[i] *= one(U) + 4 * (m - one(U)) + U(16)^it * eps(U)
            end
        end
        ok && return true
    end
    verbose && @warn "qrprecondition_rig!: ρ(B) ⊆ B not verified after maxiter" maxiter scaleV
    return false
end


"""
    reconstruct_rho(sol::TMSol3) -> Matrix{TaylorModelN}

Returns `ρ` (size `N × length(sol)`) with `ρ[:,k]` the right factor of step `k`:
`Φ_k(h_k, σ) ⊆ L_{k+1}(ρ[:,k](σ))`, `σ ∈ B`, `ρ[:,k](B) ⊆ B`. Column 1 is the identity.

It recomputes `v_k = Φ_k(h_k, ·)` with `__evaluate_rig!` (`h_k` exact, from
`domain(sol, k)`) and calls `qrprecondition_rig!` with its default keywords —
the same calls, with the same keywords, as `validated_integ3` (consistency:
if the integrator ever passes non-default keywords, they must be passed here
too). It asserts `ρ[:,k](B) ⊆ B` and that `L_{k+1}` equals, bit-wise, the τ^0
coefficient stored in `sol[k+1]`.

Composition (dependence on the initial variables): `σ_2 = ξ`,
`σ_{k+1} = ρ[:,k](σ_k)`, `x(t_k + τ) ∈ Φ_k(τ, σ_k(ξ))`. With re-centring,
`ρ[:,k](0) ≠ 0`: the constant term of `L_k` is the centre of the enclosing box,
not the trajectory from `mid(X0)`, which is `Φ_k(τ, σ_k(0))`.
"""
function reconstruct_rho(sol::TMSol3{N,T,U}) where {N,T,U}
    xTM = get_xTM(sol)
    nk = size(xTM, 2)
    proto = xTM[1, 2][0]
    vTMN  = [deepcopy(proto) for _ in 1:N]
    left  = [deepcopy(proto) for _ in 1:N]
    right = [deepcopy(proto) for _ in 1:N]
    CT = typeof(proto.pol.coeffs[1].coeffs[1])     # coefficient type (U or Interval{U})
    linTN  = zeros(CT, N, N)
    rems   = fill(zero(Interval{U}), N)
    scaleV = zeros(U, N)
    aux    = TM._aux_horner(proto)                    # workspace for __evaluate_rig!

    ρ = Matrix{typeof(proto)}(undef, N, nk)
    for i in 1:N                                   # identity in column 1
        ρ[i,1] = deepcopy(proto)
        for ordQ in eachindex(ρ[i,1].pol.coeffs), h in eachindex(ρ[i,1].pol.coeffs[ordQ].coeffs)
            ρ[i,1].pol.coeffs[ordQ].coeffs[h] = zero(CT)
        end
        ρ[i,1].pol.coeffs[2].coeffs[i] = one(CT)
        ρ[i,1].rem = zero(Interval{U})
    end

    for k in 2:nk
        dom = domain(xTM[1,k])                     # sign*interval(0, sign*δt)
        h = sup(dom) > 0 ? sup(dom) : inf(dom)     # exact δt
        for i in 1:N
            TM.__evaluate_rig!(vTMN[i], xTM[i,k], h, aux)
        end
        qrprecondition_rig!(left, right, linTN, rems, scaleV, vTMN) ||
            error("reconstruct_rho: qrprecondition_rig! failed at step $k")
        @assert _right_in_box(right) "reconstruct_rho: ρ(B) ⊄ B at step $k"
        if k < nk
            @assert all(_pol_equal(polynomial(xTM[i,k+1][0]), polynomial(left[i])) for i in 1:N) """
                reconstruct_rho: L_{k+1} ≠ τ^0 coefficient of sol[$(k+1)]
                (non-deterministic, or different evaluation/preconditioning calls)"""
        end
        for i in 1:N
            ρ[i,k] = deepcopy(right[i])
        end
    end
    return ρ
end