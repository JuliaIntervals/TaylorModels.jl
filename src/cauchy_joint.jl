# Joint Cauchy accumulation for TaylorModel1{TaylorModelN{N,Interval{S},S},S}
#
# The product (or square) of two `TaylorModel1` whose coefficients are `TaylorModelN`
# computes, for each power τ^k, the sum Σ_{i+j=k} a_i ⊗ b_j of `TaylorModelN` products. The
# generic code (`mul!` of `Taylor1`) does one TaylorModelN product per pair (i,j), each with its
# own bound of the neglected polynomial terms (degrees above `orderQ`) and of the cross terms
# with the remainders, i.e. each pair evaluates the neglected polynomials and the ranges of
# its factors over the box.
#
# Here the whole coefficient is bounded jointly (all the truncation and rounding errors of
# the coefficient are collected in one remainder):
#
#   * the retained polynomial is Σ_{i+j=k} trunc_Q(p_i q_j) (products of homogeneous pieces,
#     `_hp_muladd!`, accumulated in place);
#   * the neglected terms of degree Q+1, …, Q+ne (`ne = _TMN_TAIL_EXACT[]`) are accumulated
#     over all the pairs in one homogeneous polynomial per degree, evaluated once over the
#     box (a bound at most as large as the sum of the bounds of each pair); the higher degrees
#     use the ranges `[p_i,u][q_j,v]` of the homogeneous pieces, computed once per
#     coefficient (not once per pair);
#   * the cross terms with the remainders, p_i σ_j + ρ_i q_j + ρ_i σ_j, use the ranges of
#     the polynomials over the box (`Δp_i`, computed once per coefficient), and, as in
#     `remainder_product`, the intersection of the three equivalent groupings.
#
# The result is a rigorous enclosure, at most as wide as the one of the generic method. It is
# used if `_TM1_JOINT[]` (default `true`); the generic methods are kept (used when it is
# `false`). The kernel for the homogeneous products is selected by `_TMN_MIDRAD[]`.

const _TM1_JOINT = Ref(true)

# Accumulate each coefficient in midpoint-radius form (`_cauchy_joint_acc`, below) when
# `_TMN_MIDRAD[]` and the extension provides it; `false` uses the pair-by-pair accumulation
# (`_cauchy_joint_pairs`). Mainly for comparisons.
const _TM1_ACC = Ref(true)

# In `_cauchy_joint_acc`, with the σ-box `[-1,1]^N`, the ranges of the homogeneous pieces and of the
# neglected degrees are computed from the midpoint-radius data (`TS.midrad_range`), without interval
# arithmetic; `false` uses the interval evaluation (`_range_matrix`, `T[d](box)`). Mainly for comparisons.
const _TM1_UNIT = Ref(true)

# retained polynomial and neglected degrees Q+1, …, Q+ne of p*q accumulated in place
@inline function _pair_poly!(ck, T, p, q, Q::Int, ne::Int)
    @inbounds for u in 0:Q, v in 0:Q
        d = u + v
        if d <= Q
            _hp_muladd!(ck.pol.coeffs[d+1], p.coeffs[u+1], q.coeffs[v+1])
        elseif d <= Q + ne
            _hp_muladd!(T[d-Q], p.coeffs[u+1], q.coeffs[v+1])
        end
    end
    return nothing
end

# ranges of the homogeneous pieces of the coefficients, R[u+1, i+1] = p_{i,u}(box): one matrix per
# operand (not one vector per coefficient), and the ranges of the polynomials Δ_i = Σ_u R[u+1, i+1]
function _range_matrix(pa::Taylor1{TaylorModelN{N,Interval{S},S}}, box, n::Int, Q::Int) where {N,S}
    R = Matrix{Interval{S}}(undef, Q + 1, n + 1)
    @inbounds for i in 0:n
        pol = pa[i].pol
        for u in 0:Q
            R[u+1, i+1] = pol.coeffs[u+1](box)
        end
    end
    return R
end

# the same from the midpoint-radius data of the coefficients (box [-1,1]^N)
function _range_matrix(M, n::Int, Q::Int, ::Type{S}) where {S}
    R = Matrix{Interval{S}}(undef, Q + 1, n + 1)
    @inbounds for i in 0:n, u in 0:Q
        R[u+1, i+1] = TS.midrad_range(M[i+1], u)
    end
    return R
end

function _col_sums(R::Matrix{Interval{S}}) where {S}
    Δ = Vector{Interval{S}}(undef, size(R, 2))
    @inbounds for j in axes(R, 2)
        s = R[1, j]
        for u in 2:size(R, 1)
            s += R[u, j]
        end
        Δ[j] = s
    end
    return Δ
end

# Σ_{u+v > Q+ne} RA[u, i] RB[v, j]  (ranges of the homogeneous pieces; i, j: τ orders)
@inline function _pair_ranges(RA::Matrix{Interval{S}}, i::Int, RB::Matrix{Interval{S}}, j::Int,
        Q::Int, ne::Int) where {S}
    s = zero(Interval{S})
    @inbounds for u in 0:Q, v in 0:Q
        u + v > Q + ne && (s += RA[u+1, i+1] * RB[v+1, j+1])
    end
    return s
end

# Σ_{u+v > Q+ne} RA[u, i] RA[v, i]
@inline function _pair_ranges_sq(RA::Matrix{Interval{S}}, i::Int, Q::Int, ne::Int) where {S}
    s = zero(Interval{S})
    two = interval(S, 2)
    @inbounds for k in Q+ne+1:2*Q
        c = zero(Interval{S})
        for u in max(0, k-Q):((k-1) >> 1)            # u < k-u
            c += RA[u+1, i+1] * RA[k-u+1, i+1]
        end
        c = two * c
        iseven(k) && (c += RA[(k >> 1)+1, i+1]^2)
        s += c
    end
    return s
end

# cross terms of (p+ρ)(q+σ) other than the neglected polynomial: p σ + ρ q + ρ σ
@inline function _pair_cross(Δp, ρ, Δq, σ)
    c0 = Δq * ρ + Δp * σ + ρ * σ
    c1 = Δq * ρ + (Δp + ρ) * σ
    c2 = Δp * σ + (Δq + σ) * ρ
    return _intersect_reminders((c1, c2, c0))
end

# cross terms of (p+ρ)²
@inline function _pair_cross_sq(Δp, ρ)
    two = interval(2)
    return _intersect_reminders(((two * Δp + ρ) * ρ, two * Δp * ρ + ρ^2))
end

@inline function _double!(p::TaylorN{Interval{S}}) where {S}
    two = interval(S, 2)
    @inbounds for hp in p.coeffs, idx in eachindex(hp.coeffs)
        hp.coeffs[idx] = two * hp.coeffs[idx]
    end
    return nothing
end
@inline function _double!(hp::HomogeneousPolynomial{Interval{S}}) where {S}
    two = interval(S, 2)
    @inbounds for idx in eachindex(hp.coeffs)
        hp.coeffs[idx] = two * hp.coeffs[idx]
    end
    return nothing
end

# Taylor1{TaylorModelN} of the product a*b (or of the square, if `issquare`); the products of
# homogeneous pieces are accumulated pair by pair (`_hp_muladd!`)
function _cauchy_joint_pairs(a::TaylorModel1{TaylorModelN{N,Interval{S},S},S},
        b::TaylorModel1{TaylorModelN{N,Interval{S},S},S}, issquare::Bool) where {N,S}
    pa, pb = a.pol, b.pol
    n = TS.order(pa)
    a0 = pa[0]
    Q = TS.order(a0.pol)
    @assert 2*Q ≤ TS.order(space(a0))
    z = zero(Interval{S})
    two = interval(S, 2)
    box = centered_dom(a0)
    ne = max(0, min(_TMN_TAIL_EXACT[], Q))              # exact degrees Q+1, …, Q+ne
    # ranges of the homogeneous pieces and of the polynomials, remainders (once per coefficient)
    RA = _range_matrix(pa, box, n, Q)
    ΔA = _col_sums(RA)
    ρA = [remainder(pa[i]) for i in 0:n]
    if issquare
        RB, ΔB, ρB = RA, ΔA, ρA
    else
        RB = _range_matrix(pb, box, n, Q)
        ΔB = _col_sums(RB)
        ρB = [remainder(pb[i]) for i in 0:n]
    end
    T = [HomogeneousPolynomial(space(a0), z, Q + d) for d in 1:ne]   # accumulators
    cs = Vector{typeof(a0)}(undef, n + 1)
    for k in 0:n
        ck = zero(a0)
        for d in 1:ne
            fill!(T[d].coeffs, z)
        end
        Rk = z
        if issquare
            Roff = z
            @inbounds for i in 0:((k-1) >> 1)             # i < j = k-i
                j = k - i
                _pair_poly!(ck, T, pa[i].pol, pa[j].pol, Q, ne)
                ne < Q && (Roff += _pair_ranges(RA, i, RA, j, Q, ne))
                Roff += _pair_cross(ΔA[i+1], ρA[i+1], ΔA[j+1], ρA[j+1])
            end
            _double!(ck.pol)
            for d in 1:ne
                _double!(T[d])
            end
            Rk = two * Roff
            if iseven(k)                                   # diagonal term
                i = k >> 1
                _pair_poly!(ck, T, pa[i].pol, pa[i].pol, Q, ne)
                ne < Q && (Rk += _pair_ranges_sq(RA, i, Q, ne))
                Rk += _pair_cross_sq(ΔA[i+1], ρA[i+1])
            end
        else
            @inbounds for i in 0:k
                j = k - i
                _pair_poly!(ck, T, pa[i].pol, pb[j].pol, Q, ne)
                ne < Q && (Rk += _pair_ranges(RA, i, RB, j, Q, ne))
                Rk += _pair_cross(ΔA[i+1], ρA[i+1], ΔB[j+1], ρB[j+1])
            end
        end
        for d in 1:ne
            Rk += T[d](box)
        end
        ck.rem = Rk
        cs[k+1] = ck
    end
    return Taylor1(cs, n)
end

# Same, accumulating each τ^k coefficient in midpoint-radius form (TaylorSeriesIAExt:
# `TS.midrad_poly`, `TS.mul_midrad_acc!`, `TS.midrad_finalize!`): the polynomials are converted
# to midpoints and radii once per coefficient, the products of all the pairs and degrees are
# accumulated in hardware arithmetic (with error-free transformations) and converted to an
# `Interval` once per coefficient. Returns `nothing` if the midpoint-radius form is not
# usable (decorations, non finite values): the caller then uses `_cauchy_joint_pairs`.
function _cauchy_joint_acc(a::TaylorModel1{TaylorModelN{N,Interval{S},S},S},
        b::TaylorModel1{TaylorModelN{N,Interval{S},S},S}, issquare::Bool) where {N,S}
    pa, pb = a.pol, b.pol
    n = TS.order(pa)
    a0 = pa[0]
    Q = TS.order(a0.pol)
    @assert 2*Q ≤ TS.order(space(a0))
    z = zero(Interval{S})
    two = interval(S, 2)
    twof = S(2)
    box = centered_dom(a0)
    ne = max(0, min(_TMN_TAIL_EXACT[], Q))
    MA = [TS.midrad_poly(pa[i].pol) for i in 0:n]
    all(TS.midrad_usable, MA) || return nothing
    if issquare
        MB = MA
    else
        MB = [TS.midrad_poly(pb[i].pol) for i in 0:n]
        all(TS.midrad_usable, MB) || return nothing
    end
    # ranges over the σ-box: from the midpoint-radius data if the box is [-1,1]^N
    unit = _TM1_UNIT[] && all(x -> isequal_interval(x, interval(-one(S), one(S))), box)
    RA = unit ? _range_matrix(MA, n, Q, S) : _range_matrix(pa, box, n, Q)
    ΔA = _col_sums(RA)
    ρA = [remainder(pa[i]) for i in 0:n]
    if issquare
        RB, ΔB, ρB = RA, ΔA, ρA
    else
        RB = unit ? _range_matrix(MB, n, Q, S) : _range_matrix(pb, box, n, Q)
        ΔB = _col_sums(RB)
        ρB = [remainder(pb[i]) for i in 0:n]
    end
    acc = TS.midrad_acc(a0.pol, Q + ne)
    T = [HomogeneousPolynomial(space(a0), z, Q + d) for d in 1:(unit ? 0 : ne)]   # tails by interval evaluation
    cs = Vector{typeof(a0)}(undef, n + 1)
    for k in 0:n
        ck = zero(a0)
        TS.midrad_reset!(acc)
        Rk = z
        if issquare
            Roff = z
            @inbounds for i in 0:((k-1) >> 1)             # i < j = k-i: counted twice
                j = k - i
                TS.mul_midrad_acc!(acc, MA[i+1], MA[j+1], twof)
                ne < Q && (Roff += _pair_ranges(RA, i, RA, j, Q, ne))
                Roff += _pair_cross(ΔA[i+1], ρA[i+1], ΔA[j+1], ρA[j+1])
            end
            Rk = two * Roff
            if iseven(k)                                   # diagonal term
                i = k >> 1
                TS.mul_midrad_acc!(acc, MA[i+1], MA[i+1])
                ne < Q && (Rk += _pair_ranges_sq(RA, i, Q, ne))
                Rk += _pair_cross_sq(ΔA[i+1], ρA[i+1])
            end
        else
            @inbounds for i in 0:k
                j = k - i
                TS.mul_midrad_acc!(acc, MA[i+1], MB[j+1])
                ne < Q && (Rk += _pair_ranges(RA, i, RB, j, Q, ne))
                Rk += _pair_cross(ΔA[i+1], ρA[i+1], ΔB[j+1], ρB[j+1])
            end
        end
        # retained polynomial and neglected degrees Q+1, …, Q+ne: one interval per monomial
        ok = true
        for d in 0:Q
            ok &= TS.midrad_finalize!(ck.pol.coeffs[d+1], acc, d)
        end
        for d in 1:ne
            if unit                                   # range of the accumulated degree Q+d
                Rk += TS.midrad_range(acc, Q + d)
            else
                fill!(T[d].coeffs, z)
                ok &= TS.midrad_finalize!(T[d], acc, Q + d)
                Rk += T[d](box)
            end
        end
        ok || return nothing
        ck.rem = Rk
        cs[k+1] = ck
    end
    return Taylor1(cs, n)
end

# Taylor1{TaylorModelN} of the product a*b (or of the square): joint accumulation with the
# midpoint-radius form if `_TM1_ACC[]`, `_TMN_MIDRAD[]` and it is usable, else pair by pair
function _cauchy_joint(a::TaylorModel1{TaylorModelN{N,Interval{S},S},S},
        b::TaylorModel1{TaylorModelN{N,Interval{S},S},S}, issquare::Bool) where {N,S}
    if _TM1_ACC[] && _TMN_MIDRAD[] && isdefined(TS, :midrad_poly)
        res = _cauchy_joint_acc(a, b, issquare)
        res === nothing || return res
    end
    return _cauchy_joint_pairs(a, b, issquare)
end

# Multiplication
function *(a::TaylorModel1{TaylorModelN{N,Interval{S},S},S},
        b::TaylorModel1{TaylorModelN{N,Interval{S},S},S}) where {N,S}
    @assert all(isequal_interval.(tmdata(a), tmdata(b)))
    if _TM1_JOINT[] && TS.order(a.pol) == TS.order(b.pol)
        res = _cauchy_joint(a, b, false)
    else
        res = a.pol * b.pol
    end
    Δ = remainder_product(a, b, TS.order(res))
    return unsafe_TaylorModel1(res, Δ, expansion_point(a), domain(a))
end

# Square
function TS.square(a::TaylorModel1{TaylorModelN{N,Interval{S},S},S}) where {N,S}
    res = _TM1_JOINT[] ? _cauchy_joint(a, a, true) : TS.square(a.pol)
    Δ = remainder_square(a, TS.order(res))
    return unsafe_TaylorModel1(res, Δ, expansion_point(a), domain(a))
end