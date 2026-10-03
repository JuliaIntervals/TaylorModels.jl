# Tests for validated_integ3 and its building blocks.
#
#   1. _exact_step            representable time steps
#   2. picard_lindelof        no aliasing, exact Picard iterates, argument order
#   3. monotonicity_bounder   range enclosure (soundness, tightness)
#   4. qrprecondition_rig!    (⋆) and v ∈ L∘ρ, exact in BigFloat, T ∈ {Float64, Interval}
#   5. integration            completion, pointwise containment x(t;x0) ∈ Φ_k(τ, σ_k(ξ))
#                             with exact times, tdrift == 0, tightness
#
# Notation as in report_validated_integ3.md.

using TaylorModels, TaylorModels.ValidatedInteg
using Test, Random, StaticArrays, Logging, LinearAlgebra

const VI = TaylorModels.ValidatedInteg
const TI = TaylorIntegration

# ============================================================================
# Shared helpers
# ============================================================================

stepsize(d) = sup(d) > 0 ? sup(d) : inf(d)          # exact h_k from a τ-domain

# Run f(sp′) where sp′ is the stored default space (pass sp′, not sp, to the integrator).
# This is needed for non-autonomous f!, due to a promotion in t (::Taylor1) to Taylor1{TaylorN}
# which happens to use the default space. On exit restores the previous default.
function with_default_space(f, sp)
    old = TS.default_space[]
    TS.set_default_space!(sp; nowarn=true)
    try
        return f(TS.default_space[])
    finally
        TS.set_default_space!(old; nowarn=true)
    end
end

# validated_integ3 with no warnings allowed; `nothing` if it throws
run_integ3(f!, X0, tini, tend, orderQ, orderT, abstol, sp; kwargs...) =
    with_default_space(sp) do splocal
        @test_logs min_level=Logging.Warn validated_integ3(
            f!, X0, tini, tend, orderQ, orderT, abstol, splocal;
            parse_eqs = true, maxsteps = 2000, kwargs...)
    end

# Cache objects for one validated step, as validated_integ3 builds them:
# x1N (= L, constant in τ), dx1N, t1N (= t0 + τ), on τ ∈ [0, δt].
function setup(f!, X0, t0, δt, orderQ, orderT, sp)
    c = VI.init_cache_VI3(f!, t0, X0, 2, orderT, orderQ, sp, nothing; parse_eqs = true)
    δtI = interval(0.0, δt)
    x1N, dx1N, t1N = c.x1N, c.dx1N, c.t1N
    t1N.dom = δtI
    t1N.pol.coeffs[1].pol.coeffs[1].coeffs[1] = t0
    for i in eachindex(x1N)
        x1N[i].dom = δtI
        dx1N[i].dom = δtI
    end
    return c, x1N, dx1N, t1N
end

# --- coefficient access (T = Float64 or Interval{Float64}) -----------------
coeflist(p) = [p.coeffs[k].coeffs[h] for k in eachindex(p.coeffs) for h in eachindex(p.coeffs[k].coeffs)]
# isfin(x::Interval) = isfinite(inf(x)) && isfinite(sup(x)); isfin(x::Real) = isfinite(x)
isfin(x::Interval) = isbounded(x); isfin(x::Real) = isfinite(x)
cpoly(tm1, k) = polynomial(tm1[k])           # TaylorN of the τ^k coefficient
constpoly(x, ref) = x + zero(ref)            # constant TaylorN like ref

# Float TaylorN with coefficient list `vals` and the structure of `tmpl`
function float_poly(vals, tmpl)
    p = zero(vals[1]) * tmpl
    n = 0
    for k in eachindex(p.coeffs)
        for h in eachindex(p.coeffs[k].coeffs)
            n += 1
            p.coeffs[k].coeffs[h] = vals[n]
        end
    end
    return p
end

# Enclosure [lo, hi] (BigFloat) of {p(ξ) + Δ} over all coefficient values:
#   p(ξ) ∈ pmid(ξ) ± prad(|ξ|)   (assumes [mid ± radius] ⊇ coefficient)
function value_box(tm, ξb, tmpl)
    cl = coeflist(polynomial(tm))
    m = evaluate(float_poly(mid.(cl), tmpl), ξb)
    r = evaluate(float_poly(radius.(cl), tmpl), abs.(ξb))
    Δ = remainder(tm)
    return (m - r + big(inf(Δ)), m + r + big(sup(Δ)))
end

corners(boxes) = [collect(w) for w in Iterators.product(boxes...)][:]
box_samples(rng, N, n) = vcat(corners(fill((-1.0, 1.0), N)), [2 .* rand(rng, N) .- 1 for _ in 1:n])

# --- aliasing ----------------------------------------------------------------
function no_alias(v)
    for i in eachindex(v)
        for j in eachindex(v)
            i == j && continue
            v[i] === v[j] && return false
            v[i].pol.coeffs === v[j].pol.coeffs && return false
            for k in eachindex(v[i].pol.coeffs)
                v[i].pol.coeffs[k] === v[j].pol.coeffs[k] && return false
                v[i].pol.coeffs[k].pol.coeffs === v[j].pol.coeffs[k].pol.coeffs && return false
            end
        end
    end
    return true
end
disjoint(a, b) = all(a[i] !== b[j] && a[i].pol.coeffs !== b[j].pol.coeffs
                     for i in eachindex(a), j in eachindex(b))

# --- integration checks ------------------------------------------------------
# last covered time equals tend (exact, given representable steps)
function reached_end(sol, tend)
    n = lastindex(sol)
    return expansion_point(sol, n) + stepsize(domain(sol, n)) == tend
end

# exact start times t_k = t_1 + Σ_{j<k} h_j (BigFloat; columns 1 and 2 start at t_1)
function exact_starts(sol)
    tTM = expansion_point(sol)
    K = lastindex(sol)
    ts = Vector{BigFloat}(undef, K)
    ts[1] = ts[2] = big(tTM[1])
    for k in 3:K
        ts[k] = ts[k-1] + big(stepsize(domain(sol, k-1)))
    end
    return ts
end

# ξ ∈ B (vertices favoured), the exact initial point x0(ξ), τ samples in a step
sample_ξ(rng, N; pcorner = 0.125) =
    rand(rng) < pcorner ? rand(rng, [-1.0, 1.0], N) : 2 .* rand(rng, N) .- 1
x0_of(X0, ξ) = [interval(big(mid(X0[i])) + big(radius(X0[i])) * big(ξ[i])) for i in eachindex(X0)]
tau_samples(rng, d, ntau) =                         # both ends + ntau random points
    vcat(zero(sup(d)), stepsize(d), [inf(d) + rand(rng) * (sup(d) - inf(d)) for _ in 1:ntau])

# xe ⊆ xk -> :pass;  xe ∩ xk = ∅ -> :fail (with the gap);  else :inconclusive
function classify(xe, xk)
    xkB = [interval(BigFloat, inf(x), sup(x)) for x in xk]
    all(issubset_interval.(xe, xkB)) && return (:pass, 0.0)
    any(isdisjoint_interval.(xe, xkB)) || return (:inconclusive, 0.0)
    gap = maximum(Float64(max(inf(xe[i]) - sup(xkB[i]), inf(xkB[i]) - sup(xe[i])))
                  for i in eachindex(xe))
    return (:fail, gap)
end

# Pointwise containment along the whole trajectory, at every step and at
# τ ∈ {0, h_k, ntau random points}:  x(t_k + τ; x0(ξ)) ∈ Φ_k(τ, σ_k(ξ)),
# with σ_2 = ξ, σ_{k+1} = ρ_k(σ_k) (via reconstruct_rho), exact times, BigFloat xe.
# fexact(t::Interval, x0::Vector{<:Interval}) -> exact solution (enclosure).
function pointwise_check(sol, fexact, X0; rng, nsamples = 100, ntau = 10)
    N = length(X0)
    Bx = symmetric_box(N, Float64)
    ρ = reconstruct_rho(sol)
    ts = exact_starts(sol)
    tdrift = Float64(maximum(abs.(ts .- big.(expansion_point(sol)))))
    maxdiam, nfail, ninc = 0.0, 0, 0
    for _ in 1:nsamples
        ξ = sample_ξ(rng, N)
        x0 = x0_of(X0, ξ)
        σ = interval.(ξ)
        for k in 2:lastindex(sol)
            for τ in tau_samples(rng, domain(sol, k), ntau)
                xk = [evaluate(TM.__evaluate_rig!(deepcopy(sol[k][i][0]), sol[k][i], τ), σ)
                      for i in 1:N]
                st, gap = classify(fexact(interval(ts[k] + big(τ)), x0), xk)
                st === :fail && (nfail += 1; @warn "pointwise FAILURE" k τ ξ gap)
                st === :inconclusive && (ninc += 1)
                maxdiam = max(maxdiam, maximum(diam.(xk)))
            end
            # same bounder as qrprecondition_rig!; ∩ B is sound since ρ_k(B) ⊆ B is proved
            σ = [intersect_interval(monotonicity_bounder(ρ[i,k], σ), Bx[i]) for i in 1:N]
            @assert !any(isempty_interval, σ)
        end
    end
    return (; nfail, ninc, maxdiam, tdrift)
end

# Containment in the stored flowpipe boxes (what users get from flowpipe(sol)):
# x(t; x0) ∈ flowpipe(sol, k) for t sampled in [t_k, t_k + h_k].
function flowpipe_check(sol, fexact, X0; rng, nsamples = 100, ntau = 10)
    ts = exact_starts(sol)
    nfail, ninc = 0, 0
    for _ in 1:nsamples
        x0 = x0_of(X0, sample_ξ(rng, length(X0)))
        for k in 2:lastindex(sol)
            box = flowpipe(sol, k)
            for τ in tau_samples(rng, domain(sol, k), ntau)
                st, gap = classify(fexact(interval(ts[k] + big(τ)), x0), box)
                st === :fail && (nfail += 1; @warn "flowpipe FAILURE" k τ gap)
                st === :inconclusive && (ninc += 1)
            end
        end
    end
    return (; nfail, ninc)
end

# Rigorous range of Φ_k(τ, ·) over B, per component
function step_range(sol, k, τ)
    B = symmetric_box(length(sol[k]), Float64)
    return [monotonicity_bounder(TM.__evaluate_rig!(deepcopy(sol[k][i][0]), sol[k][i], τ), B)
            for i in eachindex(sol[k])]
end

# flowpipe width at the final time
function final_width(sol)
    n = lastindex(sol)
    return diam.(step_range(sol, n, stepsize(domain(sol, n))))
end

# 1D: solutions cannot cross, exact reachable set = [φ(t, inf X0), φ(t, sup X0)]
function exact_width_1d(fex, X0, t)
    xe(x) = fex(interval(big(t)), [interval(big(x))])[1]
    return Float64(mid(xe(sup(X0[1])) - xe(inf(X0[1]))))
end

# 1D tightness along the trajectory: enclosure width / exact width at
# τ ∈ {0, h_k/2, h_k} of every step; returns the maximum and the final ratio.
function tightness_profile_1d(sol, fex, X0)
    ts = exact_starts(sol)
    ratios = Float64[]
    for k in 2:lastindex(sol)
        h = stepsize(domain(sol, k))
        for τ in (zero(h), h/2, h)
            push!(ratios, diam(step_range(sol, k, τ)[1]) / exact_width_1d(fex, X0, ts[k] + big(τ)))
        end
    end
    return (; maxratio = maximum(ratios), finalratio = ratios[end])
end

# Completion, pointwise and flowpipe containment, tdrift == 0, and regression
# bounds on the number of steps and on the pointwise enclosure width.
# Returns sol (or nothing). Sampling is seeded (reproducible bounds).
function check_problem(name, f!, fex, X0, tini, tend, orderQ, orderT, abstol, names;
        maxlen, maxdiam, seed = 42, nsamples = 100, printinfo = false, kwargs...)
    sol = run_integ3(f!, X0, tini, tend, orderQ, orderT, abstol,
                     JetSpace(2*orderQ, names); kwargs...)
    sol === nothing && return nothing
    rng = MersenneTwister(seed)
    @test reached_end(sol, tend)
    pw = pointwise_check(sol, fex, X0; rng, nsamples)
    fp = flowpipe_check(sol, fex, X0; rng, nsamples)
    @test pw.nfail == 0
    @test fp.nfail == 0
    @test pw.tdrift == 0
    @test length(sol) <= maxlen
    @test pw.maxdiam <= maxdiam
    printinfo && @info name length(sol) pw.nfail pw.ninc pw.maxdiam pw.tdrift fp.nfail fp.ninc
    return sol
end

# 1D tightness test (prints the profile, to set the bound)
function check_tightness_1d(name, sol, fex, X0; maxratio, printinfo=false)
    sol === nothing && return
    tp = tightness_profile_1d(sol, fex, X0)
    printinfo && @info "$name tightness" tp.maxratio tp.finalratio
    @test tp.maxratio <= maxratio
end

# --- exact step ---------------------------------------------------------------
function check_step(t0, δt, tmax, sgn)
    δt1 = VI._exact_step(t0, δt, tmax, sgn)
    t1 = t0 + δt1
    s, e = VI._two_sum(t0, δt1)
    @test s == t1 && iszero(e)                       # exact
    @test t1 == tmax || abs(δt1) <= abs(δt) ||       # conservative, except landing on
          t1 == (sgn > 0 ? nextfloat(t0) : prevfloat(t0))  # tmax or a sub-ulp request
    @test sign(δt1) == sgn                           # direction
    @test sgn * t1 <= sgn * tmax                     # no overshoot
    return δt1
end

# Repeated steps with the same requested δt until tmax is reached exactly;
# returns the number of steps (every step checked with check_step).
function reach_tmax(t0, δt, tmax, sgn; maxsteps = 64)
    t, n = t0, 0
    while t != tmax
        t += check_step(t, δt, tmax, sgn)
        n += 1
        n < maxsteps || error("reach_tmax: tmax not reached")
    end
    return n
end

# --- QR preconditioning ------------------------------------------------------
# Random v = c + Aξ + h(ξ) + Δ on B (interval coefficients widened by `wid`)
function random_vTMN(rng, X, A, ::Type{T}; nl = 0.05, remr = 1e-6, wid = 1e-9) where {T}
    N = length(X); ordQ = TS.order(X[1])
    B = symmetric_box(N, Float64); zB = zero(B)
    v = Vector{TaylorModelN{N,T,Float64}}(undef, N)
    for i in 1:N
        p = randn(rng) + sum(A[i,j] * X[j] for j in 1:N)
        if nl > 0
            for d in 2:ordQ, _ in 1:3
                p += nl * randn(rng) * X[rand(rng, 1:N)]^(d-1) * X[rand(rng, 1:N)]
            end
        end
        Δ = interval(-remr * rand(rng), remr * rand(rng))
        if T <: Interval
            pI = interval(1.0) * p
            for k in eachindex(pI.coeffs), h in eachindex(pI.coeffs[k].coeffs)
                cf = pI.coeffs[k].coeffs[h]
                w = iszero(mid(cf)) ? 0.0 : wid * rand(rng)
                pI.coeffs[k].coeffs[h] = cf + interval(-w, w)
            end
            v[i] = TM.unsafe_TaylorModelN(pI, Δ, zB, B)
        else
            v[i] = TM.unsafe_TaylorModelN(p, Δ, zB, B)
        end
    end
    return v
end

# :ok iff finite, ρ(B) ⊆ B, (⋆) M_L^{-1}(w - c) ∈ B and (c) M_L^{-1}(w - c) ∈ ρ(ξ)
# for all corners w of V(ξ) at all samples (exact in BigFloat)
function check_precond(v, rng, tmpl; nsamples = 50)
    N = length(v)
    left, right = zero.(v), zero.(v)
    CT = typeof(v[1].pol.coeffs[1].coeffs[1])
    VI.qrprecondition_rig!(left, right, zeros(CT, N, N), zeros(Interval{Float64}, N), zeros(N), v) ||
        return :declined
    all(isfin, vcat((coeflist(polynomial(tm)) for tm in vcat(left, right))...)) || return :nonfinite
    VI._right_in_box(right) || return :rho_not_in_B
    ML = big.([mid(left[i].pol.coeffs[2].coeffs[j]) for i in 1:N, j in 1:N])
    c  = big.([mid(left[i].pol.coeffs[1].coeffs[1]) for i in 1:N])
    for ξ in box_samples(rng, N, nsamples)
        ξb = big.(ξ)
        Vb = [value_box(v[i], ξb, tmpl) for i in 1:N]
        Rb = [value_box(right[i], ξb, tmpl) for i in 1:N]
        for w in corners(Vb)
            s = ML \ (w .- c)
            maximum(abs.(s)) <= 1 || return :star_violated
            all(Rb[i][1] <= s[i] <= Rb[i][2] for i in 1:N) || return :compose_violated
        end
    end
    return :ok
end

rotation(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


# ============================================================================
# Tests
# ============================================================================

@testset "validated_integ3" begin
    setprecision(BigFloat, 256)
    local printinfo = false

    @testset "_exact_step" begin
        s, e = VI._two_sum(1.0, 1e-17)
        @test big(s) + big(e) == big(1.0) + big(1e-17)
        @test check_step(0.0, 0.1, 10.0, 1) == 0.1                   # t0 = 0
        check_step(4.0, 0.8802595366450658, 10.0, 1)                 # |t0| ≥ |δt|
        check_step(10.0, -0.3, 0.0, -1)
        check_step(1e-3, 0.7, 10.0, 1)                              # 0 < |t0| < |δt|
        check_step(-1e-3, -0.7, -10.0, -1)
        @test check_step(7565.0, 8e-14, 8000.0, 1) == eps(7565.0)   # request < ½ ulp(t0)
        # requested step larger than what is left (TI's estimate is often larger):
        # the step is clipped and tmax is reached exactly, possibly in a few steps
        # when tmax - t0 is not representable (tmax not within a factor 2 of t0)
        @test reach_tmax(0.0, 0.5, 0.7, 1) == 2                     # 0.5, then 0.2 → 0.7
        @test reach_tmax(4.0, 9.0, 7.3, 1) == 1                     # Sterbenz: one step
        @test reach_tmax(5.0, -9.0, 0.1, -1) <= 8                   # 0.1 - 5 not representable
        @test reach_tmax(5.0, -20.0, -3.0, -1) <= 8                 # across 0
        @test reach_tmax(1e-3, 10.0, 7.3, 1) <= 16                  # 0 < |t0| ≪ |δt|
        rng = MersenneTwister(2026)
        for _ in 1:2_000
            sgn = rand(rng, (-1, 1))
            t0 = sgn * 10 * rand(rng) * rand(rng, (0.0, 1e-6, 1.0, 1e3))
            δt = sgn * rand(rng) * rand(rng, (1e-8, 1e-2, 1.0))
            check_step(t0, δt, t0 + sgn * 100, sgn)
        end
    end

    @testset "picard_lindelof" begin
        function fb!(dx, x, p, t)                  # autonomous
            dx[1] = x[2]
            dx[2] = -one(x[1])
            nothing
        end
        function nonaut!(dx, x, p, t)              # non-autonomous (catches arg swaps)
            dx[1] = t + zero(x[1])
            dx[2] = x[1]
            nothing
        end
        orderQ, orderT = 2, 4
        X0 = SVector(interval(9.75, 10.25), interval(-0.25, 0.25))
        with_default_space(JetSpace(2*orderQ, ["ξ₁", "ξ₂"])) do sp
            c, x1N, dx1N, t1N = setup(fb!, X0, 0.0, 0.5, orderQ, orderT, sp)
            @test no_alias(c.x1N) && no_alias(c.dx1N) && no_alias(c.x2N)
            @test disjoint(c.x1N, c.dx1N) && disjoint(c.x1N, c.x2N) && disjoint(c.dx1N, c.x2N)

            # no aliasing, independence; autonomous iterate P = (x0₁ + x0₂τ, x0₂ - τ)
            x10, x20 = cpoly(x1N[1], 0), cpoly(x1N[2], 0)
            x2N = VI.picard_lindelof(fb!, dx1N, x1N, t1N, nothing)
            @test no_alias(x2N) && disjoint(x2N, x1N) && disjoint(x2N, dx1N)
            @test cpoly(x2N[1], 0) == x10 && cpoly(x2N[1], 1) == x20
            @test cpoly(x2N[2], 0) == x20 && cpoly(x2N[2], 1) == constpoly(-1.0, x10)
            @test all(iszero(cpoly(x2N[i], k)) for i in 1:2 for k in 2:orderT)
            before = deepcopy(cpoly(x2N[2], 0))
            x2N[1].pol.coeffs[1].pol.coeffs[1].coeffs[1] += 1.0
            @test cpoly(x2N[2], 0) == before

            # non-autonomous iterate, t0 = 2: P = (x0₁ + t0τ + τ²/2, x0₂ + x0₁τ)
            t0 = 2.0
            _, x1N, dx1N, t1N = setup(nonaut!, X0, t0, 0.5, orderQ, orderT, sp)
            x10, x20 = cpoly(x1N[1], 0), cpoly(x1N[2], 0)
            x1N_copy = deepcopy(x1N)
            x2N = VI.picard_lindelof(nonaut!, dx1N, x1N, t1N, nothing)
            @test cpoly(x2N[1], 0) == x10
            @test cpoly(x2N[1], 1) == constpoly(t0, x10)
            @test cpoly(x2N[1], 2) == constpoly(0.5, x10)
            @test cpoly(x2N[2], 0) == x20 && cpoly(x2N[2], 1) == x10
            @test all(iszero(cpoly(x2N[2], k)) for k in 2:orderT)

            # allocating == in-place; input not modified
            y2N = zero.(x1N)
            VI.picard_lindelof!(nonaut!, y2N, dx1N, x1N, t1N, nothing)
            @test all(cpoly(x2N[i], k) == cpoly(y2N[i], k) for i in 1:2, k in 0:orderT)
            @test all(isequal_interval(remainder(x2N[i]), remainder(y2N[i])) for i in 1:2)
            @test all(cpoly(x1N[i], k) == cpoly(x1N_copy[i], k) for i in 1:2, k in 0:orderT)
        end
    end

    @testset "monotonicity_bounder" begin
        with_default_space(JetSpace(6, ["σ"])) do sp
            σ = variables(sp; order = 3)[1]
            B = symmetric_box(1); zB = zero(B)
            # linear-dominated 1D: exact range [p(-1), p(1)] + Δ
            tm = TM.unsafe_TaylorModelN(0.5 + 0.36σ - 0.03σ^2 - 0.008σ^3,
                                        interval(-1e-10, 2e-10), zB, B)
            r = monotonicity_bounder(tm, B)
            @test issubset_interval(r, evaluate(tm, B))
            @test isapprox(inf(r), 0.118 - 1e-10; atol = 1e-14)
            @test isapprox(sup(r), 0.822 + 2e-10; atol = 1e-14)
            # non-monotone: falls back to (and contains the range of) plain evaluation
            tm2 = TM.unsafe_TaylorModelN(0.0 + σ^2, interval(0.0), zB, B)
            r2 = monotonicity_bounder(tm2, B)
            @test issubset_interval(interval(0.0, 1.0), r2)
            @test issubset_interval(r2, evaluate(tm2, B))
        end
        # 2D soundness on random polynomials (BigFloat samples)
        with_default_space(JetSpace(8, ["σ₁", "σ₂"])) do sp
            X = variables(sp; order = 4)
            rng = MersenneTwister(7)
            ok = true
            for _ in 1:50
                v = random_vTMN(rng, X, randn(rng, 2, 2), Float64; nl = 0.3)
                for tm in v
                    r = monotonicity_bounder(tm, symmetric_box(2))
                    ok &= issubset_interval(r, evaluate(tm, symmetric_box(2)))
                    for ξ in box_samples(rng, 2, 20)
                        val = evaluate(polynomial(tm), big.(ξ))
                        ok &= big(inf(r)) <= val + big(inf(remainder(tm))) &&
                              val + big(sup(remainder(tm))) <= big(sup(r))
                    end
                end
            end
            @test ok
        end
    end

    @testset "qrprecondition_rig!" begin
        with_default_space(JetSpace(8, ["σ₁", "σ₂"])) do sp
            X = variables(sp; order = 4)
            tmpl = 0.0 * X[1]
            cases = [
                ("well-conditioned", rng -> randn(rng, 2, 2),                         0.05, 1e-6),
                ("ill-conditioned",  rng -> rotation(rand(rng)) * Diagonal([1, 1e-8]), 0.05, 1e-6),
                ("strong nonlinear", rng -> randn(rng, 2, 2),                         0.5,  1e-6),
                ("point component",  rng -> [1.0 0.0; 0.3 0.0],                       0.0,  0.0),
            ]
            for T in (Float64, Interval{Float64}), (name, Agen, nl, remr) in cases
                @testset "$name, T = $T" begin
                    rng = MersenneTwister(1234)
                    counts = Dict{Symbol,Int}()
                    for _ in 1:100
                        st = check_precond(random_vTMN(rng, X, Agen(rng), T; nl, remr), rng, tmpl)
                        counts[st] = get(counts, st, 0) + 1
                    end
                    @test get(counts, :ok, 0) == 100
                    get(counts, :ok, 0) == 100 || @info "qrprecondition_rig!" name T counts
                end
            end
        end
    end

    @testset "integration: problems with known solution" begin
        # Regression bounds (maxlen, maxdiam, maxratio): ~15-50% above observed values.

        @testset "falling ball (2D, linear shear)" begin
            TI.@taylorize2 function falling_ball!(dx, x, p, t)
                dx[1] = x[2]
                dx[2] = -one(x[1])
                nothing
            end
            fex(t, x0) = [x0[1] + x0[2]*t - t^2/2, x0[2] - t]
            X0 = [10.0, 0.0] .+ 0.25 .* symmetric_box(2)
            check_problem("falling ball", falling_ball!, fex, X0, 0.0, 10.0, 2, 4, 1e-20,
                          ["ξₓ", "ξᵥ"]; maxlen = 8, maxdiam = 1e-13, printinfo)
        end

        @testset "x_square (1D)" begin
            TI.@taylorize2 x_square!(dx, x, p, t) = (dx[1] = x[1]^2;)
            fex(t, x0) = [1 / (1/x0[1] - t)]               # blow-up at t = 1/x0 ≥ 0.4848
            X0 = [2.0] .+ 0.0625 .* symmetric_box(1)
            sol = check_problem("x_square", x_square!, fex, X0, 0.0, 0.45, 2, 10, 1e-15, ["ξ"];
                                maxlen = 160, maxdiam = 0.02, printinfo)
            check_tightness_1d("x_square", sol, fex, X0; maxratio = 1 + 1e-4, printinfo)
        end

        @testset "x_cube (1D)" begin
            TI.@taylorize2 x_cube!(dx, x, p, t) = (dx[1] = - x[1]^3;)
            fex(t, x0) = [x0[1] / sqrt(1 + 2*x0[1]^2*t)]
            X0 = [0.5] .+ 0.4 .* symmetric_box(1)
            sol = check_problem("x_cube", x_cube!, fex, X0, 0.0, 3.0, 3, 20, 1e-20, ["ξ"];
                                maxlen = 13, maxdiam = 0.04,
                                adaptive = true, minabstol = 1e-50, printinfo)
            check_tightness_1d("x_cube", sol, fex, X0; maxratio = 1.1, printinfo)
        end

        @testset "Christian's test (1D)" begin
            TI.@taylorize2 ff!(dx, x, p, t) = (dx[1] = x[1]*(one(x[1])-x[1]^2);)
            fex(t, x0) = [x0[1]*exp(t) / sqrt(1 + x0[1]^2*(exp(2t) - 1))]
            X0 = [0.5] .+ 0.3 .* symmetric_box(1)
            sol = check_problem("Christian", ff!, fex, X0, 0.0, 3.0, 3, 18, 1e-20, ["ξ"];
                                maxlen = 35, maxdiam = 0.014,
                                adaptive = true, minabstol = 1e-50, printinfo)
            check_tightness_1d("Christian", sol, fex, X0; maxratio = 1.05, printinfo)
        end

        @testset "cos(t) (1D, non-autonomous)" begin
            TI.@taylorize2 function cost!(du, u, params, t)
                tt = t + zero(u[1])
                du[1] = cos(tt)
                return nothing
            end
            fex(t, x0) = [sin(t) + x0[1]]
            X0 = [0.0] .+ 0.1 .* symmetric_box(1)
            sol = check_problem("cos(t)", cost!, fex, X0, 0.0, 5.0, 1, 15, 1e-15, ["ξ"];
                                maxlen = 13, maxdiam = 1e-14,
                                adaptive = true, minabstol = 1e-50, printinfo)
            check_tightness_1d("cos(t)", sol, fex, X0; maxratio = 1 + 1.e-12, printinfo)
        end

        @testset "limit cycle (2D, rotation + nonlinear)" begin
            # x' = x - y - x r², y' = x + y - y r²; in polar r' = r(1-r²), θ' = 1
            function limcyc!(dx, x, p, t)
                r2 = x[1]^2 + x[2]^2
                dx[1] = x[1] - x[2] - x[1]*r2
                dx[2] = x[1] + x[2] - x[2]*r2
                nothing
            end
            function fex(t, x0)
                r02 = x0[1]^2 + x0[2]^2
                s = 1 / sqrt(r02 + (1 - r02) * exp(-2t))       # r(t)/r0
                return [s * (x0[1]*cos(t) - x0[2]*sin(t)), s * (x0[1]*sin(t) + x0[2]*cos(t))]
            end
            X0 = [1.2, 0.0] .+ 0.05 .* symmetric_box(2)
            sol = check_problem("limit cycle", limcyc!, fex, X0, 0.0, 4π, 3, 15, 1e-20,
                                ["ξ₁", "ξ₂"]; maxlen = 140, maxdiam = 1e-4, nsamples = 30,
                                printinfo)
            if sol !== nothing
                w = final_width(sol)
                printinfo && @info "limit cycle final width" w
                @test w[1] < 2.2e-3 && w[2] < 0.12
            end
        end
    end
end