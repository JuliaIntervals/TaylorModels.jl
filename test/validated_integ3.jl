# Tests for validated_integ3 and its building blocks.
#
#   1. picard_lindelof        no aliasing, exact Picard iterates, argument order
#   2. cache                  the Picard time of the first step is t0
#   3. monotonicity_bounder   range enclosure (soundness, tightness)
#   4. qrprecondition!    (⋆) and v ∈ L∘ρ, exact in BigFloat, T ∈ {Float64, Interval}
#   5. integration            completion, pointwise containment x(t;x0) ∈ Φ_k(τ, σ_k(ξ))
#                             with exact times, tdrift == 0, tightness
#
# The arithmetic of TaylorModel1{TaylorModelN} (products, squares, numbers, integration) is
# tested in `TM1_TMN.jl`. Every problem uses its own `JetSpace`; the default space is not set.
# `TM_TEST_MIDRAD=true` runs the integrations with the midpoint-radius kernels.
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

# Settings of the integrator in the tests
const _TEST_MIDRAD = get(ENV, "TM_TEST_MIDRAD", "false") == "true"

# validated_integ3 with no warnings allowed
run_integ3(f!, X0, tini, tend, orderQ, orderT, abstol, sp; kwargs...) =
    @test_logs min_level=Logging.Warn validated_integ3(
        f!, X0, tini, tend, orderQ, orderT, abstol, sp;
        parse_eqs = true, maxsteps = 2000, midrad = _TEST_MIDRAD, kwargs...)

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

# xe ⊆ xk -> :pass;  xe ∩ xk = ∅ -> :fail (with the gap and the gap in ulps of
# the missed endpoint of xk, Float64);  else :inconclusive
function classify(xe, xk)
    xkB = [interval(BigFloat, inf(x), sup(x)) for x in xk]
    all(issubset_interval.(xe, xkB)) && return (:pass, 0.0, 0.0)
    any(isdisjoint_interval.(xe, xkB)) || return (:inconclusive, 0.0, 0.0)
    gaps = [max(inf(xe[i]) - sup(xkB[i]), inf(xkB[i]) - sup(xe[i])) for i in eachindex(xe)]
    i = argmax(gaps)
    missed = inf(xe[i]) > sup(xkB[i]) ? sup(xk[i]) : inf(xk[i])
    return (:fail, Float64(gaps[i]), Float64(gaps[i]) / eps(Float64(missed)))
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
                xk = [evaluate(TM.__evaluate!(deepcopy(sol[k][i][0]), sol[k][i], τ), σ)
                      for i in 1:N]
                st, gap, ulps = classify(fexact(interval(ts[k] + big(τ)), x0), xk)
                st === :fail && (nfail += 1; @warn "pointwise FAILURE" k τ ξ gap ulps)
                st === :inconclusive && (ninc += 1)
                maxdiam = max(maxdiam, maximum(diam.(xk)))
            end
            # same bounder as qrprecondition!; ∩ B is sound since ρ_k(B) ⊆ B is proved
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
                st, gap, ulps = classify(fexact(interval(ts[k] + big(τ)), x0), box)
                st === :fail && (nfail += 1; @warn "flowpipe FAILURE" k τ gap ulps)
                st === :inconclusive && (ninc += 1)
            end
        end
    end
    return (; nfail, ninc)
end

# Rigorous range of Φ_k(τ, ·) over B, per component
function step_range(sol, k, τ)
    B = symmetric_box(length(sol[k]), Float64)
    return [monotonicity_bounder(TM.__evaluate!(deepcopy(sol[k][i][0]), sol[k][i], τ), B)
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
        maxlen, maxdiam, seed = 42, nsamples = 100, printinfo = false,
        testbroken = false, kwargs...)
    sol = run_integ3(f!, X0, tini, tend, orderQ, orderT, abstol,
                     JetSpace(2*orderQ, names); kwargs...)
    sol === nothing && return nothing
    rng = MersenneTwister(seed)
    @test reached_end(sol, tend)
    pw = pointwise_check(sol, fex, X0; rng, nsamples)
    fp = flowpipe_check(sol, fex, X0; rng, nsamples)
    if testbroken
        @test_broken pw.nfail == 0
    else
        @test pw.nfail == 0
    end
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
    VI.qrprecondition!(left, right, zeros(CT, N, N), zeros(Interval{Float64}, N), zeros(N), v) ||
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
        let sp = JetSpace(2*orderQ, ["ξ₁", "ξ₂"])
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

    @testset "cache: the Picard time of the first step is t0" begin
        function nonaut!(dx, x, p, t)
            dx[1] = t + zero(x[1])
            dx[2] = x[1]
            nothing
        end
        orderQ, orderT = 1, 4
        X0 = SVector(interval(1.0, 1.1), interval(-0.25, 0.25))
        let sp = JetSpace(2*orderQ, ["ξ₁", "ξ₂"])
            c = VI.init_cache_VI3(nonaut!, 2.0, X0, 2, orderT, orderQ, sp, nothing;
                                  parse_eqs = true)
            @test isequal_interval(constant_term(polynomial(c.t1NI[0])), interval(2.0))   # it was 0
            @test constant_term(polynomial(c.t1N[0])) == 2.0
        end
    end

    @testset "monotonicity_bounder" begin
        let sp = JetSpace(6, ["σ"])
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
        let sp = JetSpace(8, ["σ₁", "σ₂"])
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

    @testset "qrprecondition!" begin
        let sp = JetSpace(8, ["σ₁", "σ₂"])
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
                    get(counts, :ok, 0) == 100 || @info "qrprecondition!" name T counts
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
            # elapsed time since tini
            τ(t) = t - interval(big(tini))
            fex(t, x0) = [x0[1] + x0[2]*τ(t) - τ(t)^2/2, x0[2] - τ(t)]

            # Forward integration
            tini, tend = 0.0, 10.0
            X0 = [10.0, 0.0] .+ 0.25 .* symmetric_box(2)
            check_problem("falling ball", falling_ball!, fex, X0, tini, tend, 2, 4, 1e-20,
                          ["ξₓ", "ξᵥ"]; maxlen = 8, maxdiam = 1e-13, printinfo)

            # Backward integration
            tini, tend = 10.0, 0.0
            X0 = [-40.0, -10.0] .+ 0.25 .* symmetric_box(2)                # state at t = 10
            check_problem("falling ball backward", falling_ball!, fex, X0, tini, tend,
                          2, 4, 1e-20, ["ξₓ", "ξᵥ"]; maxlen = 9, maxdiam = 1e-13, printinfo)
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
            # Exact sol writen aiming numerical stability
            fex(t, x0) = [x0[1] +
                2 * cos((t + interval(big(tini)))/2) * sin((t - interval(big(tini)))/2)]

            # Forward integration
            tini, tend = 0.0, 5.0
            X0 = [0.0] .+ 0.1 .* symmetric_box(1)
            sol = check_problem("cos(t)", cost!, fex, X0, tini, tend, 1, 15, 1e-15, ["ξ"];
                                maxlen = 13, maxdiam = 1e-14,
                                adaptive = true, minabstol = 1e-50, printinfo)
            check_tightness_1d("cos(t)", sol, fex, X0; maxratio = 1 + 1.e-12, printinfo)

            # Backward integration
            tini, tend = 5.0, 0.0
            X0 = [0.0] .+ 0.1 .* symmetric_box(1)
            # Formerly @test_broken: the Picard step ignored rounding of the retained
            # coefficients (and used t = τ instead of t₀ + τ in the first step)
            sol = check_problem("cos(t) backward", cost!, fex, X0, tini, tend, 1, 15, 1e-15,
                                ["ξ"]; maxlen = 13, maxdiam = 1e-13,
                                adaptive = true, minabstol = 1e-50, printinfo)
            check_tightness_1d("cos(t) backward", sol, fex, X0; maxratio = 1 + 1e-12, printinfo)
        end

        @testset "polynomial chain (3D): rounding coefficients (τ³/6)" begin
            # x₁' = x₂, x₂' = x₃, x₃' = -1: the flow is polynomial, but the τ³
            # coefficient -1/6 is not representable in Float64. orderT = 5: the step
            # size estimate uses the coefficients of orders orderT-1 and orderT, which
            # vanish here (with orderT = 4 the τ³ coefficient forces h ≈ 2e-5).
            TI.@taylorize2 function chain3!(dx, x, p, t)
                dx[1] = x[2]
                dx[2] = x[3]
                dx[3] = -one(x[1])
                nothing
            end
            τ(t) = t - interval(big(tini))
            fex(t, x0) = [x0[1] + x0[2]*τ(t) + x0[3]*τ(t)^2/2 - τ(t)^3/6,
                          x0[2] + x0[3]*τ(t) - τ(t)^2/2,
                          x0[3] - τ(t)]
            names = ["ξ₁", "ξ₂", "ξ₃"]
            # NOTE: maxlen/maxdiam are placeholders; set them from `printinfo = true`
            # Forward
            tini, tend = 0.0, 3.0
            X0 = [1.0, -0.5, 0.3] .+ 0.1 .* symmetric_box(3)
            check_problem("chain3", chain3!, fex, X0, tini, tend, 1, 5, 1e-15, names;
                          maxlen = 20, maxdiam = 1e-8, printinfo)
            # Backward (X0 is the state at t = 3)
            tini, tend = 3.0, 0.0
            check_problem("chain3 backward", chain3!, fex, X0, tini, tend, 1, 5, 1e-15,
                          names; maxlen = 20, maxdiam = 1e-8, printinfo)
        end

        # Issue #208
        @testset "Geodesic flow on the sphere (6D, point components)" begin
            # x' = v, v' = -(|v|²/|x|²) x. From x0 = e_z with v0 = (v1, v2, 0), c = |v0|, the solution is
            # x(t) = cos(ct) e_z + sin(ct) v0/c and v(t) = x'(t), for every point of the box. Four components
            # of X0 are points (rank-deficient directions of the preconditioning). The stored models are
            # functions of σ_k, not of ξ (σ_2 = ξ, σ_{k+1} = ρ_k(σ_k)): the pointwise check composes with ρ,
            # while evaluating the models of the later steps at ξ would (wrongly) fail.
            function spray!(dx, x, p, t)
                q = (x[4]^2 + x[5]^2 + x[6]^2) / (x[1]^2 + x[2]^2 + x[3]^2)
                dx[1] = x[4]; dx[2] = x[5]; dx[3] = x[6]
                dx[4] = -q * x[1]; dx[5] = -q * x[2]; dx[6] = -q * x[3]
                nothing
            end
            function fex(t, x0)
                v1, v2 = x0[4], x0[5]
                c = sqrt(v1^2 + v2^2)
                s, co = sin(c * t), cos(c * t)
                return [s * v1 / c, s * v2 / c, co, co * v1, co * v2, -c * s]
            end
            A, B = 0.7, 0.7 + 2π / 256
            X0 = [interval(0.0), interval(0.0), interval(1.0),
                  interval(cos(B), cos(A)), interval(sin(A), sin(B)), interval(0.0)]
            # NOTE: maxlen/maxdiam are placeholders (3 steps and remainders of 1e-8 to 1e-6 were observed);
            # set them from `printinfo = true`
            check_problem("geodesic flow", spray!, fex, X0, 0.0, 0.5, 2, 12, 1e-15,
                          ["ξ$i" for i in 1:6]; maxlen = 6, maxdiam = 1e-4, nsamples = 30,
                          printinfo, maxsteps=100)
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