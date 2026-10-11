# Tests for the arithmetic of TaylorModel1{TaylorModelN{N,Interval{S},S},S}
# (the type used by the Picard operator of `validated_integ3`).
#
# Pointwise soundness: a TM1{TMN} represents every function
#     f(τ, σ) = Σ_i c_i(σ) τ^i + r ,   c_i(σ) ∈ a_i(σ) (TMN remainder included),
#                                       r ∈ remainder(a)  (τ ∈ dom, σ ∈ [-1,1]^N).
# For an operation `op`, we draw such functions for the operands, evaluate `op`
# on the drawn values (thin intervals) and check that the result is contained in the
# evaluation of the model `op(a, b)` at (τ, σ). A wrong sign or scaling of a remainder,
# a missing tail, or a number added to the wrong coefficient show up as failures.
#
# Note: the TMN coefficients are `Interval`s, so rounding is enclosed.

using TaylorModels, TaylorModels.ValidatedInteg
using Test, Random, StaticArrays

const TM = TaylorModels
const VI = TaylorModels.ValidatedInteg
const TAUMAX = 1.0                              # the τ-domain of the models is [0, TAUMAX]
const NSAMPLES = 40                             # samples (σ, τ) per check

# Run `f()` with the given settings of the arithmetic of TaylorModel1{TaylorModelN{Interval}}
# (default: the current ones); they are restored afterwards.
function with_settings(f; midrad = TM._TMN_MIDRAD[], joint = TM._TM1_JOINT[],
        acc = TM._TM1_ACC[], unit = TM._TM1_UNIT[], tail_exact = TM._TMN_TAIL_EXACT[])
    old = (TM._TMN_MIDRAD[], TM._TM1_JOINT[], TM._TM1_ACC[], TM._TM1_UNIT[], TM._TMN_TAIL_EXACT[])
    try
        TM._TMN_MIDRAD[], TM._TM1_JOINT[], TM._TM1_ACC[], TM._TM1_UNIT[] = midrad, joint, acc, unit
        TM._TMN_TAIL_EXACT[] = tail_exact
        return f()
    finally
        TM._TMN_MIDRAD[], TM._TM1_JOINT[], TM._TM1_ACC[], TM._TM1_UNIT[], TM._TMN_TAIL_EXACT[] = old
    end
end

# fill `tm` (in place) with random interval coefficients, TMN remainders and TM1 remainder
function random_tm1!(tm, rng; scale = 1.0, remc = 1e-4, remt = 1e-3)
    tm.dom = interval(0.0, TAUMAX); tm.x0 = 0.0
    tm.rem = interval(-remt*rand(rng), remt*rand(rng))
    for k in eachindex(tm.pol.coeffs)
        coef = tm.pol.coeffs[k]
        for q in eachindex(coef.pol.coeffs), j in eachindex(coef.pol.coeffs[q].coeffs)
            coef.pol.coeffs[q].coeffs[j] = interval(scale * randn(rng))
        end
        coef.rem = interval(-remc*rand(rng), remc*rand(rng))
    end
    return tm
end

pick(rng, I) = clamp(inf(I) + rand(rng) * (sup(I) - inf(I)), inf(I), sup(I))

# draw a represented function of `tm` at σ: values c_i(σ) and a constant r
function draw(tm, σ, rng)
    σI = interval.(σ)
    v = [pick(rng, tm.pol[i](σI)) for i in 0:TS.order(tm.pol)]
    return (; v, r = pick(rng, remainder(tm)))
end
# its value (thin interval enclosure) at τ
fval(s, τI) = sum(interval(s.v[i+1]) * τI^i for i in 0:length(s.v)-1) + interval(s.r)

# evaluation of the model at (τ, σ)
model_val(res, τ, σ) = evaluate(res, interval(τ))(interval.(σ))

# samples (σ, τ): τ at both ends and random
function sample_στ(rng, n)
    out = Tuple{Vector{Float64},Float64}[]
    for k in 1:n
        τ = k == 1 ? 0.0 : (k == 2 ? TAUMAX : TAUMAX * rand(rng))
        push!(out, (2 .* rand(rng, 2) .- 1, τ))
    end
    return out
end

# number of the samples (σ, τ) with  truth ⊄ model,  (truth, model) = case(σ, τ)
function count_failures(case, rng; name = "", n = NSAMPLES)
    nf = 0
    for (σ, τ) in sample_στ(rng, n)
        t, m = case(σ, τ)
        if !issubset_interval(t, m)
            nf += 1
            nf == 1 && @warn "pointwise failure" name σ τ truth = t model = m
        end
    end
    return nf
end

# the same for  truth(f_a, f_b) ⊆ model(a, b)(τ, σ),  f_a, f_b represented functions of a, b
function count_failures_binary(model, truth, a, b, rng; name = "")
    res = model(a, b)
    return count_failures(rng; name) do σ, τ
        τI = interval(τ)
        sa, sb = draw(a, σ, rng), draw(b, σ, rng)
        return truth(fval(sa, τI), fval(sb, τI)), model_val(res, τ, σ)
    end
end

@testset "TaylorModel1{TaylorModelN{Interval}} arithmetic" begin
    orderQ, orderT = 2, 4
    X0 = SVector(interval(1.0, 1.1), interval(-0.25, 0.25))
    function f!(dx, x, p, t)
        dx[1] = x[2]
        dx[2] = -one(x[1])
        nothing
    end
    let sp = JetSpace(2*orderQ, ["ξ₁", "ξ₂"])
        rng = MersenneTwister(2024)
        c = VI.init_cache_VI3(f!, 0.0, X0, 2, orderT, orderQ, sp, nothing; parse_eqs = true)
        a = random_tm1!(c.x1NI[1], rng)
        b = random_tm1!(c.x1NI[2], rng)

        @testset "TM1 ± TM1, TM1 * TM1, square, powers" begin
            @test count_failures_binary(+, +, a, b, rng; name = "a + b") == 0
            @test count_failures_binary(-, -, a, b, rng; name = "a - b") == 0
            @test count_failures_binary(*, *, a, b, rng; name = "a * b") == 0
            @test count_failures_binary((x, y) -> TS.square(x), (x, y) -> x^2, a, b, rng; name = "a²") == 0
            @test count_failures_binary((x, y) -> x^3, (x, y) -> x * x * x, a, b, rng; name = "a³") == 0
            @test count_failures_binary((x, y) -> x^4, (x, y) -> (x * x) * (x * x), a, b, rng; name = "a⁴") == 0
        end

        @testset "TM1 ± number, number ± TM1, TM1 / number" begin
            for num in (0.7, 2, -1.3)
                nI = interval(num)
                @test count_failures_binary((x, y) -> x + num, (x, y) -> x + nI, a, b, rng; name = "a + $num") == 0
                @test count_failures_binary((x, y) -> x - num, (x, y) -> x - nI, a, b, rng; name = "a - $num") == 0
                @test count_failures_binary((x, y) -> num + x, (x, y) -> nI + x, a, b, rng; name = "$num + a") == 0
                @test count_failures_binary((x, y) -> num - x, (x, y) -> nI - x, a, b, rng; name = "$num - a") == 0
                @test count_failures_binary((x, y) -> x / num, (x, y) -> x / nI, a, b, rng; name = "a / $num") == 0
            end
            # polynomial part: the number only goes to the τ⁰ coefficient (Float64 TMN)
            y = c.x1N[1]
            for k in 0:orderT
                y.pol.coeffs[k+1].pol.coeffs[1].coeffs[1] = float(k + 1)
            end
            z, w = y - 1.0, 1.0 - y
            @test constant_term(polynomial(z[0])) == 0.0 && constant_term(polynomial(w[0])) == 0.0
            @test all(constant_term(polynomial(z[k])) == k + 1 for k in 1:orderT)
            @test all(constant_term(polynomial(w[k])) == -(k + 1) for k in 1:orderT)
        end

        @testset "TM1 * TMN" begin
            w = deepcopy(b.pol[0])                       # TMN with remainder
            res = a * w
            nf = count_failures(rng; name = "a * w") do σ, τ
                sa = draw(a, σ, rng)
                wv = pick(rng, w(interval.(σ)))
                return fval(sa, interval(τ)) * interval(wv), model_val(res, τ, σ)
            end
            @test nf == 0
        end

        @testset "integrate! with integration constant" begin
            c0 = deepcopy(b.pol[0])
            c0.rem = interval(-1e-3, 2e-3)               # c0.rem must be kept
            res = deepcopy(a)
            TM.integrate!(res, a, c0, domain(a[0]))
            nf = count_failures(rng; name = "integrate!") do σ, τ
                τI = interval(τ)
                sa = draw(a, σ, rng)
                wv = pick(rng, c0(interval.(σ)))
                t = interval(wv) + interval(sa.r) * τI +
                    sum(interval(sa.v[i+1]) * τI^(i+1) / interval(i+1) for i in 0:orderT)
                return t, model_val(res, τ, σ)
            end
            @test nf == 0
        end

        @testset "integrate! (in place) equals integrate of the polynomial" begin
            for tm in (a, c.x1N[1])                       # Interval and Float64 coefficients
                res = deepcopy(tm)
                TM.integrate!(res, tm, domain(tm[0]))
                ip = integrate(tm.pol)
                same = true
                for k in eachindex(ip.coeffs), q in eachindex(ip.coeffs[k].pol.coeffs),
                        h in eachindex(ip.coeffs[k].pol.coeffs[q].coeffs)
                    same &= isequal_interval(interval(res.pol.coeffs[k].pol.coeffs[q].coeffs[h]),
                                             interval(ip.coeffs[k].pol.coeffs[q].coeffs[h]))
                end
                @test same
                @test all(isequal_interval(res.pol.coeffs[k].rem, ip.coeffs[k].rem) for k in eachindex(ip.coeffs))
            end
        end

        @testset "_copy_tm1: an independent copy" begin
            y = c.x1N[1]                                  # Float64 coefficients
            for k in eachindex(y.pol.coeffs), q in eachindex(y.pol.coeffs[k].pol.coeffs),
                    h in eachindex(y.pol.coeffs[k].pol.coeffs[q].coeffs)
                y.pol.coeffs[k].pol.coeffs[q].coeffs[h] = 10k + q + h / 10
            end
            y.rem = interval(-1e-3, 2e-3); y.pol.coeffs[2].rem = interval(-1e-5, 1e-5)
            y.x0 = 0.25; y.dom = interval(0.0, 0.5)
            cp = VI._copy_tm1(y)
            @test typeof(cp) == typeof(y)
            @test isequal_interval(cp.rem, y.rem) && cp.x0 == y.x0 && isequal_interval(cp.dom, y.dom)
            same = true
            for k in eachindex(y.pol.coeffs)
                same &= isequal_interval(cp.pol.coeffs[k].rem, y.pol.coeffs[k].rem) &&
                        cp.pol.coeffs[k].x0 == y.pol.coeffs[k].x0 &&
                        all(isequal_interval.(cp.pol.coeffs[k].dom, y.pol.coeffs[k].dom))
                for q in eachindex(y.pol.coeffs[k].pol.coeffs)
                    same &= cp.pol.coeffs[k].pol.coeffs[q].coeffs == y.pol.coeffs[k].pol.coeffs[q].coeffs &&
                            cp.pol.coeffs[k].pol.coeffs[q].coeffs !== y.pol.coeffs[k].pol.coeffs[q].coeffs
                end
                same &= cp.pol.coeffs[k] !== y.pol.coeffs[k]
            end
            @test same
            # modifying the original does not change the copy
            before = cp.pol.coeffs[2].pol.coeffs[1].coeffs[1]
            y.pol.coeffs[2].pol.coeffs[1].coeffs[1] += 1.0
            @test cp.pol.coeffs[2].pol.coeffs[1].coeffs[1] == before
        end

        @testset "zero and one of a TaylorModelN" begin
            for w in (a.pol[1], c.x1N[1].pol[2])          # Interval and Float64 coefficients
                z, o = zero(w), one(w)
                @test typeof(z) == typeof(w) && typeof(o) == typeof(w)
                @test TS.order(z.pol) == TS.order(w.pol) && length(z.pol.coeffs) == length(w.pol.coeffs)
                @test all(iszero(mid(interval(x))) && iszero(radius(interval(x)))
                          for hp in z.pol.coeffs for x in hp.coeffs)
                @test isequal_interval(interval(o.pol.coeffs[1].coeffs[1]), interval(1.0)) &&
                      all(iszero(mid(interval(x))) && iszero(radius(interval(x)))
                          for (k, hp) in enumerate(o.pol.coeffs) for (i, x) in enumerate(hp.coeffs)
                          if !(k == 1 && i == 1))
                @test isequal_interval(z.rem, zero(Interval{Float64})) && isequal_interval(o.rem, z.rem)
                @test all(isequal_interval.(expansion_point(z), expansion_point(w))) &&
                      all(isequal_interval.(domain(z), domain(w)))
                # independent objects and arrays
                @test all(z.pol.coeffs[k] !== w.pol.coeffs[k] && z.pol.coeffs[k].coeffs !== w.pol.coeffs[k].coeffs &&
                          z.pol.coeffs[k].coeffs !== o.pol.coeffs[k].coeffs for k in eachindex(w.pol.coeffs))
                z.pol.coeffs[1].coeffs[1] = one(z.pol.coeffs[1].coeffs[1])
                @test iszero(mid(interval(o.pol.coeffs[1].coeffs[1])) - 1.0)     # `o` is not affected
            end
        end

        @testset "TMN tail: exact products up to degree m, ranges above" begin
            for m in (0, 1, typemax(Int))
                with_settings(; tail_exact = m) do
                    # TM1 products and powers go through the TMN ones
                    @test count_failures_binary(*, *, a, b, rng; name = "a * b, m = $m") == 0
                    @test count_failures_binary((x, y) -> x^3, (x, y) -> x^3, a, b, rng; name = "a³, m = $m") == 0
                    # TMN * TMN and TMN²
                    w1, w2 = deepcopy(a.pol[1]), deepcopy(b.pol[2])
                    nf = count_failures(rng; name = "TMN * TMN, m = $m") do σ, τ
                        σI = interval.(σ)
                        v1, v2 = pick(rng, w1(σI)), pick(rng, w2(σI))
                        return interval(v1) * interval(v2), (w1 * w2)(σI)
                    end
                    @test nf == 0
                    nf = count_failures(rng; name = "TMN², m = $m") do σ, τ
                        σI = interval.(σ)
                        return interval(pick(rng, w1(σI)))^2, TS.square(w1)(σI)
                    end
                    @test nf == 0
                end
            end
        end

        # Joint Cauchy accumulation (`TM._TM1_JOINT`) and explicit midpoint-radius kernels
        # (`TM._TMN_MIDRAD`, TaylorSeriesIAExt): pointwise soundness of products, squares, powers
        @testset "joint accumulation and midpoint-radius kernels" begin
            for midrad in (false, isdefined(TS, :mul_midrad)), joint in (false, true)
                with_settings(; midrad, joint) do
                    tag = "midrad = $midrad, joint = $joint"
                    @test count_failures_binary(*, *, a, b, rng; name = "a * b, $tag") == 0
                    @test count_failures_binary((x, y) -> x^2, (x, y) -> x^2, a, b, rng; name = "a², $tag") == 0
                    @test count_failures_binary((x, y) -> x^3, (x, y) -> x^3, a, b, rng; name = "a³, $tag") == 0
                    @test count_failures_binary((x, y) -> x^4, (x, y) -> x^4, a, b, rng; name = "a⁴, $tag") == 0
                end
            end
            # the joint remainders are not wider than the generic ones
            for f in (*, (x, y) -> x^2)
                rj = with_settings(() -> f(a, b); midrad = false, joint = true)
                rg = with_settings(() -> f(a, b); midrad = false, joint = false)
                @test all(diam(remainder(rj.pol[k])) <= diam(remainder(rg.pol[k])) + 1e-12
                          for k in 0:TS.order(rj.pol))
                @test diam(remainder(rj)) <= diam(remainder(rg)) + 1e-12
            end
        end

        # Accumulation of each τ^k coefficient in midpoint-radius form (`TM._TM1_ACC`), with the ranges over
        # the σ-box from the midpoint-radius data (`TM._TM1_UNIT`), against the pair-by-pair accumulation and
        # against the interval evaluation of the ranges: they enclose the same set (coefficients overlap,
        # remainders agree)
        if isdefined(TS, :midrad_poly)
            @testset "midpoint-radius accumulation of the τ coefficients" begin
                function agree(r1, r2)
                    overlap = true
                    for k in 0:TS.order(r1.pol), q in eachindex(r1.pol[k].pol.coeffs),
                            h in eachindex(r1.pol[k].pol.coeffs[q].coeffs)
                        overlap &= !isdisjoint_interval(r1.pol[k].pol.coeffs[q].coeffs[h],
                                                        r2.pol[k].pol.coeffs[q].coeffs[h])
                    end
                    return overlap &&
                        all(isapprox(diam(remainder(r1.pol[k])), diam(remainder(r2.pol[k]));
                                     rtol = 1e-6, atol = 1e-12) for k in 0:TS.order(r1.pol)) &&
                        isapprox(diam(remainder(r1)), diam(remainder(r2)); rtol = 1e-6, atol = 1e-12)
                end
                for f in (*, (x, y) -> x^2, (x, y) -> x^3)
                    racc = with_settings(() -> f(a, b); midrad = true, joint = true, acc = true, unit = true)
                    rprs = with_settings(() -> f(a, b); midrad = true, joint = true, acc = false)
                    rias = with_settings(() -> f(a, b); midrad = true, joint = true, acc = true, unit = false)
                    @test agree(racc, rprs)          # accumulation against pair by pair
                    @test agree(racc, rias)          # ranges from midpoint-radius data against interval evaluation
                end
                with_settings(; midrad = true, joint = true, acc = true, unit = true) do
                    @test count_failures_binary(*, *, a, b, rng; name = "a * b, accumulation") == 0
                    @test count_failures_binary((x, y) -> x^2, (x, y) -> x^2, a, b, rng; name = "a², accumulation") == 0
                    @test count_failures_binary((x, y) -> x^3, (x, y) -> x^3, a, b, rng; name = "a³, accumulation") == 0
                end
            end
        end

        # The ranges of the coefficients computed by the joint accumulation (polynomial part + remainder) are
        # reused in the τ-tail of products and squares: they agree with the evaluation over the σ-box
        @testset "τ-tail: ranges of the coefficients from the joint computation" begin
            auxN = TM.centered_dom(a[0])
            eva = [a.pol[i](auxN) for i in 0:orderT]
            evb = [b.pol[i](auxN) for i in 0:orderT]
            near(x, y) = !isdisjoint_interval(x, y) && isapprox(diam(x), diam(y); rtol = 1e-6, atol = 1e-12)
            for midrad in (false, isdefined(TS, :mul_midrad)), unit in (false, true)
                with_settings(; midrad, joint = true, unit) do
                    _, ra, rb = TM._cauchy_joint(a, b, false)
                    @test all(near(ra[i+1], eva[i+1]) && near(rb[i+1], evb[i+1]) for i in 0:orderT)
                    # remainders with the reused ranges against the ones that evaluate them
                    @test isapprox(diam(remainder(a * b)), diam(TM.remainder_product(a, b, orderT));
                                   rtol = 1e-6, atol = 1e-12)
                    @test isapprox(diam(remainder(TS.square(a))), diam(TM.remainder_square(a, orderT));
                                   rtol = 1e-6, atol = 1e-12)
                end
            end
        end

        @testset "rpa functions (exp, cos)" begin
            a2 = random_tm1!(c.x2NI[1], rng; scale = 0.2, remc = 1e-6, remt = 1e-5)
            @test count_failures_binary((x, y) -> exp(x), (x, y) -> exp(x), a2, b, rng; name = "exp") == 0
            @test count_failures_binary((x, y) -> cos(x), (x, y) -> cos(x), a2, b, rng; name = "cos") == 0
        end
    end
end