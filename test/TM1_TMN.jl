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

module TM1TMNTests

using TaylorModels, TaylorModels.ValidatedInteg
using Test, Random, StaticArrays

const TM = TaylorModels
const VI = TaylorModels.ValidatedInteg
const orderQ, orderT = 2, 4
const H = 1.0                                   # τ ∈ [0, H]
const NS = 40                                   # samples per check

function with_space(f, sp)
    old = TS.default_space[]
    TS.set_default_space!(sp; nowarn=true)
    try
        return f(TS.default_space[])
    finally
        TS.set_default_space!(old; nowarn=true)
    end
end

# fill `tm` (in place) with random interval coefficients, TMN remainders and TM1 remainder
function random_tm1!(tm, rng; scale = 1.0, remc = 1e-4, remt = 1e-3)
    tm.dom = interval(0.0, H); tm.x0 = 0.0
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
function samples(rng, n)
    out = Tuple{Vector{Float64},Float64}[]
    for k in 1:n
        τ = k == 1 ? 0.0 : (k == 2 ? H : H * rand(rng))
        push!(out, (2 .* rand(rng, 2) .- 1, τ))
    end
    return out
end

# number of failures of  truth(f_a, f_b) ⊆ model(a, b)(τ, σ)
function nfail_binary(model, truth, a, b, rng; name = "")
    res = model(a, b)
    nf = 0
    for (σ, τ) in samples(rng, NS)
        τI = interval(τ)
        sa, sb = draw(a, σ, rng), draw(b, σ, rng)
        t, m = truth(fval(sa, τI), fval(sb, τI)), model_val(res, τ, σ)
        if !issubset_interval(t, m)
            nf += 1
            nf == 1 && @warn "pointwise failure" name σ τ truth = t model = m
        end
    end
    return nf
end

@testset "TaylorModel1{TaylorModelN{Interval}} arithmetic" begin
    X0 = SVector(interval(1.0, 1.1), interval(-0.25, 0.25))
    function f!(dx, x, p, t)
        dx[1] = x[2]
        dx[2] = -one(x[1])
        nothing
    end
    with_space(JetSpace(2*orderQ, ["ξ₁", "ξ₂"])) do sp
        rng = MersenneTwister(2024)
        c = VI.init_cache_VI3(f!, 0.0, X0, 2, orderT, orderQ, sp, nothing; parse_eqs = true)
        a = random_tm1!(c.x1NI[1], rng)
        b = random_tm1!(c.x1NI[2], rng)

        @testset "TM1 ± TM1, TM1 * TM1, square, powers" begin
            @test nfail_binary(+, +, a, b, rng; name = "a + b") == 0
            @test nfail_binary(-, -, a, b, rng; name = "a - b") == 0
            @test nfail_binary(*, *, a, b, rng; name = "a * b") == 0
            @test nfail_binary((x, y) -> TS.square(x), (x, y) -> x^2, a, b, rng; name = "a²") == 0
            @test nfail_binary((x, y) -> x^3, (x, y) -> x * x * x, a, b, rng; name = "a³") == 0
            @test nfail_binary((x, y) -> x^4, (x, y) -> (x * x) * (x * x), a, b, rng; name = "a⁴") == 0
        end

        @testset "TM1 ± number, number ± TM1, TM1 / number" begin
            for num in (0.7, 2, -1.3)
                nI = interval(num)
                @test nfail_binary((x, y) -> x + num, (x, y) -> x + nI, a, b, rng; name = "a + $num") == 0
                @test nfail_binary((x, y) -> x - num, (x, y) -> x - nI, a, b, rng; name = "a - $num") == 0
                @test nfail_binary((x, y) -> num + x, (x, y) -> nI + x, a, b, rng; name = "$num + a") == 0
                @test nfail_binary((x, y) -> num - x, (x, y) -> nI - x, a, b, rng; name = "$num - a") == 0
                @test nfail_binary((x, y) -> x / num, (x, y) -> x / nI, a, b, rng; name = "a / $num") == 0
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
            nf = 0
            for (σ, τ) in samples(rng, NS)
                τI = interval(τ)
                sa = draw(a, σ, rng)
                wv = pick(rng, w(interval.(σ)))
                t = fval(sa, τI) * interval(wv)
                nf += !issubset_interval(t, model_val(res, τ, σ))
            end
            @test nf == 0
        end

        @testset "integrate! with integration constant" begin
            c0 = deepcopy(b.pol[0])
            c0.rem = interval(-1e-3, 2e-3)               # c0.rem must be kept
            res = deepcopy(a)
            TM.integrate!(res, a, c0, domain(a[0]))
            nf = 0
            for (σ, τ) in samples(rng, NS)
                τI = interval(τ)
                sa = draw(a, σ, rng)
                wv = pick(rng, c0(interval.(σ)))
                t = interval(wv) + interval(sa.r) * τI +
                    sum(interval(sa.v[i+1]) * τI^(i+1) / interval(i+1) for i in 0:orderT)
                nf += !issubset_interval(t, model_val(res, τ, σ))
            end
            @test nf == 0
        end

        @testset "TMN tail: exact products up to degree m, ranges above" begin
            old = TM._TMN_TAIL_EXACT[]
            try
                for m in (0, 1, typemax(Int))
                    TM._TMN_TAIL_EXACT[] = m
                    # TM1 products and squares go through the TMN ones
                    @test nfail_binary(*, *, a, b, rng; name = "a * b, m = $m") == 0
                    @test nfail_binary((x, y) -> x^3, (x, y) -> x^3, a, b, rng; name = "a³, m = $m") == 0
                    # TMN * TMN and TMN²
                    w1, w2 = deepcopy(a.pol[1]), deepcopy(b.pol[2])
                    nf = 0
                    for (σ, τ) in samples(rng, NS)
                        σI = interval.(σ)
                        v1, v2 = pick(rng, w1(σI)), pick(rng, w2(σI))
                        nf += !issubset_interval(interval(v1) * interval(v2), (w1 * w2)(σI))
                        nf += !issubset_interval(interval(v1)^2, TS.square(w1)(σI))
                    end
                    @test nf == 0
                end
            finally
                TM._TMN_TAIL_EXACT[] = old
            end
        end

        # Joint Cauchy accumulation (`TM._TM1_JOINT`) and explicit midpoint-radius kernels
        # (`TM._TMN_MIDRAD`, TaylorSeriesIAExt): pointwise soundness of products, squares, powers
        @testset "joint accumulation and midpoint-radius kernels" begin
            old_m, old_j = TM._TMN_MIDRAD[], TM._TM1_JOINT[]
            try
                for midrad in (false, isdefined(TS, :mul_midrad)), joint in (false, true)
                    TM._TMN_MIDRAD[], TM._TM1_JOINT[] = midrad, joint
                    tag = "midrad = $midrad, joint = $joint"
                    @test nfail_binary(*, *, a, b, rng; name = "a * b, $tag") == 0
                    @test nfail_binary((x, y) -> x^2, (x, y) -> x^2, a, b, rng; name = "a², $tag") == 0
                    @test nfail_binary((x, y) -> x^3, (x, y) -> x^3, a, b, rng; name = "a³, $tag") == 0
                    @test nfail_binary((x, y) -> x^4, (x, y) -> x^4, a, b, rng; name = "a⁴, $tag") == 0
                end
                # the joint remainders are not wider than the generic ones
                TM._TMN_MIDRAD[] = false
                for (f, nm) in ((*, "a*b"), ((x, y) -> x^2, "a²"))
                    TM._TM1_JOINT[] = true;  rj = f(a, b)
                    TM._TM1_JOINT[] = false; rg = f(a, b)
                    @test all(diam(remainder(rj.pol[k])) <= diam(remainder(rg.pol[k])) + 1e-12
                              for k in 0:TS.order(rj.pol))
                    @test diam(remainder(rj)) <= diam(remainder(rg)) + 1e-12
                end
            finally
                TM._TMN_MIDRAD[], TM._TM1_JOINT[] = old_m, old_j
            end
        end

        @testset "rpa functions (exp, cos)" begin
            a2 = random_tm1!(c.x2NI[1], rng; scale = 0.2, remc = 1e-6, remt = 1e-5)
            @test nfail_binary((x, y) -> exp(x), (x, y) -> exp(x), a2, b, rng; name = "exp") == 0
            @test nfail_binary((x, y) -> cos(x), (x, y) -> cos(x), a2, b, rng; name = "cos") == 0
        end
    end
end

end # module