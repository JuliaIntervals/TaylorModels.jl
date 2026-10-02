# Regression tests for picard_lindelof / picard_lindelof! (validated_integ3).
#
#  1. No aliasing: entries of the result (and of the cache vectors) are
#     independent objects, down to the TaylorModelN coefficients.
#  2. Independence under mutation: changing x2N[1] leaves x2N[2] untouched.
#  3. Correctness: for x1N constant in τ (x1N = x0), the Picard iterate is known:
#       autonomous      x' = (x₂, -1)   →  P = (x0₁ + x0₂ τ,  x0₂ - τ)
#       non-autonomous  x' = (t, x₁)    →  P = (x0₁ + t0 τ + τ²/2,  x0₂ + x0₁ τ)
#  4. Allocating and in-place versions agree, and the input is not modified.

using TaylorModels, TaylorModels.ValidatedInteg
using Test, StaticArrays
const VI = TaylorModels.ValidatedInteg

const orderQ, orderT = 2, 4
const sp = JetSpace(2*orderQ, ["ξ₁", "ξ₂"])
const X0 = SVector(interval(9.75, 10.25), interval(-0.25, 0.25))

function fb!(dx, x, p, t)                  # autonomous
    dx[1] = x[2]
    dx[2] = -one(x[1])
    nothing
end
function nonaut!(dx, x, p, t)              # non-autonomous
    dx[1] = t + zero(x[1])                 # new object, do not alias t
    dx[2] = x[1]
    nothing
end

# Build x1N (= x0, constant in τ), dx1N, t1N (= t0 + τ) as validated_integ3 does
function setup(f!, t0, δt)
    c = VI.init_cache_VI3(f!, t0, X0, 2, orderT, orderQ, sp, nothing; parse_eqs = false)
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

# true iff no two entries share an object, at any level down to TMN coefficients
function no_alias(v)
    for i in eachindex(v), j in eachindex(v)
        i == j && continue
        v[i] === v[j] && return false
        v[i].pol.coeffs === v[j].pol.coeffs && return false
        for k in eachindex(v[i].pol.coeffs)
            v[i].pol.coeffs[k] === v[j].pol.coeffs[k] && return false
            v[i].pol.coeffs[k].pol.coeffs === v[j].pol.coeffs[k].pol.coeffs && return false
        end
    end
    return true
end
# no entry of a is (or shares coefficients with) an entry of b
disjoint(a, b) = all(a[i] !== b[j] && a[i].pol.coeffs !== b[j].pol.coeffs
                     for i in eachindex(a), j in eachindex(b))

cpoly(tm1, k) = polynomial(tm1[k])                      # TaylorN of τ^k coefficient
constpoly(x, ref) = x + zero(ref)                       # constant TaylorN like ref

@testset "picard_lindelof" begin

    @testset "cache vectors are not aliased" begin
        c, x1N, dx1N, t1N = setup(fb!, 0.0, 0.5)
        @test no_alias(c.x1N) && no_alias(c.dx1N) && no_alias(c.x2N)
        @test disjoint(c.x1N, c.dx1N) && disjoint(c.x1N, c.x2N) && disjoint(c.dx1N, c.x2N)
    end

    @testset "no aliasing and independence (allocating version)" begin
        _, x1N, dx1N, t1N = setup(fb!, 0.0, 0.5)
        x2N = picard_lindelof(fb!, dx1N, x1N, t1N, nothing)
        @test no_alias(x2N)
        @test disjoint(x2N, x1N) && disjoint(x2N, dx1N)
        before = deepcopy(cpoly(x2N[2], 0))
        x2N[1].pol.coeffs[1].pol.coeffs[1].coeffs[1] += 1.0       # mutate entry 1
        @test cpoly(x2N[2], 0) == before                          # entry 2 unchanged
    end

    @testset "autonomous: x' = (x₂, -1)" begin
        _, x1N, dx1N, t1N = setup(fb!, 0.0, 0.5)
        x10, x20 = cpoly(x1N[1], 0), cpoly(x1N[2], 0)
        x2N = picard_lindelof(fb!, dx1N, x1N, t1N, nothing)
        @test cpoly(x2N[1], 0) == x10
        @test cpoly(x2N[1], 1) == x20
        @test cpoly(x2N[2], 0) == x20
        @test cpoly(x2N[2], 1) == constpoly(-1.0, x10)
        @test all(iszero(cpoly(x2N[i], k)) for i in 1:2 for k in 2:orderT)
    end

    @testset "non-autonomous: x' = (t, x₁), t0 ≠ 0" begin
        t0 = 2.0
        _, x1N, dx1N, t1N = setup(nonaut!, t0, 0.5)
        x10, x20 = cpoly(x1N[1], 0), cpoly(x1N[2], 0)
        x2N = picard_lindelof(nonaut!, dx1N, x1N, t1N, nothing)   # errors if args swapped
        @test cpoly(x2N[1], 0) == x10
        @test cpoly(x2N[1], 1) == constpoly(t0, x10)
        @test cpoly(x2N[1], 2) == constpoly(0.5, x10)
        @test cpoly(x2N[2], 0) == x20
        @test cpoly(x2N[2], 1) == x10
        @test all(iszero(cpoly(x2N[2], k)) for k in 2:orderT)
    end

    @testset "allocating == in-place; input not modified" begin
        _, x1N, dx1N, t1N = setup(nonaut!, 2.0, 0.5)
        x1N_copy = deepcopy(x1N)
        x2N = picard_lindelof(nonaut!, dx1N, x1N, t1N, nothing)
        y2N = zero.(x1N)
        VI.picard_lindelof!(nonaut!, y2N, dx1N, x1N, t1N, nothing)
        for i in eachindex(x2N), k in 0:orderT
            @test cpoly(x2N[i], k) == cpoly(y2N[i], k)
        end
        @test all(isequal_interval(remainder(x2N[i]), remainder(y2N[i])) for i in eachindex(x2N))
        @test all(cpoly(x1N[i], k) == cpoly(x1N_copy[i], k) for i in eachindex(x1N), k in 0:orderT)
    end
end