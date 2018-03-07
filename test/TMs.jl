# Tests using TM1AbsRem and TM1RelRem

using TaylorModels

if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
    eeuler = Base.e
else
    using Test
    eeuler = Base.MathConstants.e
end

@testset "Tests for TM1AbsRem" begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)
    x1 = Interval(1.0)
    ii1 = Interval(0.5, 1.5)

    @testset "TM1AbsRem constructors" begin
        tv = TM1AbsRem{Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TM1AbsRem(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TM1AbsRem(5, x0, ii0)
        @test TM1AbsRem(x1, 5, x0, ii0) == TM1AbsRem(Taylor1(x1, 5), x0, x0, ii0)

        # Test errors in construction
        @test_throws AssertionError TM1AbsRem(Taylor1(Interval{Float64},5), x1, x0, ii0)
        @test_throws AssertionError TM1AbsRem(5, x1, ii0)

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
    end

    @testset "TM1AbsRem bounds" begin
        tpol = exp( Taylor1(2) )
        @test TaylorModels.bound_taylor1( tpol, ii0) ==
            @interval(tpol(ii0.lo), tpol(ii0.hi))
        @test TaylorModels.bound_taylor1( exp( Taylor1(Interval{Float64}, 2) ),
            ii0) == @interval(tpol(ii0.lo), tpol(ii0.hi))

        # An uncomfortable example from Makino
        t = Taylor1(5)
        @test TaylorModels.bound_taylor1(1.0-t^4+t^5, interval(0,1)) ==
            Interval(0.9180799999999999, 1.0)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = TM1AbsRem(x1+Taylor1(5), Δ, x1, ii1)
        tv = TM1AbsRem(5, x1, ii1)
        @test a+x1 == TM1AbsRem(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == TM1AbsRem(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == TM1AbsRem(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == TM1AbsRem(zero(a.pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == TM1AbsRem(a.pol*tv.pol, a.rem*tv.pol(ii1-x1), x1, ii1)
        @test b/tv == TM1AbsRem(a.pol, Interval(-0.78125, 0.84375), x1, ii1)
        b = a * a.pol[0]
        @test b == a
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), TM1AbsRem(4, x0, ii0)) ==
            TM1AbsRem(interval(5.0), 4, x0, ii0)
        @test rpa(x->5*x, TM1AbsRem(4, x1, ii1)) ==
            TM1AbsRem(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^4, TM1AbsRem(4, x0, ii0)) ==
            TM1AbsRem( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, TM1AbsRem(3, x0, ii0)) ==
            TM1AbsRem( Taylor1(x0, 3), 5*ii0^4, x0, ii0)

        # Testing remainders of an RPA
        ftest = x -> exp(x)-1
        tma = rpa(ftest, TM1AbsRem(2, x0, ii0))
        tmb = ftest(TM1AbsRem(2, x0, ii0))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii0.lo)-fT(ii0.lo-ξ0),
                        ftest(ii0.hi)-fT(ii0.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii0)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> exp(x)
        tma = rpa(ftest, TM1AbsRem(2, x1, ii1))
        tmb = ftest(TM1AbsRem(2, x1, ii1))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii1.lo)-fT(ii1.lo-ξ0),
                        ftest(ii1.hi)-fT(ii1.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii1)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> sin(x)
        tma = rpa(ftest, TM1AbsRem(3, x0, ii0))
        tmb = ftest(TM1AbsRem(3, x0, ii0))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii0.lo)-fT(ii0.lo-ξ0),
                        ftest(ii0.hi)-fT(ii0.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii0)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> sqrt(x)
        tma = rpa(ftest, TM1AbsRem(2, x1, ii1))
        tmb = ftest(TM1AbsRem(2, x1, ii1))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii1.lo)-fT(ii1.lo-ξ0),
                        ftest(ii1.hi)-fT(ii1.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii1)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> inv(x)
        tma = rpa(ftest, TM1AbsRem(5, x1, ii1))
        tmb = ftest(TM1AbsRem(5, x1, ii1))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii1.hi)-fT(ii1.hi-ξ0),
                        ftest(ii1.lo)-fT(ii1.lo-ξ0)) ⊆ remainder(tma)
        x = radius(ii1)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)
    end

    @testset "Composition of functions and their inverses" begin
        tv = TM1AbsRem(2, x0, ii0)

        tma = exp(tv)
        tmb = log(tma)
        @test tmb == log(exp(tv))
        @test tmb.pol == tv.pol

        tma = sin(tv)
        tmb = asin(tma)
        @test tmb == asin(sin(tv))
        @test tmb.pol == tv.pol

        tma = asin(tv)
        tmb = sin(tma)
        @test tmb == sin(asin(tv))
        @test tmb.pol == tv.pol

        tma = acos(tv)
        tmb = cos(tma)
        @test tmb == cos(acos(tv))
        @test sup(norm(tmb.pol - tv.pol, Inf)) < 1.0e-16

        tma = tan(tv)
        tmb = atan(tma)
        @test tmb == atan(tan(tv))
        @test tmb.pol == tv.pol

        tma = atan(tv)
        tmb = tan(tma)
        @test tmb == tan(atan(tv))
        @test tmb.pol == tv.pol


        ####
        tv = TM1AbsRem(2, x1, ii1)

        tma = log(tv)
        tmb = exp(tma)
        @test tmb == exp(log(tv))
        @test tmb.pol == tv.pol

        tma = cos(tv)
        tmb = acos(tma)
        @test tmb == acos(cos(tv))
        @test sup(norm(tmb.pol - tv.pol, Inf)) < 1.0e-15
    end
end
