# Tests using TMAbsRem and TMRelRem

using TaylorModels

if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
    eeuler = Base.e
else
    using Test
    eeuler = Base.MathConstants.e
end

@testset "Tests for TMAbsRem" begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)
    x1 = Interval(1.0)
    ii1 = Interval(0.5, 1.5)

    @testset "TMAbsRem constructors" begin
        tv = TMAbsRem{Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TMAbsRem(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TMAbsRem(5, x0, ii0)
        @test TMAbsRem(x1, 5, x0, ii0) == TMAbsRem(Taylor1(x1, 5), x0, x0, ii0)

        # Test errors in construction
        @test_throws AssertionError TMAbsRem(Taylor1(Interval{Float64},5), x1, x0, ii0)
        @test_throws AssertionError TMAbsRem(5, x1, ii0)

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
    end

    @testset "TMAbsRem bounds" begin
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

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), TMAbsRem(4, x0, ii0)) ==
            TMAbsRem(interval(5.0), 4, x0, ii0)
        @test rpa(x->5*x, TMAbsRem(4, x1, ii1)) ==
            TMAbsRem(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^4, TMAbsRem(4, x0, ii0)) ==
            TMAbsRem( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, TMAbsRem(3, x0, ii0)) ==
            TMAbsRem( Taylor1(x0, 3), 5*ii0^4, x0, ii0)

        # Testing remainders of an RPA
        ftest = x -> exp(x)-1
        tma = rpa(ftest, TMAbsRem(2, x0, ii0))
        tmb = ftest(TMAbsRem(2, x0, ii0))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii0.lo)-fT(ii0.lo-ξ0),
                        ftest(ii0.hi)-fT(ii0.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii0)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> exp(x)
        tma = rpa(ftest, TMAbsRem(2, x1, ii1))
        tmb = ftest(TMAbsRem(2, x1, ii1))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii1.lo)-fT(ii1.lo-ξ0),
                        ftest(ii1.hi)-fT(ii1.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii1)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> sin(x)
        tma = rpa(ftest, TMAbsRem(3, x0, ii0))
        tmb = ftest(TMAbsRem(3, x0, ii0))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii0.lo)-fT(ii0.lo-ξ0),
                        ftest(ii0.hi)-fT(ii0.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii0)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> sqrt(x)
        tma = rpa(ftest, TMAbsRem(2, x1, ii1))
        tmb = ftest(TMAbsRem(2, x1, ii1))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii1.lo)-fT(ii1.lo-ξ0),
                        ftest(ii1.hi)-fT(ii1.hi-ξ0)) ⊆ remainder(tma)
        x = radius(ii1)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)

        ftest = x -> inv(x)
        tma = rpa(ftest, TMAbsRem(5, x1, ii1))
        tmb = ftest(TMAbsRem(5, x1, ii1))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii1.hi)-fT(ii1.hi-ξ0),
                        ftest(ii1.lo)-fT(ii1.lo-ξ0)) ⊆ remainder(tma)
        x = radius(ii1)*rand() + ξ0/2
        @test fT(x-ξ0)+sup(Δ) ≥ ftest(x)
        @test ftest(x) ≥ fT(x-ξ0)+inf(Δ)
    end

end
