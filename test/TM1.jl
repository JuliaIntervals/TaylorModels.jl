# Tests using TM1AbsRem and TM1RelRem

using TaylorModels
using TaylorSeries, IntervalArithmetic

const _num_tests = 1000

if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
    eeuler = Base.e
else
    using Test
    eeuler = Base.MathConstants.e
end

function check_inclusion(ftest, tma::T) where {T<:Union{TM1AbsRem, TM1RelRem}}
    ii = tma.iI
    xfp = diam(tma.iI)*(rand()-0.5) + mid(tma.x0)
    xbf = big(xfp)
    range = tma(@interval(xfp)-tma.x0)
    bb = ftest(xbf) ∈ range
    bb || @show(ftest, xfp, xbf, ftest(xbf), range)
    return bb
end

@testset "Test `bound_taylor1`" begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)

    tpol = exp( Taylor1(2) )
    @test TaylorModels.bound_taylor1( tpol, ii0) ==
        @interval(tpol(ii0.lo), tpol(ii0.hi))
    @test TaylorModels.bound_taylor1( exp( Taylor1(Interval{Float64}, 2) ),
        ii0) == @interval(tpol(ii0.lo), tpol(ii0.hi))

    # An uncomfortable example from Makino
    t = Taylor1(5)
    @test interval(1-4^4/5^5,1) ⊆ TaylorModels.bound_taylor1(1-t^4+t^5, 0..1)
end

@testset "Tests for TM1AbsRem " begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)
    x1 = Interval(1.0)
    ii1 = Interval(0.5, 1.5)

    @testset "TM1AbsRem constructors" begin
        tv = TM1AbsRem{Interval{Float64},Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TM1AbsRem(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TM1AbsRem(5, x0, ii0)
        @test TM1AbsRem(x1, 5, x0, ii0) == TM1AbsRem(Taylor1(x1, 5), x0, x0, ii0)

        @test isa(tv, AbstractSeries)
        @test TM1AbsRem{Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Test errors in construction
        @test_throws AssertionError TM1AbsRem(Taylor1(Interval{Float64},5), x1, x0, ii0)
        @test_throws AssertionError TM1AbsRem(5, x1, ii0)

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = TM1AbsRem(x1+Taylor1(5), Δ, x1, ii1)
        tv = TM1AbsRem(5, x1, ii1)
        @test zero(a) == TM1AbsRem(zero(a.pol), 0..0, x1, ii1)
        @test one(a) == TM1AbsRem(one(a.pol), 0..0, x1, ii1)
        @test a+x1 == TM1AbsRem(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == TM1AbsRem(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == TM1AbsRem(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == TM1AbsRem(zero(a.pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == TM1AbsRem(a.pol*tv.pol, a.rem*tv.pol(ii1-x1), x1, ii1)
        @test b/tv == TM1AbsRem(a.pol, Interval(-0.78125, 0.84375), x1, ii1)
        b = a * a.pol[0]
        @test b == a

        a = TM1AbsRem(x0, 5, x0, ii0)
        @test a^2 == TM1AbsRem(x0^2, 5, x0, ii0)
        @test a^3 == TM1AbsRem(x0^3, 5, x0, ii0)
        a = TM1AbsRem(x1, 5, x1, ii1)
        @test a^2 == TM1AbsRem(x1^2, 5, x1, ii1)
        @test a^3 == TM1AbsRem(x1^3, 5, x1, ii1)
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
        order = 2
        ii = ii0
        xx = x0
        ftest = x -> exp(x)-1
        tma = rpa(ftest, TM1AbsRem(order, xx, ii))
        tmb = ftest(TM1AbsRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> exp(x)
        tma = rpa(ftest, TM1AbsRem(order, xx, ii))
        tmb = ftest(TM1AbsRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 3
        ii = ii0
        xx = x0
        ftest = x -> sin(x)
        tma = rpa(ftest, TM1AbsRem(order, xx, ii))
        tmb = ftest(TM1AbsRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> sqrt(x)
        tma = rpa(ftest, TM1AbsRem(order, xx, ii))
        tmb = ftest(TM1AbsRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 5
        ii = ii1
        xx = x1
        ftest = x -> inv(x)
        tma = rpa(ftest, TM1AbsRem(order, xx, ii))
        tmb = ftest(TM1AbsRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0 = rpafp(tma)
        @test interval(ftest(ii.hi)-fT(ii.hi-ξ0),
                        ftest(ii.lo)-fT(ii.lo-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tma = rpa(ftest, TM1AbsRem(order, xx, ii))
        tmb = ftest(TM1AbsRem(order, xx, ii))
        @test remainder(tmb) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end
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

    @testset "Tests for integrate" begin
        order = 4
        tm = TM1AbsRem(order, x0, ii0)

        integ_res = integrate(exp(tm), 1..1)
        exact_res = exp(tm)
        @test exact_res.pol ⊆ integ_res.pol
        @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_inclusion(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol ⊆ integ_res.pol
        @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_inclusion(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol ⊆ integ_res.pol
        @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_inclusion(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol ⊆ integ_res.pol
        # @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_inclusion(atan, integ_res)
        end
    end
end

@testset "Tests for TM1RelRem " begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)
    x1 = Interval(1.0)
    ii1 = Interval(0.5, 1.5)

    @testset "TM1RelRem constructors" begin
        tv = TM1RelRem{Interval{Float64},Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TM1RelRem(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TM1RelRem(5, x0, ii0)
        @test TM1RelRem(x1, 5, x0, ii0) == TM1RelRem(Taylor1(x1, 5), x0, x0, ii0)

        @test isa(tv, AbstractSeries)
        @test TM1RelRem{Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Zero may not be contained in the remainder of a TM1RelRem
        @test 0 ∉ remainder(TM1RelRem(Taylor1(Interval{Float64},5), x1, x0, ii0))

        # Test errors in construction
        @test_throws AssertionError TM1RelRem(5, x1, ii0)

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = TM1RelRem(x1+Taylor1(5), Δ, x1, ii1)
        tv = TM1RelRem(5, x1, ii1)

        @test zero(a) == TM1RelRem(zero(a.pol), 0..0, x1, ii1)
        @test one(a) == TM1RelRem(one(a.pol), 0..0, x1, ii1)
        @test a+x1 == TM1RelRem(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == TM1RelRem(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == TM1RelRem(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == TM1RelRem(zero(a.pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == TM1RelRem(a.pol*tv.pol, a.rem*tv.pol(ii1-x1), x1, ii1)
        @test b/tv == TM1RelRem(a.pol, Interval(-2.75, 4.75), x1, ii1)
        b = a * a.pol[0]
        @test b == a

        a = TM1RelRem(x0, 5, x0, ii0)
        @test a^2 == TM1RelRem(x0^2, 5, x0, ii0)
        @test a^3 == TM1RelRem(x0^3, 5, x0, ii0)
        a = TM1RelRem(x1, 5, x1, ii1)
        @test a^2 == TM1RelRem(x1^2, 5, x1, ii1)
        @test a^3 == TM1RelRem(x1^3, 5, x1, ii1)
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), TM1RelRem(4, x0, ii0)) ==
            TM1RelRem(interval(5.0), 4, x0, ii0)
        @test rpa(x->5*x, TM1RelRem(4, x1, ii1)) ==
            TM1RelRem(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^4, TM1RelRem(4, x0, ii0)) ==
            TM1RelRem( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, TM1RelRem(3, x0, ii0)) ==
            TM1RelRem( Taylor1(x0, 3), interval(5), x0, ii0)

        # Testing remainders and inclusion of RPAs
        order = 2
        ii = ii0
        xx = x0
        ftest = x -> exp(x)-1
        tma = rpa(ftest, TM1RelRem(order, xx, ii))
        tmb = ftest(TM1RelRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0, δ = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> exp(x)
        tma = rpa(ftest, TM1RelRem(order, xx, ii))
        tmb = ftest(TM1RelRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0, δ = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 3
        ii = ii0
        xx = x0
        ftest = x -> sin(x)
        tma = rpa(ftest, TM1RelRem(order, xx, ii))
        tmb = ftest(TM1RelRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0, δ = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> sqrt(x)
        tma = rpa(ftest, TM1RelRem(order, xx, ii))
        tmb = ftest(TM1RelRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0, δ = rpafp(tma)
        @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
                        ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        order = 5
        ii = ii1
        xx = x1
        ftest = x -> inv(x)
        tma = rpa(ftest, TM1RelRem(order, xx, ii))
        tmb = ftest(TM1RelRem(order, xx, ii))
        @test tma == tmb
        fT, Δ, ξ0, δ = rpafp(tma)
        @test interval(ftest(ii.hi)-fT(ii.hi-ξ0),
                        ftest(ii.lo)-fT(ii.lo-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tma = rpa(ftest, TM1RelRem(order, xx, ii))
        tmb = ftest(TM1RelRem(order, xx, ii))
        @test remainder(tmb) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_inclusion(ftest, tma)
        end
    end

    @testset "Composition of functions and their inverses" begin
        tv = TM1RelRem(2, x0, ii0)

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
        tv = TM1RelRem(2, x1, ii1)

        tma = log(tv)
        tmb = exp(tma)
        @test tmb == exp(log(tv))
        @test tmb.pol == tv.pol

        tma = cos(tv)
        tmb = acos(tma)
        @test tmb == acos(cos(tv))
        @test sup(norm(tmb.pol - tv.pol, Inf)) < 1.0e-15
    end

    @testset "Tests for integrate" begin
        order = 4
        tm = TM1RelRem(order, x0, ii0)

        integ_res = integrate(exp(tm), 1..1)
        exact_res = exp(tm)
        @test exact_res.pol ⊆ integ_res.pol
        @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol ⊆ integ_res.pol
        @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol ⊆ integ_res.pol
        @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol ⊆ integ_res.pol
        # @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_inclusion(atan, integ_res)
        end
    end
end
