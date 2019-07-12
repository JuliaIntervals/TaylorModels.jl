# Tests using TaylorModel1 and RTaylorModel1

using TaylorModels
using LinearAlgebra: norm
using Test

const _num_tests = 1000
const α_mid = TaylorModels.α_mid

setformat(:full)


function check_containment(ftest, tma::T) where {T<:Union{TaylorModel1, RTaylorModel1}}
    ii = tma.dom
    xfp = diam(tma.dom)*(rand()-0.5) + mid(tma.x0)
    xbf = big(xfp)
    range = tma(@interval(xfp)-tma.x0)
    bb = ftest(xbf) ∈ range
    bb || @show(ftest, xfp, xbf, ftest(xbf), range)
    return bb
end

# @testset "Test `bound_taylor1`" begin
#     x0 = Interval(0.0)
#     ii0 = Interval(-0.5, 0.5)
#
#     tpol = exp( Taylor1(2) )
#     @test TaylorModels.bound_taylor1( tpol, ii0) ==
#         @interval(tpol(ii0.lo), tpol(ii0.hi))
#     @test TaylorModels.bound_taylor1( exp( Taylor1(Interval{Float64}, 2) ),
#         ii0) == @interval(tpol(ii0.lo), tpol(ii0.hi))
#
#     # An uncomfortable example from Makino
#     t = Taylor1(5)
#     @test interval(1-4^4/5^5,1) ⊆ TaylorModels.bound_taylor1(1-t^4+t^5, 0..1)
# end

@testset "Tests for TaylorModel1 " begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)
    x1 = Interval(1.0)
    ii1 = Interval(0.5, 1.5)

    @testset "TaylorModel1 constructors" begin
        tv = TaylorModel1{Interval{Float64},Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TaylorModel1(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TaylorModel1(5, x0, ii0)
        @test tv == TaylorModel1(5, ii0)
        @test tv == TaylorModel1(5, 0.0, ii0)
        @test TaylorModel1(x1, 5, x0, ii0) == TaylorModel1(Taylor1(x1, 5), x0, x0, ii0)
        @test TaylorModel1(5, 0.7, ii1) == TaylorModel1(5, interval(0.7), ii1)

        @test isa(tv, AbstractSeries)
        @test TaylorModel1{Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Test errors in construction
        @test_throws AssertionError TaylorModel1(Taylor1(Interval{Float64},5), x1, x0, ii0)
        @test_throws AssertionError TaylorModel1(5, x1, ii0)

        # Tests for get_order, remainder, polynomial and domain
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
        @test polynomial(tv) == Taylor1(Interval{Float64},5)
        @test domain(tv) == ii0
        @test constant_term(tv) == interval(0.0)
        @test linear_polynomial(tv) == Taylor1(Interval{Float64},5)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = TaylorModel1(x1+Taylor1(5), Δ, x1, ii1)
        tv = TaylorModel1(5, x1, ii1)
        @test zero(a) == TaylorModel1(zero(a.pol), 0..0, x1, ii1)
        @test one(a) == TaylorModel1(one(a.pol), 0..0, x1, ii1)
        @test a+x1 == TaylorModel1(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == TaylorModel1(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == TaylorModel1(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == TaylorModel1(zero(a.pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == TaylorModel1(a.pol*tv.pol, a.rem*tv.pol(ii1-x1), x1, ii1)
        @test b/tv == TaylorModel1(a.pol, Interval(-0.78125, 0.84375), x1, ii1)
        b = a * a.pol[0]
        @test b == a
        @test constant_term(a) == x1
        @test linear_polynomial(a) == Taylor1(5)

        a = TaylorModel1(x0, 5, x0, ii0)
        @test a^2 == TaylorModel1(x0^2, 5, x0, ii0)
        @test a^3 == TaylorModel1(x0^3, 5, x0, ii0)
        a = TaylorModel1(x1, 5, x1, ii1)
        @test a^2 == TaylorModel1(x1^2, 5, x1, ii1)
        @test a^3 == TaylorModel1(x1^3, 5, x1, ii1)
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), TaylorModel1(4, x0, ii0)) ==
            TaylorModel1(interval(5.0), 4, x0, ii0)
        @test rpa(x->5+one(x), TaylorModel1(4, x1, ii1)) ==
            TaylorModel1(5+x1, 4, x1, ii1)
        @test rpa(x->5*x, TaylorModel1(4, x1, ii1)) ==
            TaylorModel1(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^4, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, TaylorModel1(3, x0, ii0)) ==
            TaylorModel1( Taylor1(x0, 3), 5*ii0^4, x0, ii0)

        # Testing remainders of an RPA
        order = 2
        ii = ii0
        xx = x0
        ftest = x -> exp(x)-1
        tma = rpa(ftest, TaylorModel1(order, xx, ii))
        tmb = ftest(TaylorModel1(order, xx, ii))
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        # test for TM with scalar coefficients
        @test fp_rpa(tmc) == tmc

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> exp(x)
        tma = rpa(ftest, TaylorModel1(order, xx, ii))
        tmb = ftest(TaylorModel1(order, xx, ii))
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 3
        ii = ii0
        xx = x0
        ftest = x -> sin(x)
        tma = rpa(ftest, TaylorModel1(order, xx, ii))
        tmb = ftest(TaylorModel1(order, xx, ii))
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> sqrt(x)
        tma = rpa(ftest, TaylorModel1(order, xx, ii))
        tmb = ftest(TaylorModel1(order, xx, ii))
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 5
        ii = ii1
        xx = x1
        ftest = x -> inv(x)
        tma = rpa(ftest, TaylorModel1(order, xx, ii))
        tmb = ftest(TaylorModel1(order, xx, ii))
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.hi)-tmc.pol(ii.hi-ξ0),
                        ftest(ii.lo)-tmc.pol(ii.lo-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tma = rpa(ftest, TaylorModel1(order, xx, ii))
        tmb = ftest(TaylorModel1(order, xx, ii))
        @test remainder(tmb) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
    end

    @testset "Composition of functions and their inverses" begin
        tv = TaylorModel1(2, x0, ii0)

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
        tv = TaylorModel1(2, x1, ii1)

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
        tm = TaylorModel1(order, x0, ii0)

        integ_res = integrate(exp(tm), 1..1)
        exact_res = exp(tm)
        @test exact_res.pol == integ_res.pol
        @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_containment(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol == integ_res.pol
        @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_containment(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol == integ_res.pol
        @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_containment(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol == integ_res.pol
        # @test exact_res.rem ⊆ integ_res.rem
        for ind = 1:_num_tests
            @test check_containment(atan, integ_res)
        end
    end

    @testset "Display" begin
        tm = TaylorModel1(2, x1, ii1)
        use_show_default(true)
        @test string(exp(tm)) == "TaylorModel1{Interval{Float64},Float64}" *
            "(Taylor1{Interval{Float64}}(Interval{Float64}" *
            "[Interval(2.718281828459045, 2.7182818284590455), Interval(2.718281828459045, 2.7182818284590455), " *
            "Interval(1.3591409142295225, 1.3591409142295228)], 2), Interval(-0.05020487208677582, 0.06448109909211741), " *
            "Interval(1.0, 1.0), Interval(0.5, 1.5))"
        use_show_default(false)
        @test string(tm^3) == " Interval(1.0, 1.0) + Interval(3.0, 3.0) t + " *
            "Interval(3.0, 3.0) t² + Interval(-0.125, 0.125)"
        @test string(exp(tm)) == " Interval(2.718281828459045, 2.7182818284590455) + " *
            "Interval(2.718281828459045, 2.7182818284590455) t + " *
            "Interval(1.3591409142295225, 1.3591409142295228) t² + " *
            "Interval(-0.05020487208677582, 0.06448109909211741)"
    end
end

@testset "Tests for RTaylorModel1 " begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)
    x1 = Interval(1.0)
    ii1 = Interval(0.5, 1.5)

    @testset "RTaylorModel1 constructors" begin
        tv = RTaylorModel1{Interval{Float64},Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == RTaylorModel1(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == RTaylorModel1(5, x0, ii0)
        @test tv == RTaylorModel1(5, ii0)
        @test tv == RTaylorModel1(5, 0.0, ii0)
        @test RTaylorModel1(x1, 5, x0, ii0) == RTaylorModel1(Taylor1(x1, 5), x0, x0, ii0)
        @test RTaylorModel1(5, 0.7, ii1) == RTaylorModel1(5, interval(0.7), ii1)

        @test isa(tv, AbstractSeries)
        @test RTaylorModel1{Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Zero may not be contained in the remainder of a RTaylorModel1
        @test 0 ∉ remainder(RTaylorModel1(Taylor1(Interval{Float64},5), x1, x0, ii0))

        # Test errors in construction
        @test_throws AssertionError RTaylorModel1(5, x1, ii0)

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
        @test constant_term(tv) == interval(0.0)
        @test linear_polynomial(tv) == Taylor1(Interval{Float64},5)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = RTaylorModel1(x1+Taylor1(5), Δ, x1, ii1)
        tv = RTaylorModel1(5, x1, ii1)

        @test zero(a) == RTaylorModel1(zero(a.pol), 0..0, x1, ii1)
        @test one(a) == RTaylorModel1(one(a.pol), 0..0, x1, ii1)
        @test a+x1 == RTaylorModel1(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == RTaylorModel1(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == RTaylorModel1(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == RTaylorModel1(zero(a.pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == RTaylorModel1(a.pol*tv.pol, a.rem*tv.pol(ii1-x1), x1, ii1)
        @test b/tv == RTaylorModel1(a.pol, Interval(-2.75, 4.75), x1, ii1)
        b = a * a.pol[0]
        @test b == a
        @test constant_term(a) == x1
        @test linear_polynomial(a) == Taylor1(5)

        a = RTaylorModel1(x0, 5, x0, ii0)
        @test a^2 == RTaylorModel1(x0^2, 5, x0, ii0)
        @test a^3 == RTaylorModel1(x0^3, 5, x0, ii0)
        a = RTaylorModel1(x1, 5, x1, ii1)
        @test a^2 == RTaylorModel1(x1^2, 5, x1, ii1)
        @test a^3 == RTaylorModel1(x1^3, 5, x1, ii1)
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1(interval(5.0), 4, x0, ii0)
        @test rpa(x->5*x, RTaylorModel1(4, x1, ii1)) ==
            RTaylorModel1(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^4, RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, RTaylorModel1(3, x0, ii0)) ==
            RTaylorModel1( Taylor1(x0, 3), interval(5), x0, ii0)

        # Testing remainders and inclusion of RPAs
        order = 2
        ii = ii0
        xx = x0
        ftest = x -> exp(x)-1
        tma = rpa(ftest, RTaylorModel1(order, xx, ii))
        tmb = ftest(RTaylorModel1(order, xx, ii))
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> exp(x)
        tma = rpa(ftest, RTaylorModel1(order, xx, ii))
        tmb = ftest(RTaylorModel1(order, xx, ii))
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 3
        ii = ii0
        xx = x0
        ftest = x -> sin(x)
        tma = rpa(ftest, RTaylorModel1(order, xx, ii))
        tmb = ftest(RTaylorModel1(order, xx, ii))
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> sqrt(x)
        tma = rpa(ftest, RTaylorModel1(order, xx, ii))
        tmb = ftest(RTaylorModel1(order, xx, ii))
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        order = 5
        ii = ii1
        xx = x1
        ftest = x -> inv(x)
        tma = rpa(ftest, RTaylorModel1(order, xx, ii))
        tmb = ftest(RTaylorModel1(order, xx, ii))
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.hi)-tmc.pol(ii.hi-ξ0),
                        ftest(ii.lo)-tmc.pol(ii.lo-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tma = rpa(ftest, RTaylorModel1(order, xx, ii))
        tmb = ftest(RTaylorModel1(order, xx, ii))
        @test remainder(tmb) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
    end

    @testset "Composition of functions and their inverses" begin
        tv = RTaylorModel1(2, x0, ii0)

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
        tv = RTaylorModel1(2, x1, ii1)

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
        tm = RTaylorModel1(order, x0, ii0)

        integ_res = integrate(exp(tm), 1..1)
        exact_res = exp(tm)
        @test exact_res.pol == integ_res.pol
        @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol == integ_res.pol
        @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol == integ_res.pol
        @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol == integ_res.pol
        # @test exact_res.rem*(ii0-x0)^(order+1) ⊆ integ_res.rem*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(atan, integ_res)
        end
    end

    @testset "Display" begin
        tm = RTaylorModel1(3, x1, ii1)
        use_show_default(true)
        @test string(exp(tm)) == "RTaylorModel1{Interval{Float64},Float64}" *
            "(Taylor1{Interval{Float64}}(Interval{Float64}" *
            "[Interval(2.718281828459045, 2.7182818284590455), Interval(2.718281828459045, 2.7182818284590455), " *
            "Interval(1.3591409142295225, 1.3591409142295228), Interval(0.45304697140984085, 0.45304697140984096)], 3), " *
            "Interval(0.10281598943126724, 0.1256036426541982), Interval(1.0, 1.0), Interval(0.5, 1.5))"
        use_show_default(false)
        @test string(tm^3) == " Interval(1.0, 1.0) + Interval(3.0, 3.0) t + " *
            "Interval(3.0, 3.0) t² + Interval(1.0, 1.0) t³ + Interval(0.0, 0.0) t⁴"
        @test string(exp(tm)) == " Interval(2.718281828459045, 2.7182818284590455) + " *
            "Interval(2.718281828459045, 2.7182818284590455) t + Interval(1.3591409142295225, 1.3591409142295228) t² + " *
            "Interval(0.45304697140984085, 0.45304697140984096) t³ + Interval(0.10281598943126724, 0.1256036426541982) t⁴"
    end
end
