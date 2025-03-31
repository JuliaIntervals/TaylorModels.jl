# Tests using TaylorModel1

using TaylorModels

using Test, Random


const _num_tests = 1000
const α_mid = TaylorModels.α_mid

setformat(:full)

function check_containment(ftest, tma::RTaylorModel1)
    x0 = expansion_point(tma)
    xfp = diam(domain(tma))*(rand()-0.5) + mid(x0)
    xbf = big(xfp)
    range = tma((xfp .. xfp)-x0)
    bb = ftest(xbf) ∈ range
    bb || @show(ftest, xfp, xbf, ftest(xbf), range)
    return bb
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

        @test RTaylorModel1(tv, ii0) == RTaylorModel1(Taylor1(Interval{Float64}, 5), ii0, x0, ii0)
        @test RTaylorModel1(5, x0, ii0) == RTaylorModel1(tv, x0)
        @test RTaylorModel1(5, ii0) == RTaylorModel1(tv, x0)

        @test isa(tv, AbstractSeries)
        @test RTaylorModel1{Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Zero may not be contained in the remainder of a RTaylorModel1
        @test 0 ∉ remainder(RTaylorModel1(Taylor1(Interval{Float64},5), x1, x0, ii0))

        # Test errors in construction
        @test_throws AssertionError RTaylorModel1(5, x1, ii0)

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
        @test polynomial(tv) == Taylor1(Interval{Float64},5)
        @test domain(tv) == ii0
        @test expansion_point(tv) == x0
        @test constant_term(tv) == interval(0.0)
        @test linear_polynomial(tv) == Taylor1(Interval{Float64},5)
        @test nonlinear_polynomial(tv) == zero(Taylor1(Interval{Float64},5))
        @test centered_dom(tv) == ii0
        @test centered_dom(RTaylorModel1(5, 0.7, ii1)) == ii1-0.7

        # Tests related to fixorder
        a = RTaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        b = RTaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)
        aa, bb = TaylorModels.fixorder(a, b)
        @test get_order(aa) == get_order(bb) == 1
        @test isa(aa, RTaylorModel1) == isa(bb, RTaylorModel1) == true
        @test aa == a
        @test bb == RTaylorModel1(Taylor1([1.0, 1]), -1 .. 2, 0..0, -1 .. 1)
        # a and b remain the same
        @test a == RTaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        @test b == RTaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = RTaylorModel1(x1+Taylor1(5), Δ, x1, ii1)
        tv = RTaylorModel1(5, x1, ii1)
        a_pol = polynomial(a)
        tv_pol = polynomial(tv)

        @test zero(a) == RTaylorModel1(zero(a_pol), 0..0, x1, ii1)
        @test one(a) == RTaylorModel1(one(a_pol), 0..0, x1, ii1)
        @test a+x1 == RTaylorModel1(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == RTaylorModel1(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == RTaylorModel1(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == RTaylorModel1(zero(a_pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == RTaylorModel1(a_pol*tv_pol, remainder(a)*tv_pol(ii1-x1), x1, ii1)
        @test remainder(b/tv) ⊆ Interval(-2.75, 4.75)
        @test constant_term(b) == 1..1
        @test linear_polynomial(b) == 2*x1*Taylor1(5)
        @test nonlinear_polynomial(b) == x1*Taylor1(5)^2
        b = a * a_pol[0]
        @test b == a
        @test constant_term(a) == x1
        @test linear_polynomial(a) == Taylor1(5)
        @test nonlinear_polynomial(a) == Taylor1(0..0, 5)

        a = RTaylorModel1(x0, 5, x0, ii0)
        @test a^0 == RTaylorModel1(x0^0, 5, x0, ii0)
        @test a^1 == RTaylorModel1(x0^1, 5, x0, ii0)
        @test a^2 == RTaylorModel1(x0^2, 5, x0, ii0)
        @test a^3 == RTaylorModel1(x0^3, 5, x0, ii0)
        a = RTaylorModel1(x1, 5, x1, ii1)
        @test a^0 == RTaylorModel1(x1^0, 5, x1, ii1)
        @test a^1 == RTaylorModel1(x1^1, 5, x1, ii1)
        @test a^2 == RTaylorModel1(x1^2, 5, x1, ii1)
        @test a^3 == RTaylorModel1(x1^3, 5, x1, ii1)

        # Tests involving RTM1s with different orders
        a = RTaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        b = RTaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)
        aa, bb = TaylorModels.fixorder(a, b)
        @test get_order(aa) == get_order(bb)
        @test get_order(bb) == 1
        @test a + b == aa + bb
        @test a - b == aa - bb
        res1 = a * b
        res2 = aa * bb
        @test res1 == RTaylorModel1(Taylor1([1.0, 2]), -1 .. 9 , 0..0, -1 .. 1)
        @test res2 == RTaylorModel1(Taylor1([1.0, 2]), -2 .. 9 , 0..0, -1 .. 1)
        res1 = a / b
        res2 = aa / bb
        @test res1 == RTaylorModel1(Taylor1([1.0, 0]), entireinterval() , 0..0, -1 .. 1)
        @test res2 == res1
        # a and b remain the same
        @test a == RTaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        @test b == RTaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)

        @test_throws AssertionError a+RTaylorModel1(a.pol, remainder(a), 1..1, -1..1)
        @test_throws AssertionError a+RTaylorModel1(a.pol, remainder(a), 0..0, -2..2)
        f(x) = x + x^2
        tm = RTaylorModel1(5, 0.0, -0.5 .. 0.5)
        @test f(tm)/tm == 1+tm
    end

    @testset "RTM1's with TaylorN coefficients" begin
        # Tests for RTM1's with TaylorN coefficients
        orderT = 4
        orderQ = 5
        ξ = set_variables("ξ", order = 2 * orderQ, numvars=1)
        q0 = [0.5]
        δq0 = IntervalBox(-0.1 .. 0.1, Val(1))
        qaux = normalize_taylor(q0[1] + TaylorN(1, order=orderQ), δq0, true)
        symIbox = IntervalBox(-1 .. 1, Val(1))
        t = Taylor1([qaux, one(qaux)], orderT)
        dom = -0.5 .. 0.5
        x00 = mid(dom)

        f(x) = x + x^2
        g(x) = x
        h(x) = x^3*(1+x)

        tm = RTaylorModel1(deepcopy(t), 0 .. 0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test isentire(remainder(fgTM1))

        fgTM1 = f(tm) * (g(tm))^2
        hh = h(tm)
        @test polynomial(fgTM1) ≈ polynomial(hh)
        @test remainder(fgTM1) == remainder(hh)
        for ind = 1:_num_tests
            xξ = rand(dom)-x00
            qξ = rand(symIbox)
            tt = tm(xξ)(qξ)
            @test h(tt) ⊆ fgTM1(dom-x00)(symIbox)
        end

        t = Taylor1([one(qaux), qaux], orderT)
        tm = RTaylorModel1(deepcopy(t), 0 .. 0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test !isentire(remainder(fgTM1))
        for ind = 1:_num_tests
            xξ = rand(dom)-x00
            qξ = rand(symIbox)
            tt = 1+t(xξ)(qξ)
            @test tt ⊆ fgTM1(dom-x00)(symIbox)
        end

        # Testing integration
        @test integrate(tm, symIbox) == RTaylorModel1(integrate(t), 0..0, x00, dom)
        @test integrate(f(tm), symIbox) == RTaylorModel1(integrate(f(t)), 0..0, x00, dom)
        t = Taylor1([qaux, one(qaux)], orderT)
        tm = RTaylorModel1(deepcopy(t), -0.25 .. 0.25, x00, dom)
        @test integrate(tm, symIbox) == RTaylorModel1(integrate(t),
            remainder(tm)*(domain(tm)-expansion_point(tm))/(orderT+2), x00, dom)

        # Missing tests: Changing order of a RTM1 with TaylorN coeffs (see TM1.jl)
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1(interval(5.0), 4, x0, ii0)
        @test rpa(x->5+one(x), RTaylorModel1(4, x1, ii1)) ==
            RTaylorModel1(5+x1, 4, x1, ii1)
        @test rpa(x->5*x, RTaylorModel1(4, x1, ii1)) ==
            RTaylorModel1(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^0, RTaylorModel1(4, x0, ii0)) == 5*RTaylorModel1(4, x0, ii0)^0
        @test rpa(x->5*x^0, RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1( interval(5.0)*Taylor1(4)^0, x0, x0, ii0)
        @test rpa(x->5*x^1, RTaylorModel1(4, x0, ii0)) == 5*RTaylorModel1(4, x0, ii0)^1
        @test rpa(x->5*x^1, RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1( interval(5.0)*Taylor1(4)^1, x0, x0, ii0)
        @test rpa(x->5*x^2, RTaylorModel1(4, x0, ii0)) == 5*RTaylorModel1(4, x0, ii0)^2
        @test rpa(x->5*x^2, RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1( interval(5.0)*Taylor1(4)^2, x0, x0, ii0)
        @test rpa(x->5*x^4, RTaylorModel1(4, x0, ii0)) ==
            RTaylorModel1( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, RTaylorModel1(3, x0, ii0)) ==
            RTaylorModel1( Taylor1(x0, 3), interval(5), x0, ii0)

        # Testing remainders and inclusion of RPAs
        order = 2
        ii = ii0
        xx = x0
        ftest = x -> exp(x)-1
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        # test for TM with scalar coefficients
        @test fp_rpa(tmc) == tmc

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> exp(x)
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        order = 3
        ii = ii0
        xx = x0
        ftest = x -> sin(x)
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> sqrt(x)
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        order = 5
        ii = ii1
        xx = x1
        ftest = x -> inv(x)
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        # fT, Δ, ξ0, δ = fp_rpa(tma)
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test interval(ftest(ii.hi)-tmc.pol(ii.hi-ξ0),
                        ftest(ii.lo)-tmc.pol(ii.lo-ξ0)) ⊆ remainder(tma)*(ii-ξ0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test remainder(tmb) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))
    end

    @testset "RPAs with polynomial Taylor1{TaylorN{T}}" begin
        orderT = 5
        orderQ = 5
        dom = 0 .. 1
        t00 = mid(dom)
        symIbox = IntervalBox(-1 .. 1, 1)
        δq0 = IntervalBox(-0.25 .. 0.25, 1)
        qaux = normalize_taylor(TaylorN(1, order=orderQ) + t00, δq0, true)
        xT = Taylor1([qaux, one(qaux)], orderT)
        tm = RTaylorModel1(deepcopy(xT), 0 .. 0, t00, dom)

        s(x) = sin(x)
        c(x) = cos(x)
        ee(x) = exp(x)
        p5(x) = x^5
        lone(x) = log(1.0+x)
        pol(x) = x^3 / x^5

        sT = s(tm)
        cT = c(tm)
        eeT = ee(tm)
        p5T = p5(tm)
        loneT = lone(tm)
        polT = pol(tm)

        for ind = 1:_num_tests
            xξ = rand(dom)
            q0ξ = t00 + rand(δq0)[1]
            t = Taylor1(2*orderT) + q0ξ
            st = s(t)
            ct = c(t)
            eet = ee(t)
            p5t = p5(t)
            lonet = lone(t)
            polt = pol(t)

            @test s(xξ) ⊆ sT(xξ - t00)(symIbox)
            @test st(xξ - q0ξ) ⊆ sT(xξ - t00)(symIbox)
            @test c(xξ) ⊆ cT(xξ - t00)(symIbox)
            @test ct(xξ - q0ξ) ⊆ cT(xξ - t00)(symIbox)
            @test ee(xξ) ⊆ eeT(xξ - t00)(symIbox)
            @test eet(xξ - q0ξ) ⊆ eeT(xξ - t00)(symIbox)
            @test p5(xξ) ⊆ p5T(xξ - t00)(symIbox)
            @test p5t(xξ - q0ξ) ⊆ p5T(xξ - t00)(symIbox)
            @test lone(xξ) ⊆ loneT(xξ - t00)(symIbox)
            @test lonet(xξ - q0ξ) ⊆ loneT(xξ - t00)(symIbox)
            @test p5(xξ) ⊆ p5T(xξ - t00)(symIbox)
            @test p5t(xξ - q0ξ) ⊆ p5T(xξ - t00)(symIbox)
        end
    end

    @testset "Composition of functions and their inverses" begin
        tv = RTaylorModel1(2, x0, ii0)
        tv_pol = polynomial(tv)

        tma = exp(tv)
        tmb = log(tma)
        @test tmb == log(exp(tv))
        @test tmb.pol == tv_pol

        tma = sin(tv)
        tmb = asin(tma)
        @test tmb == asin(sin(tv))
        @test tmb.pol == tv_pol

        tma = asin(tv)
        tmb = sin(tma)
        @test tmb == sin(asin(tv))
        @test tmb.pol == tv_pol

        tma = acos(tv)
        tmb = cos(tma)
        @test tmb == cos(acos(tv))
        @test sup(norm(tmb.pol - tv_pol, Inf)) < 5.0e-16

        tma = tan(tv)
        tmb = atan(tma)
        @test tmb == atan(tan(tv))
        @test tmb.pol == tv_pol

        tma = atan(tv)
        tmb = tan(tma)
        @test tmb == tan(atan(tv))
        @test tmb.pol == tv_pol


        ####
        tv = RTaylorModel1(2, x1, ii1)
        tv_pol = polynomial(tv)

        tma = log(tv)
        tmb = exp(tma)
        @test tmb == exp(log(tv))
        @test tmb.pol == tv_pol

        tma = cos(tv)
        tmb = acos(tma)
        @test tmb == acos(cos(tv))
        @test sup(norm(tmb.pol - tv_pol, Inf)) < 1.0e-15
    end

    @testset "Tests for integrate" begin
        order = 4
        tm = RTaylorModel1(order, x0, ii0)

        integ_res = integrate(exp(tm), 1..1)
        exact_res = exp(tm)
        @test exact_res.pol == integ_res.pol
        @test remainder(exact_res)*(ii0-x0)^(order+1) ⊆ remainder(integ_res)*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol == integ_res.pol
        @test remainder(exact_res)*(ii0-x0)^(order+1) ⊆ remainder(integ_res)*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol == integ_res.pol
        @test remainder(exact_res)*(ii0-x0)^(order+1) ⊆ remainder(integ_res)*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol == integ_res.pol
        # @test remainder(exact_res)*(ii0-x0)^(order+1) ⊆ remainder(integ_res)*(ii0-x0)^(order+1)
        for ind = 1:_num_tests
            @test check_containment(atan, integ_res)
        end
    end

    @testset "Display" begin
        tm = RTaylorModel1(3, x1, ii1)
        use_show_default(true)
        if VERSION < v"1.6"
            @test string(exp(tm)) == "RTaylorModel1{Interval{Float64},Float64}" *
                "(Taylor1{Interval{Float64}}(Interval{Float64}" *
                "[Interval(2.718281828459045, 2.7182818284590455), Interval(2.718281828459045, 2.7182818284590455), " *
                "Interval(1.3591409142295225, 1.3591409142295228), Interval(0.45304697140984085, 0.45304697140984096)], 3), " *
                "Interval(0.10281598943126369, 0.1256036426541982), Interval(1.0, 1.0), Interval(0.5, 1.5))"
        else
            @test string(exp(tm)) == "RTaylorModel1{Interval{Float64}, Float64}" *
                "(Taylor1{Interval{Float64}}(Interval{Float64}" *
                "[Interval(2.718281828459045, 2.7182818284590455), Interval(2.718281828459045, 2.7182818284590455), " *
                "Interval(1.3591409142295225, 1.3591409142295228), Interval(0.45304697140984085, 0.45304697140984096)], 3), " *
                "Interval(0.10281598943126369, 0.1256036426541982), Interval(1.0, 1.0), Interval(0.5, 1.5))"
        end
        use_show_default(false)
        @test string(tm^3) == " Interval(1.0, 1.0) + Interval(3.0, 3.0) t + " *
            "Interval(3.0, 3.0) t² + Interval(1.0, 1.0) t³ + Interval(0.0, 0.0) t⁴"
        @test string(exp(tm)) == " Interval(2.718281828459045, 2.7182818284590455) + " *
            "Interval(2.718281828459045, 2.7182818284590455) t + Interval(1.3591409142295225, 1.3591409142295228) t² + " *
            "Interval(0.45304697140984085, 0.45304697140984096) t³ + Interval(0.10281598943126369, 0.1256036426541982) t⁴"
    end
end
