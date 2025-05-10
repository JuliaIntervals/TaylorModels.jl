# Tests using TaylorModel1

using TaylorModels

using Test, Random

const _num_tests = 1000
const α_mid = TM.α_mid

setdisplay(:full)

function check_containment(ftest, tma::RTaylorModel1)
    x0 = expansion_point(tma)
    xfp = sample(domain(tma))#diam(domain(tma))*(rand()-0.5) + mid(x0)
    xbf = big(xfp)
    range = tma(interval(xfp, xfp)-x0)
    bb = in_interval(ftest(xbf), range)
    bb || @show(ftest, xfp, xbf, ftest(xbf), range)
    return bb
end

@testset "Tests for RTaylorModel1 " begin
    x0 = interval(0.0)
    ii0 = interval(-0.5, 0.5)
    x1 = interval(1.0)
    ii1 = interval(0.5, 1.5)
    y0 = interval(0, 1)
    y1 = interval(-1, 1)

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
        @test isequal_interval(remainder(tv), interval(0.0))
        @test polynomial(tv) == Taylor1(Interval{Float64},5)
        @test isequal_interval(domain(tv), ii0)
        @test isequal_interval(expansion_point(tv), x0)
        @test isequal_interval(constant_term(tv), interval(0.0))
        @test linear_polynomial(tv) == Taylor1(Interval{Float64},5)
        @test nonlinear_polynomial(tv) == zero(Taylor1(Interval{Float64},5))
        @test isequal_interval(centered_dom(tv), ii0)
        @test isequal_interval(centered_dom(RTaylorModel1(5, 0.7, ii1)), ii1-0.7)

        # Tests related to fixorder
        a = RTaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        b = RTaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)
        aa, bb = TM.fixorder(a, b)
        @test get_order(aa) == get_order(bb) == 1
        @test isa(aa, RTaylorModel1) == isa(bb, RTaylorModel1) == true
        @test aa == a
        @test bb == RTaylorModel1(Taylor1([1.0, 1]), interval(-1, 2), x0, y1)
        # a and b remain the same
        @test a == RTaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        @test b == RTaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = RTaylorModel1(x1+Taylor1(5), Δ, x1, ii1)
        tv = RTaylorModel1(5, x1, ii1)
        a_pol = polynomial(a)
        tv_pol = polynomial(tv)

        @test zero(a) == RTaylorModel1(zero(a_pol), x0, x1, ii1)
        @test one(a) == RTaylorModel1(one(a_pol), x0, x1, ii1)
        @test a+x1 == RTaylorModel1(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == RTaylorModel1(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == RTaylorModel1(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == RTaylorModel1(zero(a_pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == RTaylorModel1(a_pol*tv_pol, remainder(a)*tv_pol(ii1-x1), x1, ii1)
        @test issubset_interval(remainder(b/tv), interval(-2.75, 4.75))
        @test isequal_interval(constant_term(b), interval(1))
        @test linear_polynomial(b) == 2*x1*Taylor1(5)
        @test nonlinear_polynomial(b) == x1*Taylor1(5)^2
        b = a * a_pol[0]
        @test b == a
        @test isequal_interval(constant_term(a), x1)
        @test linear_polynomial(a) == Taylor1(5)
        @test nonlinear_polynomial(a) == Taylor1(x0, 5)

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
        a = RTaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        b = RTaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)
        aa, bb = TM.fixorder(a, b)
        @test get_order(aa) == get_order(bb)
        @test get_order(bb) == 1
        @test a + b == aa + bb
        @test a - b == aa - bb
        res1 = a * b
        res2 = aa * bb
        @test res1 == RTaylorModel1(Taylor1([1.0, 2]), interval(-1, 9) , x0, y1)
        @test res2 == RTaylorModel1(Taylor1([1.0, 2]), interval(-2, 9) , x0, y1)
        res1 = a / b
        res2 = aa / bb
        @test res1 == RTaylorModel1(Taylor1([1.0, 0]), entireinterval() , x0, y1)
        @test res2 == res1
        # a and b remain the same
        @test a == RTaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        @test b == RTaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)

        @test_throws AssertionError a+RTaylorModel1(a.pol, remainder(a), interval(1), y1)
        @test_throws AssertionError a+RTaylorModel1(a.pol, remainder(a), x0, interval(-2,2))
        f(x) = x + x^2
        tm = RTaylorModel1(5, 0.0, interval(-0.5, 0.5))
        @test f(tm)/tm == 1+tm
    end

    @testset "RTM1's with TaylorN coefficients" begin
        # Tests for RTM1's with TaylorN coefficients
        orderT = 4
        orderQ = 5
        ξ = set_variables("ξ", order = 2 * orderQ, numvars=1)
        q0 = [0.5]
        δq0 = [interval(-0.1, 0.1)]
        qaux = normalize_taylor(q0[1] + TaylorN(1, order=orderQ), δq0, true)
        # qaux = normalize_taylor(q0[1] + TaylorN(Interval{Float64}, 1, order=orderQ), δq0, true)
        symIbox = symmetric_box(Float64, 1)
        t = Taylor1([qaux, one(qaux)], orderT)
        dom = interval(-0.5, 0.5)
        x00 = mid(dom)

        f(x) = x + x^2
        g(x) = x
        h(x) = x^3*(1+x)

        tm = RTaylorModel1(deepcopy(t), x0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test isentire_interval(remainder(fgTM1))

        fgTM1 = f(tm) * (g(tm))^2
        hh = h(tm)
        @test_skip polynomial(fgTM1) ≈ polynomial(hh)
        @test remainder(fgTM1) == remainder(hh)
        for ind = 1:_num_tests
            xξ = sample(dom)-x00
            qξ = sample.(symIbox)
            tt = tm(xξ)(qξ)
            @test issubset_interval(h(tt), fgTM1(dom-x00)(symIbox))
        end

        t = Taylor1([one(qaux), qaux], orderT)
        tm = RTaylorModel1(deepcopy(t), x0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test !isentire_interval(remainder(fgTM1))
        for ind = 1:_num_tests
            xξ = sample(dom)-x00
            qξ = sample.(symIbox)
            tt = 1+t(xξ)(qξ)
            @test issubset_interval(tt, fgTM1(dom-x00)(symIbox))
        end

        # Testing integration
        @test integrate(tm, symIbox) == RTaylorModel1(integrate(t), x0, x00, dom)
        @test integrate(f(tm), symIbox) == RTaylorModel1(integrate(f(t)), x0, x00, dom)
        t = Taylor1([qaux, one(qaux)], orderT)
        tm = RTaylorModel1(deepcopy(t), interval(-0.25, 0.25), x00, dom)
        @test integrate(tm, symIbox) == RTaylorModel1(integrate(t),
            remainder(tm)*(domain(tm)-expansion_point(tm))/(orderT+2), x00, dom)
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
        @test issubset_interval( interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
            ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma)*(ii-ξ0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+interval(1))

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
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                        ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma)*(ii-ξ0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
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
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                        ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma)*(ii-ξ0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
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
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                        ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma)*(ii-ξ0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
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
        @test issubset_interval(interval(ftest(sup(ii))-tmc.pol(sup(ii)-ξ0),
                        ftest(inf(ii))-tmc.pol(inf(ii)-ξ0)), remainder(tma)*(ii-ξ0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tm = RTaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test issubset_interval(remainder(tmb), remainder(tma))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))
    end

    @testset "RPAs with polynomial Taylor1{TaylorN{T}}" begin
        orderT = 5
        orderQ = 5
        dom = y0
        t00 = mid(dom)
        symIbox = symmetric_box(Float64, 1)
        δq0 = [interval(-0.25, 0.25)]
        qaux = normalize_taylor(TaylorN(1, order=orderQ) + t00, δq0, true)
        # qaux = normalize_taylor(TaylorN(Interval{Float64}, 1, order=orderQ) + t00, δq0, true)
        xT = Taylor1([qaux, one(qaux)], orderT)
        tm = RTaylorModel1(deepcopy(xT), x0, t00, dom)

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
            xξ = sample(dom)
            q0ξ = t00 + sample.(δq0)[1]
            t = Taylor1(2*orderT) + q0ξ
            st = s(t)
            ct = c(t)
            eet = ee(t)
            p5t = p5(t)
            lonet = lone(t)
            polt = pol(t)

            @test in_interval(s(xξ), sT(xξ - t00)(symIbox))
            @test in_interval(st(xξ - q0ξ), sT(xξ - t00)(symIbox))
            @test in_interval(c(xξ), cT(xξ - t00)(symIbox))
            @test in_interval(ct(xξ - q0ξ), cT(xξ - t00)(symIbox))
            @test in_interval(ee(xξ), eeT(xξ - t00)(symIbox))
            @test in_interval(eet(xξ - q0ξ), eeT(xξ - t00)(symIbox))
            @test in_interval(p5(xξ), p5T(xξ - t00)(symIbox))
            @test in_interval(p5t(xξ - q0ξ), p5T(xξ - t00)(symIbox))
            @test in_interval(lone(xξ), loneT(xξ - t00)(symIbox))
            @test in_interval(lonet(xξ - q0ξ), loneT(xξ - t00)(symIbox))
            @test in_interval(p5(xξ), p5T(xξ - t00)(symIbox))
            @test in_interval(p5t(xξ - q0ξ), p5T(xξ - t00)(symIbox))
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

        integ_res = integrate(exp(tm), x1)
        exact_res = exp(tm)
        @test exact_res.pol == integ_res.pol
        @test issubset_interval(remainder(exact_res)*(ii0-x0)^(order+1),
                remainder(integ_res)*(ii0-x0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol == integ_res.pol
        @test issubset_interval(remainder(exact_res)*(ii0-x0)^(order+1),
                remainder(integ_res)*(ii0-x0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), x1)
        exact_res = cos(tm)
        @test exact_res.pol == integ_res.pol
        @test issubset_interval(remainder(exact_res)*(ii0-x0)^(order+1),
                remainder(integ_res)*(ii0-x0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol == integ_res.pol
        @test_broken issubset_interval(remainder(exact_res)*(ii0-x0)^(order+1), remainder(integ_res)*(ii0-x0)^(order+1))
        for ind = 1:_num_tests
            @test check_containment(atan, integ_res)
        end
    end

    @testset "Display" begin
        tm = RTaylorModel1(3, x1, ii1)
        use_show_default(true)
        @test string(exp(tm)) == "RTaylorModel1{Interval{Float64}, Float64}" *
            "(Taylor1{Interval{Float64}}(Interval{Float64}[" *
            "Interval{Float64}(2.718281828459045, 2.7182818284590455, com), " *
            "Interval{Float64}(2.718281828459045, 2.7182818284590455, com), " *
            "Interval{Float64}(1.3591409142295225, 1.3591409142295228, com), " *
            "Interval{Float64}(0.45304697140984085, 0.45304697140984096, com)], 3), " *
            "Interval{Float64}(0.10281598943126369, 0.1256036426541982, trv)_NG, " *
            "Interval{Float64}(1.0, 1.0, com), Interval{Float64}(0.5, 1.5, com)_NG)"
        use_show_default(false)
        @test string(tm^3) == " Interval{Float64}(1.0, 1.0, com) + " *
            "Interval{Float64}(3.0, 3.0, com) t + Interval{Float64}(3.0, 3.0, com) t² + " *
            "Interval{Float64}(1.0, 1.0, com) t³ + Interval{Float64}(0.0, 0.0, com)_NG t⁴"
        @test string(exp(tm)) == " Interval{Float64}(2.718281828459045, 2.7182818284590455, com) + " *
            "Interval{Float64}(2.718281828459045, 2.7182818284590455, com) t + " *
            "Interval{Float64}(1.3591409142295225, 1.3591409142295228, com) t² + " *
            "Interval{Float64}(0.45304697140984085, 0.45304697140984096, com) t³ + " *
            "Interval{Float64}(0.10281598943126369, 0.1256036426541982, trv)_NG t⁴"
    end
end
