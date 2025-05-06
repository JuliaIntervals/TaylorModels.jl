# Tests using TaylorModel1

using TaylorModels

using Test, Random

myrand(a::Interval) = diam(a)*rand() + inf(a)

const _num_tests = 1000
const TM = TaylorModels
const α_mid = TM.α_mid

setdisplay(:full)

function check_containment(ftest, tma::TaylorModel1)
    x0 = expansion_point(tma)
    xfp = myrand(domain(tma))#diam(domain(tma))*(rand()-0.5) + mid(x0)
    xbf = big(xfp)
    range = tma(interval(xfp, xfp)-x0)
    bb = in_interval(ftest(xbf), range)
    bb || @show(ftest, xfp, xbf, ftest(xbf), range)
    return bb
end

@testset "Test `bound_taylor1`" begin
    x0 = interval(0.0)
    x1 = interval(0,1)
    ii0 = interval(-0.5, 0.5)

    tpol = exp( Taylor1(2) )
    @test isequal_interval(TaylorModels.bound_taylor1( tpol, ii0),
            interval(tpol(inf(ii0)), tpol(sup(ii0))))
    @test isequal_interval(TaylorModels.bound_taylor1(
            exp( Taylor1(Interval{Float64}, 2) ), ii0),
            interval(tpol(inf(ii0)), tpol(sup(ii0))))

    # An uncomfortable example from Makino
    t = Taylor1(5)
    f(x) = 1 - x^4 + x^5
    @test_broken issubset_interval(interval(1-4^4/5^5,1), TM.bound_taylor1(f(t), x1))
    tm = TaylorModel1(5, x0, ii0)
    @test_broken issubset_interval(interval(1-4^4/5^5,1), TM.bound_taylor1(f(tm)))
    @test_broken issubset_interval(interval(1-4^4/5^5,1), TM.bound_taylor1(f(tm), x1))
end

@testset "Tests for TaylorModel1 " begin
    x0 = interval(0.0)
    ii0 = interval(-0.5, 0.5)
    x1 = interval(1.0)
    ii1 = interval(0.5, 1.5)
    y0 = interval(0, 1)
    y1 = interval(-1, 1)

    @testset "TaylorModel1 constructors" begin
        tv = TaylorModel1{Interval{Float64},Float64}(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TaylorModel1(Taylor1(Interval{Float64},5), x0, x0, ii0)
        @test tv == TaylorModel1(5, x0, ii0)
        @test tv == TaylorModel1(5, ii0)
        @test tv == TaylorModel1(5, 0.0, ii0)
        @test TaylorModel1(x1, 5, x0, ii0) == TaylorModel1(Taylor1(x1, 5), x0, x0, ii0)
        @test TaylorModel1(5, 0.7, ii1) == TaylorModel1(5, interval(0.7), ii1)

        @test TaylorModel1(tv, ii0) == TaylorModel1(Taylor1(Interval{Float64}, 5), ii0, x0, ii0)
        @test TaylorModel1(5, x0, ii0) == TaylorModel1(tv, x0)
        @test TaylorModel1(5, ii0) == TaylorModel1(tv, x0)

        @test isa(tv, AbstractSeries)
        @test TaylorModel1{Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Test errors in construction
        @test_throws AssertionError TaylorModel1(Taylor1(Interval{Float64},5), x1, x0, ii0)
        @test_throws AssertionError TaylorModel1(5, x1, ii0)

        # Tests for get_order, remainder, polynomial and domain
        @test get_order(tv) == 5
        @test isequal_interval(remainder(tv), interval(0.0))
        @test polynomial(tv) == Taylor1(Interval{Float64},5)
        @test isequal_interval(domain(tv), ii0)
        @test isequal_interval(expansion_point(tv), x0)
        @test isequal_interval(constant_term(tv), interval(0.0))
        @test linear_polynomial(tv) == Taylor1(Interval{Float64},5)
        @test nonlinear_polynomial(tv) == zero(Taylor1(Interval{Float64},5))
        @test isequal_interval(centered_dom(tv), ii0)
        @test isequal_interval(centered_dom(TaylorModel1(5, 0.7, ii1)), ii1-0.7)

        # Tests related to fixorder
        a = TaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        b = TaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)
        aa, bb = TaylorModels.fixorder(a, b)
        @test get_order(aa) == get_order(bb) == 1
        @test isa(aa, TaylorModel1) == isa(bb, TaylorModel1) == true
        @test aa == a
        @test bb == TaylorModel1(Taylor1([1.0, 1, 0]), interval(-1, 2), x0, y1)
        # a and b remain the same
        @test a == TaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        @test b == TaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = TaylorModel1(x1+Taylor1(5), Δ, x1, ii1)
        tv = TaylorModel1(5, x1, ii1)
        a_pol = polynomial(a)
        tv_pol = polynomial(tv)
        @test zero(a) == TaylorModel1(zero(a_pol), x0, x1, ii1)
        @test one(a) == TaylorModel1(one(a_pol), x0, x1, ii1)
        @test a+x1 == TaylorModel1(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == TaylorModel1(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == TaylorModel1(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == TaylorModel1(zero(a_pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == TaylorModel1(a_pol*tv_pol, remainder(a)*tv_pol(ii1-x1), x1, ii1)
        @test issubset_interval(remainder(b/tv), interval(-0.78125, 0.84375))
        @test isequal_interval(constant_term(b), x1)
        @test linear_polynomial(b) == 2*x1*Taylor1(5)
        @test nonlinear_polynomial(b) == x1*Taylor1(5)^2
        b = a * a_pol[0]
        @test b == a
        @test constant_term(a) == x1
        @test linear_polynomial(a) == Taylor1(5)
        @test nonlinear_polynomial(a) == Taylor1(x0, 5)

        a = TaylorModel1(x0, 5, x0, ii0)
        @test a^0 == TaylorModel1(x0^0, 5, x0, ii0)
        @test a^1 == TaylorModel1(x0^1, 5, x0, ii0)
        @test a^2 == TaylorModel1(x0^2, 5, x0, ii0)
        @test a^3 == TaylorModel1(x0^3, 5, x0, ii0)
        a = TaylorModel1(x1, 5, x1, ii1)
        @test a^0 == TaylorModel1(x1^0, 5, x1, ii1)
        @test a^1 == TaylorModel1(x1^1, 5, x1, ii1)
        @test a^2 == TaylorModel1(x1^2, 5, x1, ii1)
        @test a^3 == TaylorModel1(x1^3, 5, x1, ii1)

        # Tests involving TM1s with different orders
        a = TaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        b = TaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)
        aa, bb = TaylorModels.fixorder(a, b)
        @test get_order(aa) == get_order(bb)
        @test get_order(bb) == 1
        @test a + b == aa + bb
        @test a - b == aa - bb
        res1 = a * b
        res2 = aa * bb
        @test res1 == TaylorModel1(Taylor1([1.0, 2]), interval(-2, 9), x0, y1)
        @test res2 == TaylorModel1(Taylor1([1.0, 2]), interval(-3, 9), x0, y1)
        res1 = a / b
        res2 = aa / bb
        @test res1 == TaylorModel1(Taylor1([1.0, 0]), entireinterval() , x0, y1)
        @test res2 == res1
        # a and b remain the same
        @test a == TaylorModel1(Taylor1([1.0, 1]), y0, x0, y1)
        @test b == TaylorModel1(Taylor1([1.0, 1, 0, 1]), y0, x0, y1)

        @test_throws AssertionError a+TaylorModel1(a.pol, remainder(a), x1, y1)
        @test_throws AssertionError a+TaylorModel1(a.pol, remainder(a), x0, 2*y1)
        f(x) = x + x^2
        tm = TaylorModel1(5, x0, ii0)
        @test_throws ArgumentError f(tm)/tm
    end

    @testset "TM1's with TaylorN coefficients" begin
        # Tests for TM1's with TaylorN coefficients
        orderT = 4
        orderQ = 5
        ξ = set_variables("ξ", order = 2 * orderQ, numvars=1)
        q0 = [0.5]
        δq0 = [interval(-0.1, 0.1)]
        qaux = normalize_taylor(q0[1] +
            TaylorN(typeof(δq0[1]), 1, order=orderQ), δq0, true)
        symIbox = symmetric_box(Float64)
        t = Taylor1([qaux, one(qaux)], orderT)
        dom = y0
        x00 = mid(dom)

        f(x) = x + x^2
        g(x) = x
        h(x) = x^3*(one(x)+x)

        tm = TaylorModel1(deepcopy(t), x0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test isentire_interval(remainder(fgTM1))

        fgTM1 = f(tm) * (g(tm))^2
        hh = h(tm)
        @test_skip polynomial(fgTM1) ≈ polynomial(hh)
        @test isequal(remainder(fgTM1), remainder(hh))
        for ind = 1:_num_tests
            xξ = myrand(dom)-x00
            qξ = myrand.(symIbox)
            tt = t(xξ)(qξ)
            @test issubset_interval(h(tt), fgTM1(dom-x00)(symIbox))
        end

        t = Taylor1([one(qaux), qaux], orderT)
        tm = TaylorModel1(deepcopy(t), x0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test !isentire_interval(remainder(fgTM1))
        for ind = 1:_num_tests
            xξ = myrand(dom)-x00
            qξ = myrand.(symIbox)
            tt = 1 + t(xξ)(qξ)
            @test issubset_interval(tt, fgTM1(dom-x00)(symIbox))
        end

        # Testing integration
        @test integrate(tm, symIbox) == TaylorModel1(integrate(t), 0..0, x00, dom)
        @test integrate(f(tm), symIbox) == TaylorModel1(integrate(f(t)), 0..0, x00, dom)
        t = Taylor1([qaux,one(qaux)], orderT)
        tm = TaylorModel1(deepcopy(t), -0.25 .. 0.25, x00, dom)
        @test integrate(tm, symIbox) == TaylorModel1(integrate(t),
            remainder(tm)*(domain(tm)-expansion_point(tm)), x00, dom)

        # Changing order of a TM1 with TaylorN coeffs
        t = Taylor1([qaux, one(qaux)], orderT)
        tm = TaylorModel1(deepcopy(t), x0, x00, dom)   # order 4
        t8 = Taylor1([qaux, one(qaux)], 2*orderT)
        tm8 = TaylorModel1(deepcopy(t8), x0, x00, dom)
        tm4, _ = TaylorSeries.fixorder(tm8, tm)
        @test get_order(tm4) == get_order(tm)
        @test polynomial(tm4) == polynomial(tm)
        #
        exp_tm = exp(tm)
        exp_tm8 = exp(tm8)
        exp_tm4, _ = TaylorSeries.fixorder(exp_tm8, tm)
        @test get_order(exp_tm4) == get_order(exp_tm)
        @test polynomial(exp_tm4) == polynomial(exp_tm)
        for ind = 1:_num_tests
            xξ = myrand(dom)-x00
            qξ = myrand.(symIbox)
            tt = t(xξ)(qξ)
            @test issubset_interval(exp(tt), exp_tm4(dom-x00)(symIbox))
            @test issubset_interval(exp(tt), exp_tm(dom-x00)(symIbox))
        end
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), TaylorModel1(4, x0, ii0)) ==
            TaylorModel1(interval(5.0), 4, x0, ii0)
        @test rpa(x->5+one(x), TaylorModel1(4, x1, ii1)) ==
            TaylorModel1(5+x1, 4, x1, ii1)
        @test rpa(x->5*x, TaylorModel1(4, x1, ii1)) ==
            TaylorModel1(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^0, TaylorModel1(4, x0, ii0)) ==
            5*TaylorModel1(4, x0, ii0)^0
        @test rpa(x->5*x^0, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^0, x0, x0, ii0)
        @test rpa(x->5*x^1, TaylorModel1(4, x0, ii0)) ==
            5*TaylorModel1(4, x0, ii0)^1
        @test rpa(x->5*x^1, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^1, x0, x0, ii0)
        @test rpa(x->5*x^2, TaylorModel1(4, x0, ii0)) ==
            5*TaylorModel1(4, x0, ii0)^2
        @test rpa(x->5*x^2, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^2, x0, x0, ii0)
        @test rpa(x->5*x^4, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^4, x0, x0, ii0)
        @test rpa(x->5*x^4, TaylorModel1(3, x0, ii0)) ==
            TaylorModel1( Taylor1(x0, 3), 5*ii0^4, x0, ii0)

        # Testing remainders of an RPA
        order = 2
        ii = ii0
        xx = x0
        ftest = x -> exp(x)-1
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                            ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma))
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
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                            ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+interval(1))

        order = 3
        ii = ii0
        xx = x0
        ftest = x -> sin(x)
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                            ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+interval(1))

        order = 2
        ii = ii1
        xx = x1
        ftest = x -> sqrt(x)
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test issubset_interval(interval(ftest(inf(ii))-tmc.pol(inf(ii)-ξ0),
                            ftest(sup(ii))-tmc.pol(sup(ii)-ξ0)), remainder(tma))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+interval(1))

        order = 5
        ii = ii1
        xx = x1
        ftest = x -> inv(x)
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test tma == tmb
        ξ0 = mid(xx, α_mid)
        tmc = fp_rpa(tma)
        @test issubset_interval(interval(ftest(sup(ii))-tmc.pol(sup(ii)-ξ0),
                            ftest(inf(ii))-tmc.pol(inf(ii)-ξ0)), remainder(tma))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+interval(1))

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test issubset_interval(remainder(tmb), remainder(tma))
        for ind = 1:_num_tests
            @test check_containment(ftest, tma)
        end
        @test_throws AssertionError tmb(sup(ii)+1.0)
        @test_throws AssertionError tmb(ii+interval(1))
    end

    @testset "RPAs with polynomial Taylor1{TaylorN{T}}" begin
        orderT = 5
        orderQ = 5
        dom = y0
        t00 = mid(dom)
        symIbox = symmetric_box(Float64)
        δq0 = [interval(-0.25, 0.25)]
        qaux = normalize_taylor(TaylorN(Interval{Float64}, 1, order=orderQ) + t00, δq0, true)
        xT = Taylor1([qaux, one(qaux)], orderT)
        tm = TaylorModel1(deepcopy(xT), x0, t00, dom)

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
            xξ = myrand(dom)
            q0ξ = t00 + myrand.(δq0)[1]
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
        tv = TaylorModel1(2, x0, ii0)
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
        tv = TaylorModel1(2, x1, ii1)
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
        tm = TaylorModel1(order, x0, ii0)

        integ_res = integrate(exp(tm), 1..1)
        exact_res = exp(tm)
        @test exact_res.pol == integ_res.pol
        @test issubset_interval(remainder(exact_res), remainder(integ_res))
        for ind = 1:_num_tests
            @test check_containment(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol == integ_res.pol
        @test issubset_interval(remainder(exact_res), remainder(integ_res))
        for ind = 1:_num_tests
            @test check_containment(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol == integ_res.pol
        @test issubset_interval(remainder(exact_res), remainder(integ_res))
        for ind = 1:_num_tests
            @test check_containment(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol == integ_res.pol
        @test_broken issubset_interval(remainder(exact_res), remainder(integ_res))
        for ind = 1:_num_tests
            @test check_containment(atan, integ_res)
        end
    end

    @testset "Display" begin
        tm = TaylorModel1(2, x1, ii1)
        use_show_default(true)
        @test string(exp(tm)) == "TaylorModel1{Interval{Float64}, Float64}" *
            "(Taylor1{Interval{Float64}}(Interval{Float64}" *
            "[Interval{Float64}(2.718281828459045, 2.7182818284590455, com), " *
            "Interval{Float64}(2.718281828459045, 2.7182818284590455, com), " *
            "Interval{Float64}(1.3591409142295225, 1.3591409142295228, com)], 2), " *
            "Interval{Float64}(-0.05020487208677604, 0.06448109909211741, trv), " *
            "Interval{Float64}(1.0, 1.0, com), Interval{Float64}(0.5, 1.5, com))"
        use_show_default(false)
        @test string(tm^3) == " Interval{Float64}(1.0, 1.0, com) + " *
            "Interval{Float64}(3.0, 3.0, com) t + " *
            "Interval{Float64}(3.0, 3.0, com) t² + " *
            "Interval{Float64}(-0.125, 0.125, trv)_NG"
        @test string(exp(tm)) ==
            " Interval{Float64}(2.718281828459045, 2.7182818284590455, com) + " *
            "Interval{Float64}(2.718281828459045, 2.7182818284590455, com) t + " *
            "Interval{Float64}(1.3591409142295225, 1.3591409142295228, com) t² + " *
            "Interval{Float64}(-0.05020487208677604, 0.06448109909211741, trv)"
    end

    # @testset "Tests for bounders" begin
    #     @testset "Tests for linear dominated bounder" begin
    #         order = 3

    #         f = x -> 1 + x^5 - x^4
    #         D = interval(0.9375, 1)
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_interval = f(D)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         @test diam(bound_ldb) <= diam(bound_interval)
    #         @test issubset_interval(bound_ldb, bound_naive_tm)

    #         D = interval(0.75, 0.8125)
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_interval = f(D)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         @test diam(bound_ldb) <= diam(bound_interval)
    #         @test issubset_interval(bound_ldb, bound_naive_tm)

    #         f = x -> x^2 * sin(x)
    #         D = interval(-1.875, -1.25)
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_interval = f(D)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         @test diam(bound_ldb) <= diam(bound_interval)
    #         @test issubset_interval(bound_ldb, bound_naive_tm)

    #         D = interval(1.25, 1.875)
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_interval = f(D)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         @test diam(bound_ldb) <= diam(bound_interval)
    #         @test issubset_interval(bound_ldb, bound_naive_tm)
    #     end

    #     @testset "Tests for quadratic fast bounder" begin
    #         order = 3

    #         f = x -> 1 + x^5 - x^4
    #         D = 0.75 .. 0.8125
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         bound_qfb = quadratic_fast_bounder(fT)
    #         @test bound_qfb ⊆ bound_naive_tm
    #         # @test diam(bound_qfb) <= diam(bound_ldb)

    #         f = x -> x^2 * sin(x)
    #         D = -2.5 .. -1.875
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         bound_qfb = quadratic_fast_bounder(fT)
    #         @test bound_qfb ⊆ bound_naive_tm
    #         @test diam(bound_qfb) <= diam(bound_ldb)

    #         f = x -> x^3 * cos(x) + x
    #         D = 3.75 .. 4.375
    #         x0 = mid(D)
    #         tm = TaylorModel1(order, x0, D)
    #         fT = f(tm)
    #         bound_naive_tm = fT(centered_dom(fT))
    #         bound_ldb = linear_dominated_bounder(fT)
    #         bound_qfb = quadratic_fast_bounder(fT)
    #         @test bound_qfb ⊆ bound_naive_tm
    #         # @test diam(bound_qfb) <= diam(bound_ldb)
    #     end
    # end
end
