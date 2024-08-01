# Tests using TaylorModel1

using TaylorModels
# using LinearAlgebra: norm
using Test

const _num_tests = 1000
const α_mid = TaylorModels.α_mid

setformat(:full)

function check_containmentTM1(ftest, tma::T) where {T<:Union{TaylorModel1, RTaylorModel1}}
    x0 = expansion_point(tma)
    xfp = diam(domain(tma))*(rand()-0.5) + mid(x0)
    xbf = big(xfp)
    range = tma((xfp .. xfp)-x0)
    bb = ftest(xbf) ∈ range
    bb || @show(ftest, xfp, xbf, ftest(xbf), range)
    return bb
end

@testset "Test `bound_taylor1`" begin
    x0 = Interval(0.0)
    ii0 = Interval(-0.5, 0.5)

    tpol = exp( Taylor1(2) )
    @test TaylorModels.bound_taylor1( tpol, ii0) == tpol(ii0.lo) .. tpol(ii0.hi)
    @test TaylorModels.bound_taylor1( exp( Taylor1(Interval{Float64}, 2) ),
        ii0) == tpol(ii0.lo) .. tpol(ii0.hi)

    # An uncomfortable example from Makino
    t = Taylor1(5)
    f(x) = 1 - x^4 + x^5
    @test interval(1-4^4/5^5,1) ⊆ TaylorModels.bound_taylor1(f(t), 0..1)
    tm = TaylorModel1(5, x0, ii0)
    @test interval(1-4^4/5^5,1) ⊆ TaylorModels.bound_taylor1(f(tm))
    @test interval(1-4^4/5^5,1) ⊆ TaylorModels.bound_taylor1(f(tm), 0..1)
end

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
        @test remainder(tv) == interval(0.0)
        @test polynomial(tv) == Taylor1(Interval{Float64},5)
        @test domain(tv) == ii0
        @test expansion_point(tv) == x0
        @test constant_term(tv) == interval(0.0)
        @test linear_polynomial(tv) == Taylor1(Interval{Float64},5)
        @test nonlinear_polynomial(tv) == zero(Taylor1(Interval{Float64},5))
        @test centered_dom(tv) == ii0
        @test centered_dom(TaylorModel1(5, 0.7, ii1)) == ii1-0.7

        # Tests related to fixorder
        a = TaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        b = TaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)
        aa, bb = TaylorModels.fixorder(a, b)
        @test get_order(aa) == get_order(bb) == 1
        @test isa(aa, TaylorModel1) == isa(bb, TaylorModel1) == true
        @test aa == a
        @test bb == TaylorModel1(Taylor1([1.0, 1, 0]), -1 .. 2, 0..0, -1 .. 1)
        # a and b remain the same
        @test a == TaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        @test b == TaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        a = TaylorModel1(x1+Taylor1(5), Δ, x1, ii1)
        tv = TaylorModel1(5, x1, ii1)
        a_pol = polynomial(a)
        tv_pol = polynomial(tv)
        @test zero(a) == TaylorModel1(zero(a_pol), 0..0, x1, ii1)
        @test one(a) == TaylorModel1(one(a_pol), 0..0, x1, ii1)
        @test a+x1 == TaylorModel1(2*x1+Taylor1(5), Δ, x1, ii1)
        @test a+a == TaylorModel1(2*(x1+Taylor1(5)), 2*Δ, x1, ii1)
        @test a-x1 == TaylorModel1(zero(x1)+Taylor1(5), Δ, x1, ii1)
        @test a-a == TaylorModel1(zero(a_pol), 2*Δ, x1, ii1)
        b = a * tv
        @test b == TaylorModel1(a_pol*tv_pol, remainder(a)*tv_pol(ii1-x1), x1, ii1)
        @test remainder(b/tv) ⊆ Interval(-0.78125, 0.84375)
        @test constant_term(b) == 1..1
        @test linear_polynomial(b) == 2*x1*Taylor1(5)
        @test nonlinear_polynomial(b) == x1*Taylor1(5)^2
        b = a * a_pol[0]
        @test b == a
        @test constant_term(a) == x1
        @test linear_polynomial(a) == Taylor1(5)
        @test nonlinear_polynomial(a) == Taylor1(0..0, 5)

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
        a = TaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        b = TaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)
        aa, bb = TaylorModels.fixorder(a, b)
        @test get_order(aa) == get_order(bb)
        @test get_order(bb) == 1
        @test a + b == aa + bb
        @test a - b == aa - bb
        res1 = a * b
        res2 = aa * bb
        @test res1 == TaylorModel1(Taylor1([1.0, 2]), -2 .. 9 , 0..0, -1 .. 1)
        @test res2 == TaylorModel1(Taylor1([1.0, 2]), -3 .. 9 , 0..0, -1 .. 1)
        res1 = a / b
        res2 = aa / bb
        @test res1 == TaylorModel1(Taylor1([1.0, 0]), entireinterval() , 0..0, -1 .. 1)
        @test res2 == res1
        # a and b remain the same
        @test a == TaylorModel1(Taylor1([1.0, 1]), 0..1, 0..0, -1 .. 1)
        @test b == TaylorModel1(Taylor1([1.0, 1, 0, 1]), 0..1, 0..0, -1 .. 1)

        @test_throws AssertionError a+TaylorModel1(a.pol, remainder(a), 1..1, -1..1)
        @test_throws AssertionError a+TaylorModel1(a.pol, remainder(a), 0..0, -2..2)
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
        δq0 = IntervalBox(-0.1 .. 0.1, Val(1))
        qaux = normalize_taylor(q0[1] + TaylorN(1, order=orderQ), δq0, true)
        symIbox = IntervalBox(-1 .. 1, Val(1))
        t = Taylor1([qaux, one(qaux)], orderT)
        dom = 0 .. 1
        x00 = mid(dom)

        f(x) = x + x^2
        g(x) = x
        h(x) = x^3*(1+x)

        tm = TaylorModel1(deepcopy(t), 0 .. 0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test isentire(remainder(fgTM1))

        fgTM1 = f(tm) * (g(tm))^2
        hh = h(tm)
        @test polynomial(fgTM1) ≈ polynomial(hh)
        @test remainder(fgTM1) == remainder(hh)
        for ind = 1:_num_tests
            xξ = rand(dom)-x00
            qξ = rand(symIbox)
            tt = t(xξ)(qξ)
            @test h(tt) ⊆ fgTM1(dom-x00)(symIbox)
        end

        t = Taylor1([one(qaux), qaux], orderT)
        tm = TaylorModel1(deepcopy(t), 0 .. 0, x00, dom)
        fgTM1 = f(tm) / g(tm)
        @test !isentire(remainder(fgTM1))
        for ind = 1:_num_tests
            xξ = rand(dom)-x00
            qξ = rand(symIbox)
            tt = 1+t(xξ)(qξ)
            @test tt ⊆ fgTM1(dom-x00)(symIbox)
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
        tm = TaylorModel1(deepcopy(t), 0 .. 0, x00, dom)   # order 4
        t8 = Taylor1([qaux, one(qaux)], 2*orderT)
        tm8 = TaylorModel1(deepcopy(t8), 0 .. 0, x00, dom)
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
            xξ = rand(dom)-x00
            qξ = rand(symIbox)
            tt = t(xξ)(qξ)
            @test exp(tt) ⊆ exp_tm4(dom-x00)(symIbox)
            @test exp(tt) ⊆ exp_tm(dom-x00)(symIbox)
        end
    end

    @testset "RPAs, functions and remainders" begin
        @test rpa(x->5+zero(x), TaylorModel1(4, x0, ii0)) ==
            TaylorModel1(interval(5.0), 4, x0, ii0)
        @test rpa(x->5+one(x), TaylorModel1(4, x1, ii1)) ==
            TaylorModel1(5+x1, 4, x1, ii1)
        @test rpa(x->5*x, TaylorModel1(4, x1, ii1)) ==
            TaylorModel1(5.0*(x1+Taylor1(4)), x0, x1, ii1)
        @test rpa(x->5*x^0, TaylorModel1(4, x0, ii0)) == 5*TaylorModel1(4, x0, ii0)^0
        @test rpa(x->5*x^0, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^0, x0, x0, ii0)
        @test rpa(x->5*x^1, TaylorModel1(4, x0, ii0)) == 5*TaylorModel1(4, x0, ii0)^1
        @test rpa(x->5*x^1, TaylorModel1(4, x0, ii0)) ==
            TaylorModel1( interval(5.0)*Taylor1(4)^1, x0, x0, ii0)
        @test rpa(x->5*x^2, TaylorModel1(4, x0, ii0)) == 5*TaylorModel1(4, x0, ii0)^2
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
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containmentTM1(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

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
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containmentTM1(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

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
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containmentTM1(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

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
        @test interval(ftest(ii.lo)-tmc.pol(ii.lo-ξ0),
                        ftest(ii.hi)-tmc.pol(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containmentTM1(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

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
        @test interval(ftest(ii.hi)-tmc.pol(ii.hi-ξ0),
                        ftest(ii.lo)-tmc.pol(ii.lo-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containmentTM1(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))

        # Example of Makino's thesis (page 98 and fig 4.2)
        order = 8
        ii = interval(-0.5, 1.0)
        xx = interval(mid(ii))
        ftest = x -> x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
        tm = TaylorModel1(order, xx, ii)
        tma = rpa(ftest, tm)
        tmb = ftest(tm)
        @test remainder(tmb) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containmentTM1(ftest, tma)
        end
        @test_throws AssertionError tmb(ii.hi+1.0)
        @test_throws AssertionError tmb(ii+Interval(1))
    end

    @testset "RPAs with polynomial Taylor1{TaylorN{T}}" begin
        orderT = 5
        orderQ = 5
        dom = 0 .. 1
        x00 = mid(dom)
        q0 = [0.5]
        symIbox = IntervalBox(-1 .. 1, 1)
        δq0 = IntervalBox(-0.2 .. 0.2, 1)
        qaux = normalize_taylor(TaylorN(1, order=orderQ) + q0[1], δq0, true)
        xT = Taylor1([qaux, one(qaux)], orderT)
        tm = TaylorModel1(deepcopy(xT), 0 .. 0, x00, dom)

        f(x) = sin(x)
        ff(x) = cos(x)
        g(x) = exp(x)
        gg(x) = x^5
        h(x) = log(1.0+x)
        hh(x) = x^3 / x^5

        fT = f(tm)
        ffT = ff(tm)
        gT = g(tm)
        ggT = gg(tm)
        hT = h(tm)
        hhT = hh(tm)

        for ind = 1:_num_tests
            xξ = rand(domain(fT))
            q0ξ = (q0 .+ rand(δq0))[1]
            t = Taylor1(orderT) + q0ξ
            ft = f(t)
            fft = ff(t)
            gt = g(t)
            ggt = gg(t)
            ht = h(t)
            hht = hh(t)

            @test ft(xξ - q0ξ) ⊆ fT(xξ - fT.x0)(symIbox)
            @test fft(xξ - q0ξ) ⊆ ffT(xξ - ffT.x0)(symIbox)
            @test gt(xξ - q0ξ) ⊆ gT(xξ - gT.x0)(symIbox)
            @test ggt(xξ - q0ξ) ⊆ ggT(xξ - ggT.x0)(symIbox)
            @test ht(xξ - q0ξ) ⊆ hT(xξ - hT.x0)(symIbox)
            @test hht(xξ - q0ξ) ⊆ hhT(xξ - hhT.x0)(symIbox)
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
        @test remainder(exact_res) ⊆ remainder(integ_res)
        for ind = 1:_num_tests
            @test check_containmentTM1(exp, integ_res)
        end

        integ_res = integrate(cos(tm))
        exact_res = sin(tm)
        @test exact_res.pol == integ_res.pol
        @test remainder(exact_res) ⊆ remainder(integ_res)
        for ind = 1:_num_tests
            @test check_containmentTM1(sin, integ_res)
        end

        integ_res = integrate(-sin(tm), 1..1)
        exact_res = cos(tm)
        @test exact_res.pol == integ_res.pol
        @test remainder(exact_res) ⊆ remainder(integ_res)
        for ind = 1:_num_tests
            @test check_containmentTM1(cos, integ_res)
        end

        integ_res = integrate(1/(1+tm^2))
        exact_res = atan(tm)
        @test exact_res.pol == integ_res.pol
        # @test remainder(exact_res) ⊆ remainder(integ_res)
        for ind = 1:_num_tests
            @test check_containmentTM1(atan, integ_res)
        end
    end

    @testset "Display" begin
        tm = TaylorModel1(2, x1, ii1)
        use_show_default(true)
        if VERSION < v"1.6"
            @test string(exp(tm)) == "TaylorModel1{Interval{Float64},Float64}" *
                "(Taylor1{Interval{Float64}}(Interval{Float64}" *
                "[Interval(2.718281828459045, 2.7182818284590455), Interval(2.718281828459045, 2.7182818284590455), " *
                "Interval(1.3591409142295225, 1.3591409142295228)], 2), Interval(-0.05020487208677604, 0.06448109909211741), " *
                "Interval(1.0, 1.0), Interval(0.5, 1.5))"
        else
            @test string(exp(tm)) == "TaylorModel1{Interval{Float64}, Float64}" *
                "(Taylor1{Interval{Float64}}(Interval{Float64}" *
                "[Interval(2.718281828459045, 2.7182818284590455), Interval(2.718281828459045, 2.7182818284590455), " *
                "Interval(1.3591409142295225, 1.3591409142295228)], 2), Interval(-0.05020487208677604, 0.06448109909211741), " *
                "Interval(1.0, 1.0), Interval(0.5, 1.5))"
        end
        use_show_default(false)
        @test string(tm^3) == " Interval(1.0, 1.0) + Interval(3.0, 3.0) t + " *
            "Interval(3.0, 3.0) t² + Interval(-0.125, 0.125)"
        @test string(exp(tm)) == " Interval(2.718281828459045, 2.7182818284590455) + " *
            "Interval(2.718281828459045, 2.7182818284590455) t + " *
            "Interval(1.3591409142295225, 1.3591409142295228) t² + " *
            "Interval(-0.05020487208677604, 0.06448109909211741)"
    end

    @testset "Tests for bounders" begin
        @testset "Tests for linear dominated bounder" begin
            order = 3

            f = x -> 1 + x^5 - x^4
            D = 0.9375 .. 1
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_interval = f(D)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            @test diam(bound_ldb) <= diam(bound_interval)
            @test bound_ldb ⊆ bound_naive_tm

            D = 0.75 .. 0.8125
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_interval = f(D)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            @test diam(bound_ldb) <= diam(bound_interval)
            @test bound_ldb ⊆ bound_naive_tm

            f = x -> x^2 * sin(x)
            D = -1.875 .. -1.25
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_interval = f(D)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            @test diam(bound_ldb) <= diam(bound_interval)
            @test bound_ldb ⊆ bound_naive_tm

            D = 1.25 .. 1.875
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_interval = f(D)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            @test diam(bound_ldb) <= diam(bound_interval)
            @test bound_ldb ⊆ bound_naive_tm
        end

        @testset "Tests for quadratic fast bounder" begin
            order = 3

            f = x -> 1 + x^5 - x^4
            D = 0.75 .. 0.8125
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            bound_qfb = quadratic_fast_bounder(fT)
            @test bound_qfb ⊆ bound_naive_tm
            # @test diam(bound_qfb) <= diam(bound_ldb)

            f = x -> x^2 * sin(x)
            D = -2.5 .. -1.875
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            bound_qfb = quadratic_fast_bounder(fT)
            @test bound_qfb ⊆ bound_naive_tm
            @test diam(bound_qfb) <= diam(bound_ldb)

            f = x -> x^3 * cos(x) + x
            D = 3.75 .. 4.375
            x0 = mid(D)
            tm = TaylorModel1(order, x0, D)
            fT = f(tm)
            bound_naive_tm = fT(centered_dom(fT))
            bound_ldb = linear_dominated_bounder(fT)
            bound_qfb = quadratic_fast_bounder(fT)
            @test bound_qfb ⊆ bound_naive_tm
            # @test diam(bound_qfb) <= diam(bound_ldb)
        end
    end
end
