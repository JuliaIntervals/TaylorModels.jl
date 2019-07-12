# Tests using TaylorModel1 and RTaylorModel1

using TaylorModels
using LinearAlgebra: norm
using Test

const _num_tests = 1000

setformat(:full)


function check_containment(ftest, xx::TaylorModelN{N,T,S}, tma::TaylorModelN{N,T,S}) where {N,T,S}
    xfp = diam.(tma.dom) .* (rand(N) .- 0.5) .+ mid(tma.x0)
    xbf = [big(xfp[i]) for i=1:N]
    ib = IntervalBox([@interval(xfp[i]) for i=1:N]...)
    range = evaluate(tma, ib-tma.x0)
    bb = all(ftest(xx(xbf .- mid(tma.x0))) ⊆ range)
    bb || @show(ftest, ib, xbf, ftest(xbf...), range)
    return bb
end

const _order = 2
const _order_max = 2*(_order+1)
set_variables(Interval{Float64}, [:x, :y], order=_order_max)

@testset "Tests for TaylorModelN " begin
    b0 = Interval(0.0) × Interval(0.0)
    ib0 = Interval(-0.5, 0.5) × Interval(-0.5, 0.5)
    b1 =  Interval(0.0) × Interval(1.0)
    ib1 = Interval(-0.5, 0.5) × Interval(0.5, 1.5)

    zi = 0..0
    xT = TaylorN(Interval{Float64}, 1, order=_order)
    yT = TaylorN(Interval{Float64}, 2, order=_order)

    @testset "TaylorModelN constructors" begin
        xm = TaylorModelN{2, Interval{Float64}, Float64}(xT, zi, b0, ib0)
        ym = TaylorModelN{2, Interval{Float64}, Float64}(yT, zi, b0, ib0)
        @test xm == TaylorModelN(xT, zi, b0, ib0)
        @test ym == TaylorModelN(yT, zi, b0, ib0)
        @test xm == TaylorModelN(1, _order, b0, ib0)
        @test ym == TaylorModelN(2, _order, b0, ib0)
        @test TaylorModelN( b1[1], 2, b0, ib0) ==
                TaylorModelN(TaylorN(b1[1], _order), zi, b0, ib0)

        @test isa(xm, AbstractSeries)
        @test TaylorModelN{2, Interval{Float64},Float64} <: AbstractSeries{Interval{Float64}}

        # Test errors in construction
        @test_throws AssertionError TaylorModelN(xT, zi, IntervalBox(1..1), IntervalBox(1..1))
        @test_throws AssertionError TaylorModelN(xT, zi, b0, ib1)
        @test_throws AssertionError TaylorModelN(xT, 1..1, b0, ib0)
        @test_throws BoundsError TaylorModelN(5, _order, b0, ib0) # wrong variable number

        # Tests for get_order and remainder
        @test get_order() == 6
        @test get_order(xm) == 2
        @test remainder(ym) == zi
        @test constant_term(xm) == interval(0.0)
        @test constant_term(ym) == interval(0.0)
        @test linear_polynomial(xm) == xT
        @test linear_polynomial(ym) == yT
        @test linear_polynomial(xm^2) == zero(xT)
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        xm = TaylorModelN(xT, zi, b1, ib1)
        ym = TaylorModelN(yT, zi, b1, ib1)
        a = TaylorModelN( b1[1]+xT, Δ, b1, ib1)
        @test zero(a) == TaylorModelN(zero(a.pol), 0..0, b1, ib1)
        @test one(a) == TaylorModelN(one(a.pol), 0..0, b1, ib1)
        @test a + a == TaylorModelN(2*(b1[1]+ xT), 2*Δ, b1, ib1)
        @test -a == TaylorModelN( -(b1[1]+xT), -Δ, b1, ib1)
        @test a - a == TaylorModelN(zero(a.pol), 2*Δ, b1, ib1)
        @test b1[2] + ym == TaylorModelN(b1[2] + yT, zi, b1, ib1)
        @test a - b1[1] == TaylorModelN(zero(b1[1])+xT, Δ, b1, ib1)
        @test constant_term(a) == b1[1]
        @test constant_term(b1[2] + ym) == b1[2]
        @test linear_polynomial(a) == xT
        @test linear_polynomial(b1[2] + ym) == yT
        @test linear_polynomial(ym^2) == zero(yT)

        @test xm * ym == TaylorModelN( xT * yT, zi, b1, ib1)
        b = a * ym
        @test b == TaylorModelN( xT * yT, Δ*ym.pol(ib1-b1), b1, ib1)
        b = ym * TaylorModelN(xT^2, zi, b1, ib1)
        @test b == TaylorModelN( zero(xT), (ib1[1]-b1[1])^2*(ib1[2]-b1[2]), b1, ib1 )
        b = b1[2] * a
        @test b == TaylorModelN( b1[2]*a.pol, Δ*b1[2], b1, ib1 )
        @test b / b1[2] == a
        @test_throws AssertionError TaylorModelN(TaylorN(1, order=_order_max), zi, b1, ib1) *
            TaylorModelN(TaylorN(2, order=_order_max), zi, b1, ib1)

        remt = remainder(1/(1-TaylorModel1(_order, b1[1], ib1[1])))
        @test remainder(1 / (1-xm)) == remt
        @test remainder(ym / (1-xm)) == Interval(-0.25, 0.25)

        @test remainder(xm^2) == remainder(ym^2)
        @test (xm.dom[1]-xm.x0[1])^3 == remainder(xm^3)
        @test (ym.dom[2]-ym.x0[2])^4 ⊆ remainder(ym^4)
    end

    @testset "RPAs, functions and remainders" begin
        xm = TaylorModelN(1, _order, b1, ib1)
        ym = TaylorModelN(2, _order, b1, ib1)

        @test rpa(x->5+zero(x), xm) == 5+zero(xm)
        @test rpa(x->5+one(x), ym) == 5+one(ym)
        @test rpa(x->5*x, ym) == 5*ym
        remT = remainder(5*TaylorModel1(2, b1[1], ib1[1])^4)
        @test rpa(x->5*x^4, xm) == TaylorModelN(zero(xT), remT, b1, ib1)
        remT = 5 * (ib1[1]-b1[1])^2 * (2*(ib1[2]-b1[2])+(ib1[2]-b1[2])^2)
        @test rpa(x->5*x^2, xm*ym) == TaylorModelN( 5*xT^2, remT, b1, ib1)

        # Testing remainders of an RPA
        ftest = x -> exp(x)-1
        xx = xm + ym
        tma = rpa(ftest, xx)
        tmb = ftest(xx)
        @test tma == tmb
        # fT, Δ, ξ0 = rpafp(tma)
        # @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
        #                 ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, xx, tma)
        end

        # test for TM with scalar coefficients
        tmc = fp_rpa(tma)
        @test fp_rpa(tmc) == tmc

        ftest = x -> exp(x)
        xx = xm + 2*ym
        tma = rpa(ftest, xx)
        tmb = ftest(xx)
        @test tma == tmb
        # fT, Δ, ξ0 = rpafp(tma)
        # @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
        #                 ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, xx, tma)
        end

        ftest = x -> sin(x)
        xx = 1- xm^2 + ym
        tma = rpa(ftest, xx)
        tmb = ftest(xx)
        @test tma == tmb
        # fT, Δ, ξ0 = rpafp(tma)
        # @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
        #                 ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, xx, tma)
        end

        ftest = x -> sqrt(x)
        xx = xm^2 + ym
        tma = rpa(ftest, xx)
        tmb = ftest(xx)
        @test tma == tmb
        # fT, Δ, ξ0 = rpafp(tma)
        # @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
        #                 ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, xx, tma)
        end

        ftest = x -> inv(x)
        xx = 1 + xm + ym
        tma = rpa(ftest, xx)
        tmb = ftest(xx)
        @test tma == tmb
        # fT, Δ, ξ0 = rpafp(tma)
        # @test interval(ftest(ii.lo)-fT(ii.lo-ξ0),
        #                 ftest(ii.hi)-fT(ii.hi-ξ0)) ⊆ remainder(tma)
        for ind = 1:_num_tests
            @test check_containment(ftest, xx, tma)
        end
    end

    @testset "Composition of functions and their inverses" begin
        xm = TaylorModelN(xT, zi, b1, ib1)
        ym = TaylorModelN(yT, zi, b1, ib1)

        tma = exp(ym)
        tmb = log(tma)
        @test tmb == log(exp(ym))
        @test tmb.pol == ym.pol

        tma = sin(xm)
        tmb = asin(tma)
        @test tmb == asin(sin(xm))
        @test tmb.pol == xm.pol

        tma = asin(xm)
        tmb = sin(tma)
        @test tmb == sin(asin(xm))
        @test tmb.pol == xm.pol

        tma = acos(ym)
        tmb = cos(tma)
        @test tmb == cos(acos(ym))
        @test sup(norm(tmb.pol - ym.pol, Inf)) < 1.0e-16

        tma = tan(xm)
        tmb = atan(tma)
        @test tmb == atan(tan(xm))
        @test tmb.pol == xm.pol

        tma = atan(xm)
        tmb = tan(tma)
        @test tmb == tan(atan(xm))
        @test tmb.pol == xm.pol


        ####
        tma = log(1+ym)
        tmb = exp(tma)
        @test tmb == exp(log(1+ym))
        @test tmb.pol == (1+ym).pol

        tma = cos(1+ym)
        tmb = acos(tma)
        @test tmb == acos(cos(1+ym))
        @test sup(norm(tmb.pol - (1+ym).pol, Inf)) < 1.0e-15
    end

    @testset "Display" begin
        xm = TaylorModelN(1, _order, b1, ib1)
        ym = TaylorModelN(2, _order, b1, ib1)
        use_show_default(true)
        @test string(xm+ym) == "TaylorModelN{2,Interval{Float64},Float64}" *
            "(TaylorN{Interval{Float64}}" *
            "(HomogeneousPolynomial{Interval{Float64}}" *
            "[HomogeneousPolynomial{Interval{Float64}}" *
            "(Interval{Float64}[Interval(1.0, 1.0)], 0), " *
            "HomogeneousPolynomial{Interval{Float64}}" *
            "(Interval{Float64}[Interval(1.0, 1.0), Interval(1.0, 1.0)], 1), " *
            "HomogeneousPolynomial{Interval{Float64}}" *
            "(Interval{Float64}[Interval(0.0, 0.0), Interval(0.0, 0.0), " *
            "Interval(0.0, 0.0)], 2)], 2), Interval(0.0, 0.0), IntervalBox(Interval(0.0, 0.0), " *
            "Interval(1.0, 1.0)), IntervalBox(Interval(-0.5, 0.5), Interval(0.5, 1.5)))"
        use_show_default(false)
        @test string((xm+ym)^2) == " Interval(1.0, 1.0) + Interval(2.0, 2.0) x + " *
            "Interval(2.0, 2.0) y + Interval(1.0, 1.0) x² + Interval(2.0, 2.0) x y + " *
            "Interval(1.0, 1.0) y² + Interval(0.0, 0.0)"
    end

end
