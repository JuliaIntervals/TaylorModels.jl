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

const _order = 2
const _order_max = 2*(_order+1)
set_variables(Interval{Float64}, [:x, :y], order=_order_max)

@testset "Tests for TMNAbsRem " begin
    b0 = Interval(0.0) × Interval(0.0)
    ib0 = Interval(-0.5, 0.5) × Interval(-0.5, 0.5)
    b1 =  Interval(0.0) × Interval(1.0)
    ib1 = Interval(-0.5, 0.5) × Interval(0.5, 1.5)

    zi = 0..0
    xT = TaylorN(Interval{Float64}, 1, order=_order)
    yT = TaylorN(Interval{Float64}, 2, order=_order)

    @testset "TMNAbsRem constructors" begin
        xm = TMNAbsRem{2, Interval{Float64}, Float64}(xT, zi, b0, ib0)
        ym = TMNAbsRem{2, Interval{Float64}, Float64}(yT, zi, b0, ib0)
        @test xm == TMNAbsRem(xT, zi, b0, ib0)
        @test ym == TMNAbsRem(yT, zi, b0, ib0)
        @test xm == TMNAbsRem(1, _order, b0, ib0)
        @test ym == TMNAbsRem(2, _order, b0, ib0)
        @test TMNAbsRem( b1[1], 2, b0, ib0) ==
                TMNAbsRem(TaylorN(b1[1], _order), zi, b0, ib0)

        # Test errors in construction
        @test_throws AssertionError TMNAbsRem(xT, zi, IntervalBox(1..1), IntervalBox(1..1))
        @test_throws AssertionError TMNAbsRem(xT, zi, b0, ib1)
        @test_throws AssertionError TMNAbsRem(xT, 1..1, b0, ib0)
        @test_throws BoundsError TMNAbsRem(5, _order, b0, ib0) # wrong variable number

        # Tests for get_order and remainder
        @test get_order() == 6
        @test get_order(xm) == 2
        @test remainder(ym) == zi
    end

    @testset "Arithmetic operations" begin
        Δ = interval(-0.25, 0.25)
        xm = TMNAbsRem(xT, zi, b1, ib1)
        ym = TMNAbsRem(yT, zi, b1, ib1)
        a = TMNAbsRem( b1[1]+xT, Δ, b1, ib1)
        @test a + a == TMNAbsRem(2*(b1[1]+ xT), 2*Δ, b1, ib1)
        @test -a == TMNAbsRem( -(b1[1]+xT), -Δ, b1, ib1)
        @test a - a == TMNAbsRem(zero(a.pol), 2*Δ, b1, ib1)
        @test b1[2] + ym == TMNAbsRem(b1[2] + yT, zi, b1, ib1)
        @test a - b1[1] == TMNAbsRem(zero(b1[1])+xT, Δ, b1, ib1)

        @test xm * ym == TMNAbsRem( xT * yT, zi, b1, ib1)
        b = a * ym
        @test b == TMNAbsRem( xT * yT, Δ*ym.pol(ib1-b1), b1, ib1)
        b = ym * TMNAbsRem(xT^2, zi, b1, ib1)
        @test b == TMNAbsRem( zero(xT), (ib1[1]-b1[1])^2*(ib1[2]-b1[2]), b1, ib1 )
        b = b1[2] * a
        @test b == TMNAbsRem( b1[2]*a.pol, Δ*b1[2], b1, ib1 )
        @test b / b1[2] == a
        @test_throws AssertionError TMNAbsRem(TaylorN(1, order=_order_max), zi, b1, ib1) *
            TMNAbsRem(TaylorN(2, order=_order_max), zi, b1, ib1)

        remt = remainder(1/(1-TM1AbsRem(_order, b1[1], ib1[1])))
        @test remainder(1 / (1-xm)) == remt
        @test remainder(ym / (1-xm)) == Interval(-0.25, 0.25)

        @test remainder(xm^2) == remainder(ym^2)
        @test (xm.iI[1]-xm.x0[1])^3 == remainder(xm^3)
        @test (ym.iI[2]-ym.x0[2])^4 ⊆ remainder(ym^4)
    end

    @testset "RPAs, functions and remainders" begin
        xm = TMNAbsRem(xT, zi, b1, ib1)
        ym = TMNAbsRem(yT, zi, b1, ib1)

        @test rpa(x->5+zero(x), xm) == 5+zero(xm)
        @test rpa(x->5*x, ym) == 5*ym
        remT = remainder(5*TM1AbsRem(2, b1[1], ib1[1])^4)
        @test rpa(x->5*x^4, xm) == TMNAbsRem(zero(xT), remT, b1, ib1)
        @test rpa(x->5*x^2, xm*ym) ==
                TMNAbsRem( zero(xT), 5*(ib1[1]-b1[1])^4, b1, ib1)
    end

    @testset "Composition of functions and their inverses" begin
        xm = TMNAbsRem(xT, zi, b1, ib1)
        ym = TMNAbsRem(yT, zi, b1, ib1)

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
end
