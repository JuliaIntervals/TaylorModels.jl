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

    @testset "TMAbsRem constructors" begin
        tv = TMAbsRem(Taylor1(Interval{Float64},5), interval(0.0),
            interval(0.0), interval(-0.5, 0.5))
        @test tv == TMAbsRem{Float64}(Taylor1(Interval{Float64},5), interval(0.0),
            interval(0.0), interval(-0.5, 0.5))
        @test tv == TMAbsRem(5, interval(0.0), interval(-0.5, 0.5))
        @test TMAbsRem(interval(1.0), 5, interval(0.0), interval(-0.5, 0.5)) ==
            TMAbsRem{Float64}(Taylor1(interval(1.0), 5), interval(0.0),
                interval(0.0), interval(-0.5, 0.5))

        # Test errors in construction
        @test_throws AssertionError TMAbsRem(Taylor1(Interval{Float64},5),
            interval(1.0), interval(0.0), interval(-0.5, 0.5))
        @test_throws AssertionError TMAbsRem(5, interval(1.0), interval(-0.5, 0.5))

        # Tests for get_order and remainder
        @test get_order(tv) == 5
        @test remainder(tv) == interval(0.0)
    end

    @testset "TMAbsRem bounds" begin
        tpol = exp( Taylor1(2) )
        @test TaylorModels.bound_taylor1( tpol, interval(-0.5, 0.5)) ==
            @interval(tpol(-0.5), tpol(0.5))
        @test TaylorModels.bound_taylor1( exp( Taylor1(Interval{Float64}, 2) ),
            interval(-0.5, 0.5)) == @interval(tpol(-0.5), tpol(0.5))
        # An uncomfortable example from Makino
        t = Taylor1(5)
        TaylorModels.bound_taylor1(1.0-t^4+t^5, interval(0,1)) ==
            Interval(0.9180799999999999, 1.0)
    end

    #

end
