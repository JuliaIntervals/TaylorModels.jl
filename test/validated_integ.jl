# Tests for validated_integ

using TaylorModels
using LinearAlgebra: norm
using Test

setformat(:full)

# NOTE: IntervalArithmetic v0.16.0 includes this function; but
# IntervalRootFinding is bounded to use v0.15.x
interval_rand(X::Interval{T}) where {T} = X.lo + rand(T) * (X.hi - X.lo)
interval_rand(X::IntervalBox) = interval_rand.(X)

@testset "Tests for `validated_integ`" begin
    @testset "Forward integration" begin
        @taylorize function falling_ball!(dx, x, p, t)
            dx[1] = x[2]
            dx[2] = -one(x[1])
            nothing
        end
        exactsol(t, t0, x0) =
            (x0[1] + x0[2]*(t-t0) - 0.5*(t-t0)^2, x0[2] - (t-t0))

        # Initial conditions
        tini, tend = 0.0, 10.0
        q0 = [10.0, 0.0]
        δq0 = IntervalBox(-0.25 .. 0.25, Val(2))

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))
        normalized_box = IntervalBox(-1 .. 1, Val(orderQ))

        tTM, qTM, qq = validated_integ(falling_ball!, q0, δq0,
            tini, tend, orderQ, orderT, abstol)

        @test size(qq) == size(qTM)
        @test size(qq)[2] == size(qTM)[2] == length(tTM)
        @test length(tTM) < 501

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1,n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:,n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, tini, q0 .+ q0ξ) .∈ q)
            end
        end
    end

    @testset "Backward integration" begin
        @taylorize function falling_ball!(dx, x, p, t)
            dx[1] = x[2]
            dx[2] = -one(x[1])
            nothing
        end
        exactsol(t, t0, x0) = (x0[1] + x0[2]*(t-t0) - 0.5*(t-t0)^2, x0[2] - (t-t0))

        # Initial conditions
        tini, tend = 10.0, 0.0
        q0 = [10.0, 0.0]
        δq0 = IntervalBox(-0.25 .. 0.25, Val(2))

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))
        normalized_box = IntervalBox(-1 .. 1, Val(orderQ))

        tTM, qTM, qq = validated_integ(falling_ball!, q0, δq0,
            tini, tend, orderQ, orderT, abstol)

        @test size(qq) == size(qTM)
        @test size(qq)[2] == size(qTM)[2] == length(tTM)
        @test length(tTM) < 501

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1,n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:,n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, tini, q0 .+ q0ξ) .∈ q)
            end
        end
    end
end
