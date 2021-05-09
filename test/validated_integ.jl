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
        δq0 = IntervalBox(-0.25 .. 0.25, 2)
        X0 = IntervalBox(q0 .+ δq0)

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))
        normalized_box = IntervalBox(-1 .. 1, Val(orderQ))

        tTM, qv, qTM = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol)

        @test length(qv) == length(qTM[1, :]) == length(tTM)
        @test length(tTM) < 501

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1,n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:,n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, tini, q0 .+ q0ξ) .∈ q)
            end
        end

        tTMf, qvf, qTMf = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol,
            adaptive=false)
        @test length(qvf) == length(qv)
        @test all(qTM .== qTMf)

        # initializaton with a Taylor model
        X0tm = qTM[:, 1]
        @assert X0tm isa Vector{TaylorModel1{TaylorN{Float64}, Float64}}
        tTM2, qv2, qTM2 = validated_integ(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol)
        @test all(iszero, (qTM - qTM2))
        
        @taylorize function x_square!(dx, x, p, t)
            dx[1] = x[1]^2
            nothing
        end

        exactsol(t, x0) = x0[1] / (1 - x0[1] * t)

        tini, tend = 0., 0.45
        abstol = 1e-15
        orderQ = 5
        orderT = 20
        q0 = [2.]
        δq0 = IntervalBox(-0.1 .. 0.1, Val(1))
        X0 = IntervalBox(q0 .+ δq0)
        ξ = set_variables("ξₓ", numvars=1, order=2*orderQ)
        
        tTM, qv, qTM = validated_integ(x_square!, X0, tini, tend, orderQ, orderT, abstol)

        @test length(qv) == length(qTM[1, :]) == length(tTM)
        @test length(tTM) < 501
        normalized_box = IntervalBox(-1 .. 1, Val(1))

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1, n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:, n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, q0 .+ q0ξ) .∈ q)
            end
        end

        tTMf, qvf, qTMf = validated_integ(x_square!, X0, tini, tend, orderQ, orderT, abstol,
            adaptive=false)
        @test length(qvf) == length(qv)
        @test all(qTMf .== qTM)
    end

    @testset "Forward integration for validated_integ2" begin
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
        δq0 = IntervalBox(-0.25 .. 0.25, 2)
        X0 = IntervalBox(q0 .+ δq0)

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))
        normalized_box = IntervalBox(-1 .. 1, Val(orderQ))

        tTM, qv, qTM = validated_integ2(falling_ball!, X0,
            tini, tend, orderQ, orderT, abstol)

        @test length(qv) == length(qTM[1, :]) == length(tTM)
        @test length(tTM) < 501

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1,n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:,n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, tini, q0 .+ q0ξ) .∈ q)
            end
        end

        @taylorize function x_square(dx, x, p, t)
            dx[1] = x[1]^2
            nothing
        end

        exactsol(t, x0) = x0[1] / (1 - x0[1] * t)

        tini, tend = 0., 0.45
        abstol = 1e-15
        orderQ = 5
        orderT = 20
        q0 = [2.]
        δq0 = IntervalBox(-0.1 .. 0.1, Val(1))
        X0 = IntervalBox(q0 .+ δq0)
        ξ = set_variables("ξₓ", numvars=1, order=2*orderQ)
        tTM, qv, qTM = validated_integ2(x_square, X0, tini, tend, orderQ, orderT, abstol)

        @test length(qv) == length(qTM[1, :]) == length(tTM)
        @test length(tTM) < 501
        normalized_box = IntervalBox(-1 .. 1, Val(1))

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1, n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:, n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, q0 .+ q0ξ) .∈ q)
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
        δq0 = IntervalBox(-0.25 .. 0.25, 2)
        X0 = IntervalBox(q0 .+ δq0)

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))
        normalized_box = IntervalBox(-1 .. 1, Val(orderQ))

        tTM, qv, qTM = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol)

        @test length(qv) == length(qTM[1, :]) == length(tTM)
        @test length(tTM) < 501

        for n = 2:length(tTM)
            for it = 1:10
                δt = interval_rand(domain(qTM[1,n]))
                q0ξ = interval_rand(δq0)
                q = evaluate.(evaluate.(qTM[:,n], δt), (normalized_box,))
                @test all(exactsol(tTM[n-1]+δt, tini, q0 .+ q0ξ) .∈ q)
            end
        end

        tTMf, qvf, qTMf = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol,
            adaptive=false)
        @test length(qvf) == length(qv)
        @test all(qTM .== qTMf)

        # initializaton with a Taylor model
        X0tm = qTM[:, 1]
        @assert X0tm isa Vector{TaylorModel1{TaylorN{Float64}, Float64}}
        tTM2, qv2, qTM2 = validated_integ(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol)
        @test all(iszero, (qTM - qTM2))

        tTM2f, qv2f, qTM2f = validated_integ(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol,
            adaptive=false)
        @test length(qv2f) == length(qv2)
        @test all(qTM2 .== qTM2f)
    end

    @testset "Backward integration for validated_integ2" begin
        @taylorize function falling_ball!(dx, x, p, t)
            dx[1] = x[2]
            dx[2] = -one(x[1])
            nothing
        end
        exactsol(t, t0, x0) = (x0[1] + x0[2]*(t-t0) - 0.5*(t-t0)^2, x0[2] - (t-t0))

        # Initial conditions
        tini, tend = 10.0, 0.0
        q0 = [10.0, 0.0]
        δq0 = IntervalBox(-0.25 .. 0.25, 2)
        X0 = IntervalBox(q0 .+ δq0)

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))
        normalized_box = IntervalBox(-1 .. 1, Val(orderQ))

        tTM, qv, qTM = validated_integ2(falling_ball!, X0,
                                        tini, tend, orderQ, orderT, abstol)

        @test length(qv) == length(qTM[1, :]) == length(tTM)
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

    @testset "Pendulum with constant torque" begin
        @taylorize function pendulum!(dx, x, p, t)
            si = sin(x[1])
            aux = 2 *  si
            dx[1] = x[2]
            dx[2] = aux + 8*x[3]
            dx[3] = zero(x[1])
            nothing
        end
        # Conserved quantity
        ene_pendulum(x) = x[2]^2/2 + 2 * cos(x[1]) - 8 * x[3]

        # Initial conditions
        tini, tend = 0.0, 12.0
        q0 = [1.1, 0.1, 0.0]
        δq0 = IntervalBox(-0.1 .. 0.1, -0.1 .. 0.1, 0..0)
        X0 = IntervalBox(q0 .+ δq0)
        ene0 = ene_pendulum(X0)

        # Parameters
        abstol = 1e-10
        orderQ = 3
        orderT = 10
        ξ = set_variables("ξ", order=2*orderQ, numvars=length(q0))

        tTM, qv, qTM = validated_integ(pendulum!, X0, tini, tend, orderQ, orderT, abstol,
            maxsteps=1800);
        @test all(ene0 .⊆ ene_pendulum.(qv))

        tTM, qv, qTM = validated_integ2(pendulum!, X0, tini, tend, orderQ, orderT, abstol,
            maxsteps=1800);
        @test all(ene0 .⊆ ene_pendulum.(qv))

        # Initial conditions 2
        q0 = [1.1, 0.1, 0.0]
        δq0 = IntervalBox(-0.1 .. 0.1, -0.1 .. 0.1, -0.01 .. 0.01)
        X0 = IntervalBox(q0 .+ δq0)
        ene0 = ene_pendulum(X0)

        tTM, qv, qTM = validated_integ(pendulum!, X0, tini, tend, orderQ, orderT, abstol,
            maxsteps=1800);
        @test all(ene0 .⊆ ene_pendulum.(qv))

        tTM, qv, qTM = validated_integ2(pendulum!, X0, tini, tend, orderQ, orderT, abstol,
            maxsteps=1800);
        @test all(ene0 .⊆ ene_pendulum.(qv))
    end
end
