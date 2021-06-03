# Tests for validated_integ

using TaylorModels
# using LinearAlgebra: norm
using Test
using Random

const _num_tests = 1_000

setformat(:full)

# NOTE: IntervalArithmetic v0.16.0 includes this function; but
# IntervalRootFinding is bounded to use v0.15.x
interval_rand(X::Interval{T}) where {T} = X.lo + rand(T) * (X.hi - X.lo)
interval_rand(X::IntervalBox) = interval_rand.(X)

# Function to check, against an exact solution of the ODE, the computed
# validted solution
function test_integ(fexact, t0, qTM, q0, δq0)
    normalized_box = symmetric_box(length(q0), Float64)
    # Time domain
    domt = domain(qTM[1])
    # Random time (within time domain) and random initial condition
    δt = rand(domt)
    δtI = (δt .. δt) ∩ domt
    q0ξ = interval_rand(δq0)
    q0ξB = IntervalBox([(q0ξ[i] .. q0ξ[i]) ∩ δq0[i] for i in eachindex(q0ξ)])
    # Box computed to overapproximate the solution at time δt
    q = evaluate.(evaluate.(qTM, δtI), (normalized_box,))
    # Box computed from the exact solution must be within q
    bb = all(fexact(t0+δtI, q0 .+ q0ξB) .⊆ q)
    # Display details if bb is false
    bb || @show(t0, domt, remainder.(qTM), 
            δt, δtI, q0ξ, q0ξB, q,
            fexact(t0+δtI, q0 .+ q0ξB))
    return bb
end


@testset "Tests for `validated_integ`" begin
    @testset "falling_ball!" begin
        @taylorize function falling_ball!(dx, x, p, t)
            dx[1] = x[2]
            dx[2] = -one(x[1])
            nothing
        end
        exactsol(t, t0, x0) = (x0[1] + x0[2]*(t-t0) - 0.5*(t-t0)^2, x0[2] - (t-t0))

        # Initial conditions
        tini, tend = 0.0, 10.0
        normalized_box = symmetric_box(2, Float64)
        q0 = [10.0, 0.0]
        δq0 = 0.25 * normalized_box
        X0 = IntervalBox(q0 .+ δq0)

        # Parameters
        abstol = 1e-20
        orderQ = 2
        orderT = 4
        ξ = set_variables("ξₓ ξᵥ", order=2*orderQ, numvars=length(q0))

        @testset "Forward integration 1" begin
            sol = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol)
            tTM = expansion_point(sol)
            qv  = flowpipe(sol)
            qTM = get_xTM(sol)
            @test length(qv) == length(qTM[1,:]) == length(sol)
            @test firstindex(sol) == 1
            @test sol[2] == get_xTM(sol,2)
            @test domain(sol,1) == 0..0

            end_idx = lastindex(sol)
            Random.seed!(1)
            for it = 1:_num_tests
                n = rand(2:end_idx)
                @test test_integ((t,x)->exactsol(t,tini,x), tTM[n], sol[n], q0, δq0)
            end

            solf = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol,
                adaptive=false)
            qvf, qTMf = getfield.((solf,), 2:3)

            @test length(qvf) == length(qv)
            @test qTM == qTMf

            # initializaton with a Taylor model
            X0tm = sol[1]
            sol2 = validated_integ(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol)
            qTM2 = getfield(sol2, 3)
            @test qTM == qTM2
            @test sol2[1,2] == X0tm[2]

            # Tests for TaylorModels.mince_in_time
            domT = TaylorModels.mince_in_time(sol)
            @test domT == expansion_point(sol) .+ domain(sol)
            timesdiv = TaylorModels.mince_in_time(sol, var=0, timediv=2)
            fpdiv = TaylorModels.mince_in_time(sol, var=1, timediv=2)
            @test timesdiv[3] ⊂ domT[2]
            @test hull(timesdiv[1],timesdiv[2]) == domT[1]
            @test fpdiv[3] ⊂ qv[2][1]
            @test hull(fpdiv[3],fpdiv[4]) ⊂ qv[2][1]
        end

        @testset "Forward integration 2" begin
            sol = validated_integ2(falling_ball!, X0, tini, tend, orderQ, orderT, abstol)
            tTM = expansion_point(sol)
            qv  = flowpipe(sol)
            qTM = get_xTM(sol)
            @test length(qv) == length(qTM[1,:]) == length(sol)
            @test firstindex(sol) == 1
            @test sol[2] == get_xTM(sol,2)
            @test domain(sol,1) == 0..0

            Random.seed!(1)
            end_idx = lastindex(sol)
            for it = 1:_num_tests
                n = rand(2:end_idx)
                @test test_integ((t,x)->exactsol(t,tini,x), tTM[n], sol[n], q0, δq0)
            end

            # initializaton with a Taylor model
            X0tm = get_xTM(sol,1)
            sol2 = validated_integ2(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol)
            qTM2 = getfield(sol2, 3)
            @test qTM == qTM2
            @test sol2[1,2] == X0tm[2]
        end

        # Initial conditions
        tini, tend = 10.0, 0.0
        q0 = [10.0, 0.0]
        δq0 = IntervalBox(-0.25 .. 0.25, 2)
        X0 = IntervalBox(q0 .+ δq0)

        @testset "Backward integration 1" begin
            sol = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol)
            tTM, qv, qTM = getfield.((sol,), 1:3)
            @test length(qv) == length(qTM[1,:]) == length(sol)
            @test firstindex(sol) == 1
            @test sol[2] == get_xTM(sol,2)
            @test domain(sol,1) == 0..0

            Random.seed!(1)
            end_idx = lastindex(sol)
            for it = 1:_num_tests
                n = rand(2:end_idx)
                @test test_integ((t,x)->exactsol(t,tini,x), expansion_point(sol,n), sol[n], q0, δq0)
            end

            solf = validated_integ(falling_ball!, X0, tini, tend, orderQ, orderT, abstol,
                adaptive=false)
            tTMf, qvf, qTMf = getfield.((solf,), 1:3)

            @test length(qvf) == length(qv)
            @test all(qTM .== qTMf)

            # initializaton with a Taylor model
            X0tm = sol[1]
            sol2 = validated_integ(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol)
            qv2, qTM2 = getfield.((sol2,), 2:3)
            @test qTM == qTM2

            sol2f = validated_integ(falling_ball!, X0tm, tini, tend, orderQ, orderT, abstol,
                adaptive=false)
            qv2f, qTM2f = getfield.((sol2f,), 2:3)
            @test length(qv2f) == length(qv2)
            @test all(qTM .== qTM2f)
        end

        @testset "Backward integration 2" begin
            sol = validated_integ2(falling_ball!, X0, tini, tend, orderQ, orderT, abstol)
            tTM, qv, qTM = getfield.((sol,), 1:3)
            @test length(qv) == length(qTM[1,:]) == length(sol)
            @test firstindex(sol) == 1
            @test sol[2] == get_xTM(sol,2)
            @test domain(sol,1) == 0..0

            Random.seed!(1)
            end_idx = lastindex(sol)
            for it = 1:_num_tests
                n = rand(2:end_idx)
                @test test_integ((t,x)->exactsol(t,tini,x), expansion_point(sol,n), sol[n], q0, δq0)
            end
        end
    end

    @testset "x_square!" begin
        @taylorize function x_square!(dx, x, p, t)
            dx[1] = x[1]^2
            nothing
        end

        exactsol(t, x0) = 1 / (1/x0[1] - t)

        tini, tend = 0., 0.45
        normalized_box = symmetric_box(1, Float64)
        abstol = 1e-15
        orderQ = 5
        orderT = 20
        q0 = [2.]
        δq0 = 0.0625 * normalized_box
        X0 = IntervalBox(q0 .+ δq0)
        ξ = set_variables("ξₓ", numvars=1, order=2*orderQ)

        @testset "Forward integration 1" begin
            sol = validated_integ(x_square!, X0, tini, tend, orderQ, orderT, abstol)
            tTM, qv, qTM = getfield.((sol,), 1:3)
            @test length(qv) == length(qTM[1, :]) == length(tTM)
            @test domain(sol,1) == 0..0

            Random.seed!(1)
            end_idx = lastindex(tTM)
            for it = 1:_num_tests
                n = rand(1:end_idx)
                @test test_integ((t,x)->exactsol(t,x), tTM[n], sol[n], q0, δq0)
            end

            solf = validated_integ(x_square!, X0, tini, tend, orderQ, orderT, abstol,
                adaptive=false)
            tTMf, qvf, qTMf = getfield.((solf,), 1:3)

            @test length(qvf) == length(qv)
            @test all(qTMf .== qTM)

            # initializaton with a Taylor model
            X0tm = copy(qTM[:, 1])
            sol2 = validated_integ(x_square!, X0tm, tini, tend, orderQ, orderT, abstol)
            tTM2, qv2, qTM2 = getfield.((sol2,), 1:3)
            @test qTM == qTM2
        end

        @testset "Forward integration 2" begin
            sol = validated_integ2(x_square!, X0, tini, tend, orderQ, orderT, abstol)
            tTM, qv, qTM = getfield.((sol,), 1:3)
            @test domain(sol,1) == 0..0

            @test length(qv) == length(qTM[1, :]) == length(tTM)

            Random.seed!(1)
            end_idx = lastindex(tTM)
            for it = 1:_num_tests
                n = rand(1:end_idx)
                @test test_integ((t,x)->exactsol(t,x), tTM[n], sol[n], q0, δq0)
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

        sol = validated_integ(pendulum!, X0, tini, tend, orderQ, orderT, abstol);
        @test all(ene0 .⊆ ene_pendulum.(flowpipe(sol)))

        sol = validated_integ2(pendulum!, X0, tini, tend, orderQ, orderT, abstol,
            validatesteps=32);
        @test all(ene0 .⊆ ene_pendulum.(flowpipe(sol)))

        # Initial conditions 2
        q0 = [1.1, 0.1, 0.0]
        δq0 = IntervalBox(-0.1 .. 0.1, -0.1 .. 0.1, -0.01 .. 0.01)
        X0 = IntervalBox(q0 .+ δq0)
        ene0 = ene_pendulum(X0)

        sol = validated_integ(pendulum!, X0, tini, tend, orderQ, orderT, abstol);
        @test all(ene0 .⊆ ene_pendulum.(flowpipe(sol)))

        sol = validated_integ2(pendulum!, X0, tini, tend, orderQ, orderT, abstol,
            validatesteps=32);
        @test all(ene0 .⊆ ene_pendulum.(flowpipe(sol)))
    end
end
