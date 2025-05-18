# Tests related to validated integration

using TaylorModels, TaylorModels.ValidatedInteg
using LinearAlgebra: norm
using Test

# const _num_tests = 1000
# const α_mid = TaylorModels.α_mid

setdisplay(:full)

# This is based on an example from Buenger, Numer Algor (2018) 78:1001–1017,
# which shows the efectivity of shrink-wrapping
@testset "Testing `shrink_wrapping` 1" begin
    _order = 2
    set_variables("x", numvars=1, order=2*_order)
    x0 = SVector{1}(0..0)
    dom = symmetric_box(Float64, 1)#IntervalBox(-1..1, 1)
    for δ in 1/16:1/16:1
        Δ = interval(-δ, δ)
        p = TaylorModelN(TaylorN(1, order=_order), Δ, x0, dom)
        # Usual TM case
        result = p * p
        rem_result = δ*(2+δ)*dom[1]
        range_result = dom[1]^2 + rem_result
        @test isequal_interval(result(dom), range_result)
        @test isequal_interval(remainder(result), rem_result)
        # Shrink-wrapping:
        q = [ p ]
        shrink_wrapping!(q)
        @test isequal_interval(p(dom), q[1](dom))
        @test issubset_interval(remainder(q[1]), Δ)
        @test isequal_interval(remainder(q[1]), x0[1])
        # Result using shrink-wrapped variables
        result_sw = q[1] * q[1]
        @test issubset_interval(remainder(result_sw), rem_result)
        @test issubset_interval(result_sw(dom), range_result)
    end
end

# This is based on the example of Sect 4 of M. Berz, K. Makino,
# in https://bt.pa.msu.edu/pub/papers/VIShrink06/VIShrink06.pdf
@testset "Testing `shrink_wrapping` 2" begin
    function two_state1!(v)
        x, y = v
        x_2 = x^2
        y_2 = y^2
        r_2 = x_2 + y_2
        s = sqrt(1 + r_2)
        v[1] = x * s
        v[2] = y * s
        return v
    end

    function two_state2!(v)
        x1, y1 = v
        x1_2 = x1^2
        y1_2 = y1^2
        r1_2 = x1_2 + y1_2
        r1 = sqrt( 1 + 4*r1_2 )
        s1 = 2 / (1+r1)
        ss = sqrt( s1 )
        v[1] = x1 * ss
        v[2] = y1 * ss
        return v
    end

    local B = symmetric_box(Float64, 2)
    local δ = 0.05
    local ib = δ * B
    local mib = fill( 0..0, SVector{2})
    local zi = interval(zero(Float64))

    # Diverges using naive Interval arithmetic methods
    ib0 = [ib...,]
    for iter = 1:210 # iter is twice the number of iterates
        two_state1!(ib0)
        two_state2!(ib0)
    end
    @test minimum(diam.(ib0)) > 1e6

    # Taylor model controls grow and behaves *essentially* as the identity
    # (due to the normalization to the symmetric box) after two iterates
    order = 10
    set_variables("x y", order=2*order)
    vm = [TaylorModelN(normalize_taylor(TaylorN(i, order=order), ib, true),
        zi, mib, B) for i = 1:2]
    vm0 = deepcopy(vm)
    for iter = 1:1_000
        two_state1!(vm0)
        two_state2!(vm0)
    end
    @test maximum(mag.(remainder.(vm0))) < 1e-6 # 2.08e-7
    # Constant and linear terms are identic to the initial ones
    @test all(isequal_interval.(constant_term.(vm0), constant_term.(vm)))
    @test all(linear_polynomial.(vm0) .== linear_polynomial.(vm))
    # Maximum difference of polynomials is very small
    @test mag(norm(polynomial(vm0[1] - vm[1]), Inf)) < 1e-14

    # Taylor model with shrink-wrapping after two iterates
    vm0 .= deepcopy.(vm)
    for iter = 1:10#00
        two_state1!(vm0)
        two_state2!(vm0)
        shrink_wrapping!(vm0)
    end
    @test maximum(mag.(remainder.(vm0))) < 1e-8 #2.2e-13
    # The dominating difference is in the linear term
    @test isequal_interval(norm(polynomial(vm0[1] - vm[1]), Inf),
            norm(vm0[1][1] - vm[1][1], Inf))


    # Taylor model with shrink-wrapping after each iterate
    vm0 .= deepcopy.(vm)
    for iter = 1:1#00
        two_state1!(vm0)
        shrink_wrapping!(vm0)
        two_state2!(vm0)
        shrink_wrapping!(vm0)
    end
    @test maximum(mag.(remainder.(vm0))) < 2.2e-13
    # The dominating difference is in the linear term
    @test isequal_interval(norm(polynomial(vm0[1] - vm[1]), Inf),
            norm(vm0[1][1] - vm[1][1], Inf))

    # Test AssertionError
    vm = [TaylorModelN(TaylorN(i, order=order), zi, mib, ib) for i = 1:2]
    @test_throws AssertionError shrink_wrapping!(vm)
end
