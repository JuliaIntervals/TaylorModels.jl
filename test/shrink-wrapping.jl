# Tests related to validated integration

using TaylorModels
# using LinearAlgebra: norm
using Test

# const _num_tests = 1000
# const α_mid = TaylorModels.α_mid

setformat(:full)

# This is based on the examples from Buenger, Numer Algor (2018) 78:1001–1017,
# showing efectivity of shrink-wrapping
@testset "Testing `shrink_wrapping` 1" begin
    _order = 2
    set_variables("x", numvars=1, order=2*_order)
    x0 = IntervalBox(0..0, 1)
    dom = IntervalBox(-1..1, 1)
    for δ in 1/16:1/16:1
        rem = δ * Interval(-1, 1)
        p1 = TaylorModelN(TaylorN(1, order=_order), rem, x0, dom)
        p2 = deepcopy(p1)
        res = p1 * p2
        bound_prod = dom[1]^2 + δ*(2+δ)*dom[1]
        @test res(dom) ⊆ polynomial(res)(dom) + bound_prod
        @test remainder(res) ⊆ bound_prod
        q = [ res ]
        TaylorModels.shrink_wrapping!(q)
        @test q[1](dom) ⊆ res.pol(dom) + bound_prod
        @test remainder(q[1]) ⊆ bound_prod
        @test remainder(q[1]) == 0..0
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

    local δ = 0.05
    local ib = IntervalBox( -δ .. δ, 2)
    local mib = IntervalBox( 0..0, 2)

    # Diverges using naive Interval arithmetic methods
    ib0 = [ib...]
    for iter = 1:210 # iter is twice the number of iterates
        two_state1!(ib0)
        two_state2!(ib0)
    end
    @test minimum(diam.(ib0)) > 1e6

    # Taylor model controls grow and behaves as identity for after two iterates
    order = 10
    set_variables("x y", order=2*order)
    vm = [TaylorModelN(TaylorN(i, order=order), 0..0, mib, ib) for i = 1:2]
    vm0 = [deepcopy(vm)...]
    for iter = 1:1_000
        two_state1!(vm0)
        two_state2!(vm0)
    end
    @test maximum(remainder.(vm0)) < 1e-6 # 2.2e-7
    @test all(polynomial.(vm0) .== polynomial.(vm)) # only remainder grows!

    # Taylor model with shrink-wrapping after two iterates
    vm0 .= [deepcopy(vm)...]
    for iter = 1:1_000
        two_state1!(vm0)
        two_state2!(vm0)
        TaylorModels.shrink_wrapping!(vm0)
    end
    @test maximum(remainder.(vm0)) < 1e-8 #2.2e-13

    # Taylor model with shrink-wrapping after each iterate
    vm0 .= [deepcopy(vm)...]
    for iter = 1:1_000
        two_state1!(vm0)
        TaylorModels.shrink_wrapping!(vm0)
        two_state2!(vm0)
        TaylorModels.shrink_wrapping!(vm0)
    end
    @test maximum(remainder.(vm0)) < 2.2e-13
end
