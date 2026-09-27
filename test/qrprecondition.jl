# Tests using qrprecondition

using TaylorModels
using TaylorModels.ValidatedInteg
using Test
# using Random
# using StaticArrays

const _num_tests = 1_000

# const TI = TaylorIntegration
# const VI = TaylorModels.ValidatedInteg

setdisplay(:full)

@testset "Tests for `qrprecondition`" begin
    sp = JetSpace(4, ["x", "y"])
    x, y = variables(sp)
    q0 = [interval(0), interval(0)]
    dom = [interval(-1,1), interval(-1,1)]
    xTMN = TaylorModelN(1+3x-2y+x^2-y^2, interval(-1,1) ./ 8, q0, dom)
    yTMN = TaylorModelN(1+x+2y+x^2+x*y, interval(-1,1) ./ 16, q0, dom)
    vTMN = [xTMN, yTMN]
    leftTMN, rightTMN = qrprecondition(vTMN)
    res = affine_compose(leftTMN, rightTMN)
    for ind in eachindex(vTMN)
        @test isapprox(polynomial(vTMN[ind]), polynomial(res[ind]))
        @test issubset_interval(remainder(vTMN[ind]), remainder(res[ind]))
    end
    for i in eachindex(vTMN)
        p, r = polynomial(vTMN[i]), polynomial(res[i])
        for ordQ in eachindex(p.coeffs)
            @test all(isapprox.(p.coeffs[ordQ].coeffs, r.coeffs[ordQ].coeffs; atol=1e-9))
        end
    end
    #
    # Quadratic example by Neher, Jackson, Nedialkov
    xTMN = TaylorModelN(0.904667 + 0.0505*x + 0.005*y, interval(−5.09307E-5, 7.86167E-5), q0, dom)
    yTMN = TaylorModelN(−0.909333 + 0.0095*x + 0.0505*y + 0.00025*x^2,
        interval(−1.75707E-4, 1.60933E-4), q0, dom)
    vTMN = [xTMN, yTMN]
    leftTMN, rightTMN = qrprecondition(vTMN)
    affine_compose!(res, leftTMN, rightTMN)
    for ind in eachindex(vTMN)
        @test isapprox(polynomial(vTMN[ind]), polynomial(res[ind]))
        @test issubset_interval(remainder(vTMN[ind]), remainder(res[ind]))
    end
end
