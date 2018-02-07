using Base.Test, TaylorModels, IntervalArithmetic

# import Documenter
# Documenter.makedocs(
#     modules = [TaylorModels],
#     format = :html,
#     sitename = "TaylorModels.jl",
#     root = joinpath(dirname(dirname(@__FILE__)), "docs"),
#     pages = Any["Home" => "index.md"],
#     strict = true,
#     linkcheck = true,
#     checkdocs = :exports,
#     authors = "David Sanders"
# )

# write your own tests here

@testset "Constructor" begin
    t = taylor1_var(3, 0, -0.5..0.5)

    @test t isa Taylor1Model
    @test t.n == 3
    @test t.x0 == 0
    @test t.I == -0.5..0.5
end

@testset "Exponential" begin
    t = taylor1_var(3, 0, -0.5..0.5)
    s = exp(t)

    @test s isa Taylor1Model

end
