using Test
using TaylorModels
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(TaylorModels)
    ua = Aqua.detect_unbound_args_recursively(TaylorModels)
    @test length(ua) == 0

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(TaylorModels; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("TaylorModels", pkgdir(last(x).module)), ambs)
    for method_ambiguity in ambs
        @show method_ambiguity
    end
    @test length(ambs) == 0
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(TaylorModels)
    Aqua.test_deps_compat(TaylorModels)
    Aqua.test_stale_deps(TaylorModels)
    Aqua.test_piracies(TaylorModels)
    Aqua.test_unbound_args(TaylorModels)
    Aqua.test_project_extras(TaylorModels)
    Aqua.test_persistent_tasks(TaylorModels)
end

