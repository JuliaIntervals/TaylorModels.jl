import Documenter

Documenter.deploydocs(
    repo = "github.com/dpsanders/TaylorModels.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
