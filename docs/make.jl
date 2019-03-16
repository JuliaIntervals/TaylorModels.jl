using Documenter, TaylorModels

makedocs(
    modules = [TaylorModels],
    format = Documenter.HTML(),
    sitename = "TaylorModels",
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ]
)

Documenter.deploydocs(
    repo = "github.com/dpsanders/TaylorModels.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
