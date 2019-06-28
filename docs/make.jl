using Documenter, TaylorModels

makedocs(
    modules = [TaylorModels],
    format = Documenter.HTML(),
    sitename = "TaylorModels.jl",
    authors  = "Luis Benet and David P. Sanders",
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ]
)

Documenter.deploydocs(
    repo = "github.com/JuliaIntervals/TaylorModels.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)
