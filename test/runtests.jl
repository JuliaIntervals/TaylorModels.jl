using TaylorModels

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

# using Base.Test

include("TM1.jl")
include("TMN.jl")
