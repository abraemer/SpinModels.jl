push!(LOAD_PATH,"../src/")

using Documenter
using SpinModels

makedocs(
    sitename = "SpinModels",
    format = Documenter.HTML(),
    modules = [SpinModels],
    authors = "Adrian Braemer",
    pages = [
        "Home" => "index.md",
        "couplings.md",
        "quick-ref.md",
        "Full Documentation" => Any[
            "full-docs/hamiltonian.md",
            "full-docs/geometry.md",
            "full-docs/interaction.md",
            "full-docs/internal.md"
        ],
    ],
)

deploydocs(
    repo = "github.com/abraemer/SpinModels.jl.git",
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
