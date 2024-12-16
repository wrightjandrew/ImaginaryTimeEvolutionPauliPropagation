using Documenter, PauliPropagation

DocMeta.setdocmeta!(PauliPropagation, :DocTestSetup, :(using PauliPropagation); recursive=true)

makedocs(;
    modules=[PauliPropagation],
    doctest=true,
    linkcheck=true,
    authors="Manuel S. Rudolph <manuel.rudolph@epfl.ch> and contributors",
    repo="https://github.com/MSRudolph/PauliPropagation.jl/blob/{commit}{path}#{line}",
    sitename="PauliPropagation.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://MSRudolph.github.io/PauliPropagation.jl",
        assets=["assets/style.css"],
        size_threshold=400 * 2^10
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md",
    ],
    checkdocs=:none,  # :exports,
)

deploydocs(; repo="github.com/MSRudolph/PauliPropagation.jl", push_preview=false)

# TODO: Move all constructors into types, give them docstrings, and remove docstring from type itself.