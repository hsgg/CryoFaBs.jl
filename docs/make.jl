using CryoFaBs
using Documenter

DocMeta.setdocmeta!(CryoFaBs, :DocTestSetup, :(using CryoFaBs); recursive=true)

makedocs(;
    modules=[CryoFaBs],
    authors="Henry Gebhardt <henry.s.gebhardt@jpl.nasa.gov> and contributors",
    repo="https://github.com/hsgg/CryoFaBs.jl/blob/{commit}{path}#{line}",
    sitename="CryoFaBs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hsgg.github.io/CryoFaBs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hsgg/CryoFaBs.jl",
)
