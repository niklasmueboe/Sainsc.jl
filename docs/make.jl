using Documenter
using DocumenterInterLinks
using Sainsc

links = InterLinks(
    # "CategoricalArrays" => "https://categoricalarrays.juliadata.org/stable/",
    # "DataFrames" => "https://dataframes.juliadata.org/stable/",
    # "OffsetArrays" => "https://juliaarrays.github.io/OffsetArrays.jl/stable/",
    "Muon" => "https://scverse.org/Muon.jl/dev/",
    "SparseArrays" => "https://docs.julialang.org/en/v1/",
)

makedocs(;
    sitename="Sainsc.jl",
    pages=[
        "Home" => "index.md",
        "Example analysis" => "examples/ExampleAnalysis.md",
        "Reference API" => "api.md",
    ],
    authors="Niklas Müller-Bötticher",
    plugins=[links],
    format=Documenter.HTML(; size_threshold_ignore=["examples/ExampleAnalysis.md"]),
)

deploydocs(; repo="github.com/niklasmueboe/Sainsc.jl.git")
