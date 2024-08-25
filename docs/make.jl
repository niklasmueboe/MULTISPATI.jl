using Documenter
using DocumenterInterLinks
using Multispati

links = InterLinks(
# "CategoricalArrays" => "https://categoricalarrays.juliadata.org/stable/",
)

makedocs(;
    sitename="Multispati.jl",
    pages=["Home" => "index.md", "Reference API" => "api.md"],
    authors="Niklas Müller-Bötticher",
    plugins=[links],
)

deploydocs(; repo="github.com/niklasmueboe/Multispati.jl.git")
