using Documenter
using DocumenterInterLinks
using MULTISPATI

links = InterLinks(
# "CategoricalArrays" => "https://categoricalarrays.juliadata.org/stable/",
)

makedocs(;
    sitename="MULTISPATI.jl",
    pages=["Home" => "index.md", "Reference API" => "api.md"],
    authors="Niklas Müller-Bötticher",
    plugins=[links],
)

deploydocs(; repo="github.com/niklasmueboe/MULTISPATI.jl.git")
