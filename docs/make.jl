using WaveMC
using Documenter: Documenter, DocMeta, deploydocs, makedocs

DocMeta.setdocmeta!(
    WaveMC, :DocTestSetup, :(using WaveMC); recursive = true
)

# Copy the license file into docs/src
cp(joinpath(@__DIR__, "..", "LICENSE"), joinpath(@__DIR__, "src", "LICENSE"); force=true)


include("make_index.jl")

makedocs(;
    modules = [WaveMC],
    authors = "Andreas Feuerpfeil <development@manybodylab.com>",
    sitename = "WaveMC.jl",
    format = Documenter.HTML(;
        canonical = "https://manybodylab.github.io/WaveMC.jl",
        edit_link = "main",
        assets = [#"assets/logo.png", 
            "assets/extras.css"],
    ),
    pages = ["Home" => "index.md", "Reference" => "reference.md"],
)

deploydocs(;
    repo = "github.com/ManyBodyLab/WaveMC.jl", devbranch = "main", push_preview = true
)
