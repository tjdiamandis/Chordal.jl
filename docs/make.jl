using ChordalDecomp
using Documenter

makedocs(;
    modules=[ChordalDecomp],
    authors="Theo Diamandis",
    repo="https://github.com/tjdiamandis/Chordal.jl/blob/{commit}{path}#L{line}",
    sitename="Chordal.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tjdiamandis.github.io/Chordal.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tjdiamandis/Chordal.jl",
)
