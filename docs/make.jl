using Chordal
using Documenter
using Literate

# For Plots.jl
# https://discourse.julialang.org/t/plotting-errors-when-building-documentation-using-plots-jl-and-documenter-jl/67849
ENV["GKSwstype"]="100"

EXCLUDED_EXAMPLES = ["chordal_graph.jl", "sdp_standard.jl"]

# utility function from https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")

# utility functions from https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl
fix_suffix(filename) = replace(filename, ".jl" => ".md")
function postprocess(content)
      """
      The source files for all examples can be found in [/examples](https://github.com/tjdiamandis/Chordal.jl/tree/main/examples).
      """ * content
end

examples_path = joinpath(@__DIR__, "../examples/")
examples = filter(x -> endswith(x, ".jl") && !in(x, EXCLUDED_EXAMPLES), readdir(examples_path))
build_path =  joinpath(@__DIR__, "src", "examples/")

for example in examples
      Literate.markdown(
        examples_path * example, build_path;
        preprocess = fix_math_md,
        postprocess = postprocess,
        flavor = Literate.DocumenterFlavor(),
        credit = true
    )
end

examples_nav = fix_suffix.(joinpath.("examples", examples))

makedocs(;
    modules=[Chordal],
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
        "Examples" => examples_nav,
        "API Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/tjdiamandis/Chordal.jl",
    devbranch = "main"
)
