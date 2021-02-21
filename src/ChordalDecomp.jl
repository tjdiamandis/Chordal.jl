
module ChordalDecomp

using SparseArrays, SuiteSparse, LinearAlgebra
using DataStructures


include("chordal_extension.jl")
include("clique_graph.jl")
include("chordal_decomp.jl")
include("utils.jl")

export sparsity_pattern, get_chordal_extension
export generate_clique_graph, merge_cliques!, get_cliques
export get_selectors, make_selectors_from_cliques, make_selectors_from_clique_graph
export generate_random_sdp

end
