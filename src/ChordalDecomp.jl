
module ChordalDecomp

using SparseArrays, SuiteSparse, LinearAlgebra
using DataStructures


include("chordal_extension.jl")
include("clique_graph.jl")
include("chordal_decomp.jl")
include("utils.jl")

export preprocess!, sparsity_pattern, get_chordal_extension, get_cliques
export generate_clique_graph, merge_cliques!, get_cliques
export get_selectors, make_selectors_from_cliques
export build_perm_matrix, generate_random_sdp

end
