
module ChordalDecomp

using SparseArrays, SuiteSparse, LinearAlgebra
using DataStructures


include("chordal_extension.jl")
include("clique_graph.jl")
include("chordal_decomp.jl")
include("utils.jl")

export preprocess!, sparsity_pattern, get_chordal_extension, get_cliques
export generate_clique_graph, merge_cliques!, get_cliques
export make_selectors
export build_perm_matrix

end
