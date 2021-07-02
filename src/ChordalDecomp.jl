
module ChordalDecomp

using SparseArrays, SuiteSparse, LinearAlgebra
using JuMP
using DataStructures


include("types.jl")
include("utils.jl")
include("chordal_decomp.jl")
include("chordal_extension.jl")
include("chordal_graph.jl")
include("clique_graph.jl")
include("clique_tree.jl")
include("completion.jl")
include("jump.jl")

export sparsity_pattern, get_chordal_extension
export generate_clique_graph, merge_cliques!, get_cliques
export get_selectors, make_selectors_from_cliques, make_selectors_from_clique_graph
export maxdet_completion
export build_constraints_lmi!
export get_sparsity_pattern_from_cliques, reconstruct_from_sparse_varref

end
