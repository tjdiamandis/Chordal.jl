struct CliqueGraph
    membership_mat::SparseMatrixCSC{Bool, Int}  # nodes x num_cliques
    edge_mat::Matrix{Float64}                   # num_cliques x num_cliques
    # TODO: change to bitset?
    active_cliques::Set{Int}
end
