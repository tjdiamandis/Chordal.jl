# These algorithms are from
#   GARSTKA M., CANNON M., GOULART P.:
#   A clique graph based merging strategy for decomposable sdps.
#   IFAC-PapersOnLine 53, 2 (2020), 7355–7361.

import Base.length
length(cg::CliqueGraph) = size(cg.edge_mat, 1)


# Input:
#   cliques : vector of cliques
#   n       : num nodes of graph
# Output:
#   CliqueGraph : rep of relationship bt cliques
"""
    generate_clique_graph(cliques, n::Integer)

Generates CliqueGraph datastructure for an undirected graph from `cliques` and
the number of nodes `n`.

Reference
* [A clique graph based merging strategy for decomposable SDPs](https://arxiv.org/pdf/1911.05615)
by Michael Garstka, Mark Cannon, Paul Goulart
"""
function generate_clique_graph(cliques, n::Integer)
    m = length(cliques)
    membership_mat = spzeros(Bool, n, m)

    # TODO: edge mat as matrix or as list of tuples?
    edge_mat = zeros(Float64, m, m)
    active_cliques = Set(1:m)

    for (idx, clique) in enumerate(cliques)
        for node in clique
            membership_mat[node, idx] = true
        end
    end
    for i in 1:m, j in i+1:m
        ci = @view(membership_mat[:, i])
        cj = @view(membership_mat[:, j])
        l_int = sum(ci .& cj)
        if l_int > 0
            edge_mat[i,j] = edge_mat[j,i] = weight_function(ci, cj)
        end
    end

    return CliqueGraph(membership_mat, edge_mat, active_cliques)
end


# TODO: this should be a user-defined function
"""
    weight_function(ci, cj)

Defines the weight function used to determine if cliques `ci` and `cj` should be
merged. Cliques are merged if this is positive.

Reference
* [A clique graph based merging strategy for decomposable SDPs](https://arxiv.org/pdf/1911.05615)
by Michael Garstka, Mark Cannon, Paul Goulart
"""
function weight_function(ci, cj)
    l_ci = sum(ci)
    l_cj = sum(cj)
    l_union = sum(ci[k] | cj[k] for k in 1:length(ci))
    # l_union = 0.0
    # for k in 1:length(ci)
    #     # l_ci += ci[k]
    #     # l_cj += cj[k]
    #     l_union += ci[k] .| cj[k]
    # end
    return l_ci^3 + l_cj^3 - l_union^3
end


function _merge_cliques!(cg::CliqueGraph, i::Integer, j::Integer)
    # nodes X cliques
    n, m = size(cg.membership_mat)

    ci = @view(cg.membership_mat[:, i])
    cj = @view(cg.membership_mat[:, j])

    n_rem = sum(ci) < sum(cj) ? i : j
    n_keep = sum(ci) < sum(cj) ? j : i
    c_keep = @view(cg.membership_mat[:,n_keep])
    c_rem =  @view(cg.membership_mat[:,n_rem])

    # Update membership matrix
    c_keep .= ci .| cj
    c_rem .= false
    cg.edge_mat[n_keep,n_rem] = cg.edge_mat[n_rem,n_keep] = 0.0

    # TODO: change to only cliques with nonzero intersection
    for k in 1:m
        k == n_keep && continue
        # Remove edges from n_rem
        cg.edge_mat[k,n_rem] = cg.edge_mat[n_rem,k] = 0.0

        # Update edges to n_keep
        wij = weight_function(c_keep, @view(cg.membership_mat[:,k]))
        cg.edge_mat[k,n_keep] = cg.edge_mat[n_keep,k] = wij
    end

    pop!(cg.active_cliques, n_rem)

end


"""
    merge_cliques!(cg::CliqueGraph; verbose=false)

Merges cliques in `cg` in a greedy fashion starting with the pair with the largest
`weight_function`. Stops when all pairs of cliques have a non-positive weight.

Reference
* [A clique graph based merging strategy for decomposable SDPs](https://arxiv.org/pdf/1911.05615)
by Michael Garstka, Mark Cannon, Paul Goulart
"""
function merge_cliques!(cg::CliqueGraph; verbose=false)
    # nodes X cliques
    n, m = size(cg.membership_mat)
    if verbose
        max_clique_size = maximum(nnz, @view(cg.membership_mat[:,i]) for i in 1:m)
        @info "Starting merging; num cliques = $m. Max size is $max_clique_size."
    end

    max_val, max_ind = findmax(cg.edge_mat)
    while max_val > 0
        #TODO: verbose printing
        i = max_ind[1]
        j = max_ind[2]
        _merge_cliques!(cg, i, j)
        max_val, max_ind = findmax(cg.edge_mat)
    end
    dropzeros!(cg.membership_mat)
    if verbose
        max_clique_size = maximum(nnz, @view(cg.membership_mat[:,i]) for i in cg.active_cliques)
        @info "Finished merging; num cliques = $(length(cg.active_cliques)). Max size is $max_clique_size."
    end
end


"""
    get_cliques(cg::CliqueGraph)

Returns a list of the cliques in CliqueGraph `cg`.
"""
function get_cliques(cg::CliqueGraph)
    return [findnz(cg.membership_mat[:,i])[1] for i in sort(collect(cg.active_cliques))]
end
