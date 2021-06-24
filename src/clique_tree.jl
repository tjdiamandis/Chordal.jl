struct CliqueTree
    par::Vector{Int}
    child::Vector{Set{Int}}
    cliques::Vector{Vector{Int}}
    postordering::Vector{Int}
end


"""
    etree(A)
Compute the elimination tree of a symmetric sparse matrix `A` from `triu(A)`.
"""
function etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && throw(ArgumentError("Matrix must be symmetric"))
    m, n = size(A)
    par = zeros(Ti, n)
    ancestor = zeros(Ti, n)

    _etree(A, par, ancestor)
    return par
end

# (almost) nonallocating verion of above
function _etree(A::SparseMatrixCSC{Tv, Ti}, par::Vector{Ti}, ancestor::Vector{Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    for col::Ti in 1:n, p in nzrange(A, col) #A.colptr[col]:(A.colptr[col+1]-1)
        i = rowvals(A)[p]
        while !iszero(i) && i < col
            i_next = ancestor[i]
            ancestor[i] = col
            if iszero(i_next)
                par[i] = col
            end
            i = i_next
        end
    end
    return nothing
end


# Input:
#   L       : sparse lower triangular matrix that rep chordal graph
#   cliques : vector of cliques
function generate_clique_tree(L, cliques)
    n = size(L, 1)
    c_int_graph = generate_clique_graph(cliques, n, int_graph=true)
    edges_tree = max_weight_span_tree(c_int_graph, cliques)
    root, par, child = det_par_child(L, cliques, edges_tree)
    postordering = get_postordering(root, par, child)

    return CliqueTree(par, child, cliques, postordering)
end


function max_weight_span_tree(cg::CliqueGraph, cliques)
    m = length(cg)
    edges_tree = spzeros(Int, m, m)

    cliques_added = [argmax([length(c) for c in cliques])]
    cliques_to_add = collect(1:m)
    deleteat!(cliques_to_add, cliques_added[1])

    @views for r in 1:m-1
        max_int_ind = argmax(cg.edge_mat[cliques_added, cliques_to_add])
        C = cliques_added[max_int_ind[1]]
        C_new = cliques_to_add[max_int_ind[2]]

        edges_tree[C, C_new] = edges_tree[C_new, C] = 1
        append!(cliques_added, C_new)
        deleteat!(cliques_to_add, max_int_ind[2])
    end

    return edges_tree
end


# Determine the parent-child structure of the tree
function det_par_child(L, cliques, edges_tree)
    root_vert = argmax(vec(sum(L, dims=1)))
    root = findfirst(x -> root_vert in x, cliques)

    par = zeros(Int, size(edges_tree, 1))
    child = [Set{Int}() for _ in 1:length(par)]
    par[root] = 0

    det_children!(par, child, root, edges_tree)

    return root, par, child
end


# Traverses tree to fill parent and child data structures
function det_children!(par, child, node, edges_tree)
    children = filter(x->(x != par[node]), findnz(edges_tree[:,node])[1])

    if isempty(children)
        return nothing
    end

    child[node] = Set(children)
    for c in children
        par[c] = node
        det_children!(par, child, c, edges_tree)
    end
    return nothing
end


function get_postordering(root, par, child)
    m = length(par)
    @assert par[root] == 0
    visit_order = zeros(Int, m)

    # post order via DFS (https://en.wikipedia.org/wiki/Depth-first_search)
    curr = m
    S = [root]
    while !isempty(S)
        v = pop!(S)
        visit_order[v] = curr
        curr -= 1
        push!(S, child[v]...)
    end
    postordering = sort(1:m, by=x->visit_order[x])

    return postordering
end
