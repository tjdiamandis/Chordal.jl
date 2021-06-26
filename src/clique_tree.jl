struct CliqueTree
    par::Vector{Int}
    child::Vector{Set{Int}}
    cliques::Vector{Vector{Int}}
    postordering::Vector{Int}
end


# Input:
#   L       : sparse lower triangular matrix that rep chordal graph
#   cliques : vector of cliques
function generate_clique_tree(L, cliques)
    n = size(L, 1)
    c_int_graph = generate_clique_graph(cliques, n, int_graph=true)
    edges_tree = max_weight_span_tree(c_int_graph, cliques)
    root, par, child = det_par_child(L, cliques, edges_tree)
    postordering = get_postordering(par, child)

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


function get_postordering(par, child)
    n = length(par)
    root = findfirst(x->par[x] == 0, 1:n)
    isnothing(root) && error(ArgumentError("Can't find root; set par[root] := 0"))
    visit_order = zeros(Int, n)

    # post order via DFS (https://en.wikipedia.org/wiki/Depth-first_search)
    curr = n
    S = [root]
    while !isempty(S)
        v = pop!(S)
        visit_order[v] = curr
        curr -= 1
        push!(S, child[v]...)
    end
    postordering = sort(1:n, by=x->visit_order[x])

    return postordering
end


"""
    etree(A)
Compute the elimination tree of a lower triangular sparse matrix `L`
"""
function etree(L::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    n = size(L, 1)
    par = zeros(Ti, n)
    ancestor = zeros(Ti, n)

    _etree(sparse(L'), par, ancestor)
    count(x->x==0, par) > 1 && error(ArgumentError("L must be connected"))
    return par
end

# (almost) nonallocating verion of above
# Computes using the upper triangular part
# TODO: change to lower triangular??
# Davis. Direct Methods for Sparse Linear Systems, pg 42
function _etree(U::SparseMatrixCSC{Tv, Ti}, par::Vector{Ti}, ancestor::Vector{Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    n = size(U, 1)
    for col::Ti in 1:n, p in nzrange(U, col) #A.colptr[col]:(A.colptr[col+1]-1)
        i = rowvals(U)[p]
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


# Algorithm from Pothen and Sun
#   restated in VA, Algorithm 4.1
function max_supernode_etree(L::SparseMatrixCSC, etree_par::Vector{Int})
    n = size(A, 1)
    etree_children = get_children_from_par(etree_par)
    post_ord = get_postordering(etree_par, etree_children)
    deg⁺ = get_higher_deg(L)

    # Vc = representative vertices
    Vc = Vector{Int}(undef, 0)
    sizehint!(Vc, n)
    snd_membership = zeros(Int, n)
    snd_par = zeros(Int, n)

    for v = post_ord
        ŵ = -1
        for w in etree_children[v]
            deg⁺[v] > deg⁺[w] - 1 && continue
            ŵ = w
            break
        end

        # v is a rep vertex
        if ŵ == -1
            push!(Vc, v)
            snd_membership[v] = v
            u = v
        # otherwise, add v to snd(ŵ)
        else
            u = snd_membership[ŵ]
            snd_membership[v] = u
        end

        for w in etree_children[v]
            w == ŵ && continue
            # q(z) = u, a(z) = v; z is the vertex in Vc that satisfies w ∈ snd(z)
            z = snd_membership[w]
            snd_par[z] = u
            # a[z] = v
        end
    end

    return Vc, snd_par, snd_membership

end


# input: parent vector
# returns vector v st v[i] = {k | k is a child of i}
function get_children_from_par(par::Vector{Int})
    n = length(par)
    children = [Set{Int}() for _ in 1:n]

    for i in 1:n
        par_i = par[i]
        if par_i != 0
            push!(children[par_i], i)
        end
    end

    return children
end


# NOTE: Assumes that L has order σ = 1:n & is lower triangular
function get_higher_deg(L::SparseMatrixCSC)
    n = size(L, 1)
    deg⁺ = zeros(Int, n)

    for i in 1:n
        deg⁺[i] = L.colptr[i+1] - L.colptr[i]
    end

    return deg⁺
end
