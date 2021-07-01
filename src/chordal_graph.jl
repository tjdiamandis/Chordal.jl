using LightGraphs: AbstractGraph, adjacency_matrix

struct Graph
    neighbors::Vector{Set{Int}}
    active::BitSet
    peo::Vector{Int}
end



# # SIMPLE LINEAR-TIME ALGORITHMS TO TEST CHORDALITY OF GRAPHS, TEST ACYCLICITY OF HYPERGRAPHS, AND SELECTIVELY REDUCE ACYCLIC HYPERGRAPHS
# # ROBERT E. TARJAN AND MIHALIS YANNAKAKIS
# Algorithm 4.2 in VA
"""
    maximum_cardinality_search(A)
Compute the perfect elimination ordering for a chordal graph represented by
a symmetric sparse matrix `A`. This function can also be used with AbstractGraph
objects from `LightGraphs.jl`.
"""
function maximum_cardinality_search(A::SparseMatrixCSC)
    !issymmetric(A) && error(ArgumentError("Matrix must be symmetric"))
    n = size(A, 1)
    peo = zeros(Int, n)
    cache = zeros(Int, n)

    _maximum_cardinality_search!(peo, A; n = n, cache=cache)
    return peo
end
maximum_cardinality_search(G::AbstractGraph) = maximum_cardinality_search(adjacency_matrix(G))

function _maximum_cardinality_search!(peo, A; n=n, cache=cache)
    v = 1
    for i in n:-1:1
        max_val = 0
        for j in 1:n
            peo[j] != 0 && continue
            if cache[j] > max_val
                max_val = cache[j]
                v = j
            end
        end
        peo[v] = i
        @views cache[(rowvals(A)[nzrange(A, v)])] .+= 1
    end

    return peo
end


"""
    is_chordal(A)

Tests if the graph represented by symmetric sparse matrix `A` is chordal. This
function can also be used with AbstractGraph objects from `LightGraphs.jl`.

References
* [Simple linear-time algorithms to test chordality of graphs, test acyclicity
of hypergraphs, and selectively reduce acyclic hypergraphs](https://epubs.siam.org/doi/pdf/10.1137/0213035?casa_token=A22jkwrsrL0AAAAA:rx-G6F21ubTkMiJmRTH3IKqxmFTo_IVWDDZfJig5lsZxnQtNH2vUKWfZ3eZJKUv9CiKbIPt1VQs)
by Robert Tarjan and Mihalis Yannakakis
"""
function is_chordal(A::SparseMatrixCSC; peo=nothing)
    if isnothing(peo)
        peo = maximum_cardinality_search(A)
    end
    i_peo = invperm(peo)

    n = length(peo)
    f = zeros(Int, n)
    ind = zeros(Int, n)

    for i in 1:n
        w = i_peo[i]
        f[w] = w
        ind[w] = i
        for v in rowvals(A)[nzrange(A, w)]
            peo[v] >= i && continue
            ind[v] = i
            if f[v] == v
                f[v] = w
            end
        end

        for v in rowvals(A)[nzrange(A, w)]
            peo[v] >= i && continue
            ind[f[v]] < i && return false
        end
    end
    return true
end

is_chordal(G::AbstractGraph; peo=nothing) = is_chordal(adjacency_matrix(G); peo=peo)
