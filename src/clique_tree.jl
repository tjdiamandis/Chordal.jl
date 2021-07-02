struct EliminationTree
    par::Vector{Int}
    child::Vector{Vector{Int}}
    peo::Vector{Int}
    function EliminationTree(L::SparseMatrixCSC; peo::Vector{Int})
        n = size(L, 1)
        if !is_chordal(L; peo=1:n)
            error(ArgumentError("L must be chordal"))
        end
        par = get_etree(L)
        child = get_children_from_par(par)
        new(par, child, peo)
    end
end


struct CliqueTree
    etree::EliminationTree
    vreps::Vector{Int}
    snd_par::Vector{Int}
    snd_child::Vector{Vector{Int}}
    snds::Vector{Vector{Int}}
    seps::Vector{Vector{Int}}
    postordering::Vector{Int}
    perm::Vector{Int}
    function CliqueTree(L::SparseMatrixCSC; peo=nothing)
        n = size(L, 1)
        if isnothing(peo)
            peo = collect(1:n)
        end
        etree = EliminationTree(L; peo=peo)
        vreps, snd_par, snd_membership = max_supernode_etree(L, etree.par)
        snd_children = get_children_from_par(snd_par)
        postordering = get_postordering(snd_par, snd_children)

        n_snds = length(vreps)
        snds = [Vector{Int}(undef, 0) for _ in 1:n_snds]
        for v in 1:n
            push!(snds[snd_membership[v]], v)
        end

        seps = [Vector{Int}(undef, 0) for _ in 1:n_snds]
        for (ind, vrep) in enumerate(vreps)
            seps[ind] = filter(x->!(x in snds[ind]), rowvals(L)[nzrange(L, vrep)])
        end

        new(etree, vreps, snd_par, snd_children, snds, seps, postordering, collect(1:n))
    end
end

# construct a permutation as in VA section 4.6
# st the supernodes have consecutive vertices
function order_snds!(ct::CliqueTree)
    n = length(ct.perm)
    perm = zeros(Int, n)
    ind = 1
    for j in 1:length(ct.snds)
        snd_ind = ct.postordering[j]
        snd = ct.snds[snd_ind]
        len_snd = length(snd)
        perm[ind:ind+len_snd-1] .= snd
        snd .= ind:ind+len_snd-1
        ind += len_snd
    end
    ct.perm .= perm

    iperm = invperm(perm)
    for sep in ct.seps, i in 1:length(sep)
        sep[i] = iperm[sep[i]]
    end
    ct.vreps .= first.(ct.snds)
end


# Input:
#   L : sparse lower triangular matrix that rep chordal graph
function generate_clique_tree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    size(A, 1) != size(A, 2) && error(ArgumentError("A must be square"))
    if !is_chordal(A)
        error(ArgumentError("Matrix must be chordal"))
    end
    peo = !is_chordal(A; peo=1:n) ? collect(1:n) : maximum_cardinality_search(A)
    # NOTE: this sets the diagonal equal to 0
    L = tril(A[peo, peo], -1)

    return CliqueTree(L, peo)
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
    get_etree(A)
Compute the parent function of the elimination tree of a symmetric sparse matrix `A`
"""
function get_etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    (!istril(A) && !issymmetric(A)) && error(ArgumentError("A must be symmetric or lower triangular"))
    n = size(A, 1)
    par = zeros(Ti, n)
    ancestor = zeros(Ti, n)

    U = istril(A) ? sparse(A') : A
    _etree(U, par, ancestor)
    count(x->x==0, par) > 1 && error(ArgumentError("A must be connected"))
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
    n = size(L, 1)
    etree_children = get_children_from_par(etree_par)
    post_ord = get_postordering(etree_par, etree_children)
    deg⁺ = get_higher_deg(L)

    vreps = Vector{Int}(undef, 0)
    sizehint!(vreps, n)
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
            push!(vreps, v)
            snd_membership[v] = v
            u = v
        # otherwise, add v to snd(ŵ)
        else
            u = snd_membership[ŵ]
            snd_membership[v] = u
        end

        for w in etree_children[v]
            w == ŵ && continue
            # q(z) = u, a(z) = v; z is the vertex in vreps that satisfies w ∈ snd(z)
            z = snd_membership[w]
            snd_par[z] = u
            # a[z] = v
        end
    end

    # numbers the supernodes
    par = zeros(Int, length(vreps))
    @views for (ind, v) in enumerate(snd_par[vreps])
        vrep_ind = findfirst(i->vreps[i] == v, 1:length(vreps))
        if (!isnothing(vrep_ind)) par[ind] = vrep_ind end
    end

    # updates membership similarly
    # TODO: make this efficient
    for (ind, v) in enumerate(snd_membership)
        snd_membership[ind] = findfirst(i->vreps[i] == v, 1:length(vreps))
    end

    return vreps, par, snd_membership

end


# input: parent vector for tree
# returns vector v st v[i] = {k | k is a child of i}
"""
    get_children_from_par(par)
Compute the children of each vertex `v` in a tree specified by `par`
"""
function get_children_from_par(par::Vector{Int})
    n = length(par)
    children = [Vector{Int}(undef, 0) for _ in 1:n]

    for i in 1:n
        par_i = par[i]
        if par_i != 0
            push!(children[par_i], i)
        end
    end

    return children
end


# NOTE: Assumes that L has order σ = 1:n & is lower triangular w zeros on diagonal
function get_higher_deg(L::SparseMatrixCSC)
    n = size(L, 1)
    deg⁺ = zeros(Int, n)

    for i in 1:n
        deg⁺[i] = L.colptr[i+1] - L.colptr[i]
    end

    return deg⁺
end


# TODO: Iterator for cliques
function get_cliques(ct::CliqueTree)
    return [vcat(ct.snds[i], ct.seps[i]) for i in 1:length(ct.snds)]
end
