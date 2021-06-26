
# VA Algorithm 10.2, pg 362
function maxdet_completion(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && error(ArgumentError("A must be symmetric"))
    n = size(A, 1)

    sp = sparsity_pattern(A)
    perm, iperm, L = get_chordal_extension(sp)
    if isnothing(perm)
        perm = 1:n
        iperm = 1:n
    end
    etree_par = etree(L)
    vreps, snd_par, snd_membership = max_supernode_etree(L, etree_par)
    n_snds = length(vreps)
    snds = [Vector{Int}(undef, 0) for _ in 1:n_snds]
    for v in 1:n
        push!(snds[snd_membership[v]], v)
    end
    snd_children = get_children_from_par(snd_par)
    post_ord = get_postordering(snd_par, snd_children)


    # construct W
    # W = Matrix{Tv}(undef, n, n)
    # for col in 1:n, k in nzrange(A, col)
    #     W[iperm[rowvals(A)[k]], iperm[col]] = nonzeros(A)[k]
    # end
    W = Matrix(A[perm, perm])

    for j in n_snds-1:-1:1
        vrep_ind = post_ord[j]
        vrep = vreps[vrep_ind]

        ν = snds[vrep_ind]
        col_nz = rowvals(L)[nzrange(L, vrep)]
        α = filter(x->!(x in ν), col_nz)
        η = filter(x->!(x in col_nz), vrep+1:n)

        # cache = W_αα^† * W_αν
        cache = W[α, α] \ W[α, ν]
        # res = W_ηα * cache = W_ηα * W_αα^† * W_αν
        @views mul!(W[η, ν], W[η, α], cache)
        W[ν, η] .= @view(W[η, ν])'
    end

    # Apply inverse permutation
    return W[iperm, iperm]
end
