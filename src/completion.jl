
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

# Logarithmic Barriers for Sparse Matrix Cones
# Andersen, Dahl, Vandenberghe
# Algorithm 4.2
# NOTE: This is very inefficient -- needs some significant work
function maxdet_completion_etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && error(ArgumentError("A must be symmetric"))
    n = size(A, 1)

    sp = sparsity_pattern(A)
    L = get_chordal_extension(sp; perm=nothing, verbose=false)[3]
    # L[diagind(L)] .= 1
    etree_par = etree(L)
    etree_children = get_children_from_par(etree_par)
    post_ord = get_postordering(etree_par, etree_children)

    L_chol = spzeros(n,n)
    L_chol[diagind(L_chol)] .= 1
    D_chol = Diagonal(zeros(n))
    V = Vector{Matrix{Float64}}(undef, n)

    for i in n:-1:1
        j = post_ord[i]
        Ij = rowvals(L)[nzrange(L, j)]

        Vj = (i == n) ? [A[end, end]] : V[j]
        if i != n
            L_chol[Ij, j] = (-Vj) \ Vector{Float64}(A[Ij, j])
            D_chol[j,j] = 1/(A[j,j] + dot(A[Ij, j], L_chol[Ij, j]))
        else
            D_chol[j,j] = 1 / A[j,j]
        end

        for ch in etree_children[j]
            Ich = rowvals(L)[nzrange(L, ch)]

            nv = length(Ij) + 1
            tmp = zeros(nv, nv)
            tmp[1,1] = A[j,j]
            tmp[2:end, 1] .= A[Ij, j]
            tmp[1, 2:end] .= A[Ij, j]
            tmp[2:end, 2:end] .= Vj

            E_Jj_Ich = spzeros(length(Ij) + 1, length(Ich))
            Jj = vcat([j], Ij)
            for ii in 1:length(Jj), jj in 1:length(Ich)
                if Jj[ii] == Ich[jj]
                    E_Jj_Ich[ii, jj] = 1.0
                end
            end

            V[ch] = E_Jj_Ich' * tmp * E_Jj_Ich
        end
    end

    return L_chol, D_chol
end



# Logarithmic Barriers for Sparse Matrix Cones
# Andersen, Dahl, Vandenberghe
# Algorithm 7.3
# NOTE: This is very inefficient -- needs some significant work
function maxdet_completion_factors(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && error(ArgumentError("A must be symmetric"))
    n = size(A, 1)

    sp = sparsity_pattern(A)
    L = get_chordal_extension(sp; perm=nothing, verbose=false)[3]
    etree_par = etree(L)
    vreps, snd_par, snd_membership = max_supernode_etree(L, etree_par)
    n_snds = length(vreps)
    snds = [Vector{Int}(undef, 0) for _ in 1:n_snds]
    for v in 1:n
        push!(snds[snd_membership[v]], v)
    end
    snd_children = get_children_from_par(snd_par)
    post_ord = get_postordering(snd_par, snd_children)


    L_chol = spzeros(n,n)
    L_chol[diagind(L_chol)] .= 1
    D_chol = spzeros(n,n)
    V = Vector{Matrix{Float64}}(undef, n_snds)

    for i in n_snds:-1:1
        vrep_ind = post_ord[i]
        j = vreps[vrep_ind]
        Jj = vcat([j], rowvals(L)[nzrange(L, j)])
        Nj = snds[vrep_ind]
        Aj = filter(x->!(x in Nj), Jj)

        Vj = (i == n_snds) ? [A[end, end]] : V[vrep_ind]
        if i != n_snds
            L_chol[Aj, Nj] = pinv(-Vj) * A[Aj, Nj]
            D_chol[Nj,Nj] = inv(Matrix(A[Nj,Nj] + A[Aj, Nj]'*L_chol[Aj, Nj]))
        else
            D_chol[Nj,Nj] = inv(Matrix(A[Nj,Nj]))
        end

        for ch_ind in snd_children[vrep_ind]
            Nch = snds[ch_ind]
            ch = vreps[ch_ind]
            Jch = vcat([ch], rowvals(L)[nzrange(L, ch)])
            Ach = filter(x->!(x in Nch), Jch)

            nv = length(Jj)
            tmp = zeros(nv, nv)
            len_Nj = length(Nj)
            len_Aj = length(Aj)
            tmp[1:len_Nj,1:len_Nj] = A[Nj,Nj]
            if len_Aj > 0
                tmp[len_Nj+1:end, 1:len_Nj] .= A[Aj, Nj]
                tmp[1:len_Nj, len_Nj+1:end] .= A[Aj, Nj]'
                tmp[len_Nj+1:end, len_Nj+1:end] .= Vj
            end
            E_Jj_Ach = spzeros(length(Jj), length(Ach))
            for ii in 1:length(Jj), jj in 1:length(Ach)
                if Jj[ii] == Ach[jj]
                    E_Jj_Ach[ii, jj] = 1.0
                end
            end
            V[ch_ind] = E_Jj_Ach' * tmp * E_Jj_Ach
        end
    end

    return L_chol, D_chol
end


function minrank_completion(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && error(ArgumentError("A must be symmetric"))
    n = size(A, 1)

    sp = sparsity_pattern(A)
    _, _, L = get_chordal_extension(sp; perm=nothing)
    etree_par = etree(L)
    vreps, snd_par, snd_membership = max_supernode_etree(L, etree_par)
    n_snds = length(vreps)
    snds = [Vector{Int}(undef, 0) for _ in 1:n_snds]
    for v in 1:n
        push!(snds[snd_membership[v]], v)
    end
    snd_children = get_children_from_par(snd_par)
    post_ord = get_postordering(snd_par, snd_children)

    # Determine rank
    r = maximum([rank(Matrix(A[c,c]), rtol=1e-10) for c in get_cliques(L)])

    Y = zeros(n, r)
    for j in n_snds:-1:1
        vrep_ind = post_ord[j]
        vrep = vreps[vrep_ind]

        # col_j = vcat([vrep], rowvals(L)[nzrange(L, vrep)])
        ν = snds[vrep_ind]
        α = filter(x->!(x in ν), rowvals(L)[nzrange(L, vrep)])
        col_j = vcat(ν, α)
        @show ν, α

        dd, VV = eigen(Matrix(@view(A[col_j,col_j])), sortby=x->-real(x))
        r_ = min(length(dd), r)
        Z = VV[:,1:r_]*Diagonal(sqrt.(max.(real.(dd[1:r_]), 0.0)))


        if j == n_snds
            Y[ν, 1:r_] .= Z[1:length(ν), :]
            continue
        end



        U = @view(Z[1:length(ν), 1:r_])
        V = @view(Z[length(ν)+1:end, 1:r_])

        # W_Y, Σ_Y, Q_Y = svd(Y[α, :], full=true, alg = LinearAlgebra.QRIteration())
        # W_V, Σ_V, Q_V = svd(V, full=true, alg = LinearAlgebra.QRIteration())
        #
        # Q = Q_V*Q_Y'
        # # @show Y[α, :]
        # if !all(Y[α, :] .≈ V*Q)
        #     d = sign.((Y[α, :]*Q_Y)[1,:]) .* sign.((V*Q_V)[1,:])
        #     Q = Q_V*Diagonal(d)*Q_Y'
        # end
        # Q2, Σ, Q1 = svd(V'*Y[α, 1:r_], full=true, alg = LinearAlgebra.QRIteration())
        Q2, Σ, Q1 = svd(V'*Y[α, 1:r_], full=true)
        Q = Q2*Q1'
        Y[ν, 1:r_] .= U*Q

    end
    return Y
end
