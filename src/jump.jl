# dispatch on matrix inner product to handle sparse JuMP container
import LinearAlgebra.tr
function tr(X::SparseMatrixCSC, Y::JuMP.Containers.SparseAxisArray)
    return sum(v*Y[j,k] for (j,k,v) in zip(findnz(X)...))
end


# A + S ∈ PSDCone(), where S ⪰ 0
function build_constraints_lmi!(m::JuMP.Model, A::SparseMatrixCSC{JuMP.AffExpr, Int}; verbose=false)
    perm, iperm, Tls = get_selectors(A; verbose=verbose, ret_cliques=false)
    num_cliques = length(Tls)
    A = A[perm, perm]
    S = Vector{AbstractMatrix}(undef, num_cliques)
    for p=1:num_cliques
            len_p = size(Tls[p], 1)
            S[p] = sparse(@variable(m, [1:len_p, 1:len_p] in PSDCone(), base_name="Cℓ$p"))
            for (j, k, v) ∈ zip(findnz(-Tls[p]'*S[p]*Tls[p])...)
                add_to_expression!(A[j,k], v)
            end
            # @info "Build constraint $p of $num_cliques"
    end
    for p = 1:num_cliques
            @constraint(m, Tls[p]*A*Tls[p]' .== 0)
    end
end

function build_constraints_lmi!(
        m::JuMP.Model,
        y::Vector{VariableRef},
        F::AbstractVector{SparseMatrixCSC{T, S}},
        G::SparseMatrixCSC{T, S};
        verbose=false
) where {T <: Number, S <: Integer}
    build_constraints_primal!(m, build_A(y, F, G); verbose=verbose)
end


# A = F_1 y_1 + F_2 y_2 + ... + F_n y_n + G
function build_A(y::Vector{JuMP.VariableRef}, F::AbstractVector{SparseMatrixCSC{T, S}}, G::SparseMatrixCSC{T, S}) where {T <: Number, S <: Integer}
    n = length(y)
    A = spzeros(GenericAffExpr{Float64, VariableRef}, n, n)

    # This is a hack that seems to work...
    for (j, k, v) ∈ zip(findnz(sparsity_pattern([F..., G]))...)
        A[j,k] = 1
    end

    for i in 1:n
        for (j, k, v) ∈ zip(findnz(F[i])...)
            add_to_expression!(A[j,k], F[i][j,k], y[i])
        end
    end
    for (j, k, v) ∈ zip(findnz(G)...)
        add_to_expression!(A[j,k], G[j,k])
    end

    # hack part II
    for (j, k, v) ∈ zip(findnz(sparsity_pattern([F..., G]))...)
        A[j,k] -= 1
    end
    return A
end


function get_sparsity_pattern_from_cliques(cliques)
    n = maximum(maximum.(cliques))

    sp = spzeros(Float64, n, n)
    for c in cliques
        for i in c, j in c
            sp[i,j] = 1.0
        end
    end

    return sp
end


function reconstruct_from_sparse_varref(Zref, n)
    Zcv = value.(Zref)
    Z_uncomp = spzeros(Float64, n, n)
    for ind in eachindex(Zcv)
        Z_uncomp[ind[1], ind[2]] = Zcv[ind]
    end
    return Z_uncomp
end




#FIXME
function build_constraints_standard!(
        model::JuMP.Model,
        X,
        F::AbstractVector{SparseMatrixCSC{T, S}},
        G::SparseMatrixCSC{T, S},
        c::Vector{T};
        verbose=false
) where {T <: Number, S <: Integer}
    perm, iperm, Tls, cliques = get_selectors(sparsity_pattern([F..., G]); verbose=verbose, ret_cliques=true)
    num_cliques = length(Tls)

    G .= G[perm, perm]
    for i in 1:length(F)
        F[i] .= F[i][perm,perm]
    end

    for i in 1:length(c)
        @constraint(model, tr(F[i], X) == c[i])
    end

    for (p, cl) in enumerate(cliques)
        len_cl = length(cl)
        Xp = @variable(model, [1:len_cl, 1:len_cl], base_name="X$p")
        for i in 1:len_cl, j in 1:len_cl
            Xp[i,j] = Z[cl[i], cl[j]]
        end
        @constraint(model, Xp in PSDCone())
    end
end
