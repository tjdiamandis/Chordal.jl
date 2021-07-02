
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
