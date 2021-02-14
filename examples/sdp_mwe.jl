using ChordalDecomp
using MosekTools
using SparseArrays, JuMP, LinearAlgebra, Random

function print_and_return(x, i, n)
    @info "$i/$n"
    return x
end

Random.seed!(1234)

function block(i, n)
    A = spzeros(n, n)
    block_size = randn() > -2 ? 2 : 200
    A[i:min(i+block_size,n), i:min(i+block_size,n)] .= 1
    return A
end

## Problem Setup

n = 200
m = Model(Mosek.Optimizer)
A = spzeros(GenericAffExpr{Float64, VariableRef}, n, n)
@variable(m, u[1:n])

for i=1:n
    Q = block(i, n)
    for (j, k, v) ∈ zip(findnz(Q)...)
        A[j,k] = 0.0001
        add_to_expression!(A[j, k], u[i], v)
    end
end

@show nnz(A)
@show size(A)

## Standard Solve
@time begin
    @constraint(m, A - sparse(I(n)) ∈ PSDCone())
    @info "Finished adding constraints"
    @objective(m, Min, 0)
    optimize!(m)
end




## Functions to do things™  (Chordal Decomp)
function do_the_thing™(A)
    ChordalDecomp.preprocess!(A)
    sp = sparsity_pattern(A)
    P, L = get_chordal_extension(sp)
    P = build_perm_matrix(P)
    cliques = get_cliques(L)
    cg = generate_clique_graph(cliques)
    merge_cliques!(cg; verbose=true)
    Cℓs = get_cliques(cg)
    Tℓs = make_selectors(Cℓs, n)
    return P, Cℓs, Tℓs
end


function build_new_and_improved_constraint!(m, A)
    P, Cℓs, Tℓs = do_the_thing™(A)
    A_prime = P'*A*P
    S = Vector{AbstractMatrix}(undef, length(Cℓs))
    for p=1:length(Cℓs)
            len_p = length(Cℓs[p])
            S[p] = @variable(m, [1:len_p, 1:len_p], PSD)
            # A_prime += Tℓs[p]' * S[p] * Tℓs[p]

            @constraint(m, Tℓs[p]*A_prime*Tℓs[p]' + S[p] .== sparse(I(len_p)))
    end
    return nothing
end


## do things™
@time begin
    # A_prime = build_new_and_improved_constraint!(m, A)
    # @constraint(m, A_prime .== sparse(I(n)))
    build_new_and_improved_constraint!(m, A)
    @info "Made the constraint."
    optimize!(m)
end
