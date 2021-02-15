cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
using ChordalDecomp
using MosekTools
using SparseArrays, JuMP, LinearAlgebra, Random

function print_and_return(x, i, n)
    @info "$i/$n"
    return x
end

function block(i, n)
    A = spzeros(n, n)
    block_size = randn() > 0 ? 2 : 5
    A[i:min(i+block_size,n), i:min(i+block_size,n)] .= 1
    return A
end

function build_ΣAᵢuᵢ(m, u; rng=0)
    Random.seed!(rng)
    n = length(u)
    A = spzeros(GenericAffExpr{Float64, VariableRef}, n, n)
    for i=1:n
        Q = block(i, n)
        for (j, k, v) ∈ zip(findnz(Q)...)
            A[j,k] = u[i]*v
            # add_to_expression!(A[j, k], u[i], v)
        end
    end
    return A
end

rng = 123

## Problem Setup
n = 100
C = sparse(I(n))

m = Model(optimizer_with_attributes(
    Mosek.Optimizer,
    "QUIET" => false,
))
@variable(m, u[1:n])
A = build_ΣAᵢuᵢ(m, u; rng=rng)
@show nnz(A)
@show size(A)

## Standard Solve
@time begin
    @constraint(m, A - C ∈ PSDCone())
    @info "Finished adding constraints"
    @objective(m, Min, u[n])
    optimize!(m)
end

pstar = value.(u[n])
utrue = value.(u)

@info "termination status: $(termination_status(m))"
@info "solution status: $(primal_status(m))"
@info "optimal value: $pstar"
## Functions to do things™  (Chordal Decomp)

# TODO: this function should probably go into the package
function build_constraints!(m, A, C)
    n = size(A)[1]
    P, Cℓs, Tℓs = get_selectors(A)
    A_prime = P*A*P'
    C_prime = P*C*P'
    S = Vector{AbstractMatrix}(undef, length(Cℓs))
    for p=1:length(Cℓs)
            len_p = length(Cℓs[p])
            S[p] = @variable(m, [1:len_p, 1:len_p] in PSDCone(), base_name="Cℓ$p")
            A_prime -= Tℓs[p]'*S[p]*Tℓs[p]
    end
    for p = 1:length(Cℓs)
            @constraint(m, Tℓs[p]*A_prime*Tℓs[p]' .== Tℓs[p]*C_prime*Tℓs[p]')
    end
end


## do things™
m = Model(Mosek.Optimizer)
@variable(m, u[1:n])
A = build_ΣAᵢuᵢ(m, u; rng=rng)
build_constraints!(m, A, sparse(I(n)))
@info "Finished adding constraints"
##
@time begin
    # A_prime = build_new_and_improved_constraint!(m, A)
    # @constraint(m, A_prime .== sparse(I(n)))
    @objective(m, Min, u[n])
    optimize!(m)
end

u_chordal = value.(u)
pstar_chordal = value(u[end])
@info "termination status: $(termination_status(m))"
@info "solution status: $(primal_status(m))"
@info "pstar = $(pstar_chordal)"
