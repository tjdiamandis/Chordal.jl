cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
using ChordalDecomp
using JuMP, MosekTools
using SparseArrays, LinearAlgebra, Random


rand_seed = 1234
n = 200
c, F, G, xstar, D = generate_random_sdp(n, rand_seed=rand_seed)


## Without decomposition
m = Model(optimizer_with_attributes(
    Mosek.Optimizer,
    "QUIET" => false,
))
# m = Model(Hypatia.Optimizer)
@variable(m, x[1:n])
@constraint(m, sum(F[i]*x[i] for i in 1:n) + G ∈ PSDCone())
@objective(m, Min, dot(c, x))
time_non_chordal = @elapsed optimize!(m)

xv = value.(x)
pstar = dot(c,xv)
@info "termination status: $(termination_status(m))"
@info "solution status: $(primal_status(m))"
@info "Difference with true optimal: $(norm(xv - xstar))"
@info "Optimal value: $pstar"


## Chordal Decomposition Setup
# A + S ∈ PSDCone()
function build_constraints!(m, A)
    P, Tℓs = get_selectors(A; verbose=false, ret_cliques=false)
    num_cliques = length(Tℓs)
    A = P*A*P'
    S = Vector{AbstractMatrix}(undef, num_cliques)
    for p=1:num_cliques
            len_p = size(Tℓs[p], 1)
            S[p] = sparse(@variable(m, [1:len_p, 1:len_p] in PSDCone(), base_name="Cℓ$p"))
            for (j, k, v) ∈ zip(findnz(-Tℓs[p]'*S[p]*Tℓs[p])...)
                add_to_expression!(A[j,k], v)
            end
            # @info "Build constraint $p of $num_cliques"
    end
    for p = 1:num_cliques
            @constraint(m, Tℓs[p]*A*Tℓs[p]' .== 0)
    end
end


function build_A(F, G, y)
    n = length(y)
    A = spzeros(GenericAffExpr{Float64, VariableRef}, n, n)
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

    for (j, k, v) ∈ zip(findnz(sparsity_pattern([F..., G]))...)
        A[j,k] -= 1
    end
    return A
end




m2 = Model(optimizer_with_attributes(
    Mosek.Optimizer,
    "QUIET" => false,
))
@variable(m2, y[1:n])
A = build_A(F, G, y)
# sum(F[i]*y[i] for i in 1:n) + G
nnzA = nnz(A)
@info "There are $nnzA nonzeros; density = $(round(nnzA/n^2, digits=3))"
time_constraints = @elapsed build_constraints!(m2, A)
@info "Finished adding constraints, time = $(round(time_constraints, digits=3))s"
@objective(m2, Min, dot(c,y))
@info "Finished setup"


## do things™
time_chordal = @elapsed optimize!(m2)
yv = value.(y)
pstar_chordal = dot(c,yv)

@info "termination status: $(termination_status(m))"
@info "solution status: $(primal_status(m))"
@info "difference with non-chordal: $(norm(yv - xv))"
@info "difference with true variable: $(norm(yv - xstar))"
@info "Time with chordal is $(round(time_chordal, digits=3))s vs $(round(time_non_chordal, digits=3))s"
