cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
using ChordalDecomp
using MosekTools
using SparseArrays, JuMP, LinearAlgebra, Random


rand_seed = 1234
n = 200
c, F, G, xstar, D = generate_random_sdp(n, rand_seed=rand_seed)


## Without decomposition
m = Model(optimizer_with_attributes(
    Mosek.Optimizer,
    "QUIET" => false,
))
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
function build_constraints!(m, A)
    P, Cℓs, Tℓs = get_selectors(A)
    A = P*A*P'
    S = Vector{AbstractMatrix}(undef, length(Cℓs))
    for p=1:length(Cℓs)
            len_p = length(Cℓs[p])
            S[p] = @variable(m, [1:len_p, 1:len_p] in PSDCone(), base_name="Cℓ$p")
            A -= Tℓs[p]'*S[p]*Tℓs[p]
    end
    for p = 1:length(Cℓs)
            @constraint(m, Tℓs[p]*A*Tℓs[p]' .== 0)
    end
end

m2 = Model(optimizer_with_attributes(
    Mosek.Optimizer,
    "QUIET" => false,
))
@variable(m2, y[1:n])
A = dropzeros(sum(F[i]*y[i] for i in 1:n) + G)
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
@info "Time with chordal is $(round(time_chordal, digits=3))s vs $(round(time_non_chordal, digits=3))s"
