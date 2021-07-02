cd(@__DIR__)
Pkg.activate("..")
using ChordalDecomp
using JuMP, Hypatia, COSMO
using SparseArrays, LinearAlgebra, Random
const CD = ChordalDecomp

rand_seed = 1234
n = 200
# optimizer = Hypatia.Optimizer
optimizer = optimizer_with_attributes(COSMO.Optimizer, "decompose" => false)

# min c^Tx st ΣF_ix_i + G ∈ PSDCONE()
# D is dual variable
c, F, G, xstar, D = CD.generate_random_sdp(n, rand_seed=rand_seed)

# ------------------------------------------------------
# -------------------- LMI Form SDP --------------------
# ------------------------------------------------------
## Without decomposition
m = Model(optimizer)
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
m2 = Model(optimizer)
@variable(m2, y[1:n])

# A = sum(F[i]*y[i] for i in 1:n) + G
# A + S ∈ PSDCone(), where S ⪰ 0
A = CD.build_A(y, F, G)
nnzA = nnz(A)
@info "There are $nnzA nonzeros; density = $(round(nnzA/n^2, digits=3))"
time_constraints = @elapsed build_constraints_lmi!(m2, A)
# NOTE: Can also call build_constraints_lmi!(m2, y, F, G)
@info "Finished adding constraints, time = $(round(time_constraints, digits=3))s"
@objective(m2, Min, dot(c,y))
@info "Finished setup"

## Chordal Decomposition Solve
time_chordal = @elapsed optimize!(m2)
yv = value.(y)
pstar_chordal = dot(c,yv)

@info "termination status: $(termination_status(m2))"
@info "solution status: $(primal_status(m2))"
@info "difference with non-chordal: $(norm(yv - xv))"
@info "difference with true variable: $(norm(yv - xstar))"
@info "Time with chordal is $(round(time_chordal, digits=3))s vs $(round(time_non_chordal, digits=3))s"
