#=
# Semidefinite Program Decomposition
This example illustrates how to decompose a positive semidefinite constraint with
`Chordal.jl`. The decomposed constraint can then be plugged into your favorite
semidefinite program (SDP) solver.
=#

#=
## Problem Setup
We first generate a random SDP using `generate_random_sdp(n)`, which
constructs problems of the form

$$
\begin{array}{ll}
\text{minimize} &c^Tx \\
\text{subject to} & \sum_{i=1}^m F_ix_i + G \succeq 0.
\end{array}
$$
=#

#-
using Chordal
using JuMP, SCS
using SparseArrays, LinearAlgebra, Random

rand_seed = 1234
n = 500
optimizer = SCS.Optimizer

## D is the dual variable and xstar is an optimal solution
c, F, G, xstar, D = Chordal.generate_random_sdp(n, rand_seed=rand_seed);

#=
## Standard Solve
First, we build the SDP using `JuMP` and solve it using `SCS` without chordal
decomposition.
=#
#-
## Without decomposition
m = Model(optimizer)
JuMP.set_silent(m)
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

#=
## Solving after Chordal Decomposition
Next, we use chordal decomposition tools from `Chordal.jl` to decompose the PSD
constraint into smaller blocks.
=#

## Chordal Decomposition Setup
m2 = Model(optimizer)
JuMP.set_silent(m2)
@variable(m2, y[1:n])

## build_A is a helper function to construct the SDP
## A = sum(F[i]*y[i] for i in 1:n) + G
## A + S ∈ PSDCone(), where S ⪰ 0
A = Chordal.build_A(y, F, G)
nnzA = nnz(A)
@info "There are $nnzA nonzeros; density = $(round(nnzA/n^2, digits=3))"
time_constraints = @elapsed build_constraints_lmi!(m2, A)
## NOTE: Can also call build_constraints_lmi!(m2, y, F, G)
@info "Finished adding constraints, time = $(round(time_constraints, digits=3))s"
@objective(m2, Min, dot(c,y))
@info "Finished setup"

#=
We can solve the decomposed model with any SDP solver (here we use SCS).
=#
#-
## Chordal Decomposition Solve
time_chordal = @elapsed optimize!(m2)
yv = value.(y)
pstar_chordal = dot(c,yv)

@info "termination status: $(termination_status(m2))"
@info "solution status: $(primal_status(m2))"
@info "difference with non-chordal: $(norm(yv - xv))"
@info "difference with true variable: $(norm(yv - xstar))"
@info "Time with chordal is $(round(time_chordal, digits=3))s vs $(round(time_non_chordal, digits=3))s"
