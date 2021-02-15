cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
using ChordalDecomp
using Random
using LinearAlgebra, SparseArrays
using JuMP, MosekTools
import Plots, GraphPlot

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
# Example from A Clique Graph Based Merging Strategy for Decomposable SDPs
# by Michael Garstka, Mark Cannon, Paul Goulart
#       See Figures 1 and 3
#       Note: This does not include chordal extension (already chordal)

## Setup SDP
Random.seed!(0)
n, m = 9, 4
C = sparse(I(n))

nonzero_inds = vcat([
        (3,1),
        (3,2),
        (5,4),
        (6,1),
        (6,3),
        (7,3),
        (7,6),
        (8,3),
        (8,4),
        (8,5),
        (8,6),
        (8,7),
        (9,6),
        (9,7),
        (9,8)
    ],
    [(i,i) for i in 1:n]
)
u_match = rand(1:m, length(nonzero_inds))

function build_A(nonzero_inds, u, u_match)
        m = length(u)
        A = spzeros(GenericAffExpr{Float64, VariableRef}, n, n)
        for k in 1:length(nonzero_inds)
                i, j = nonzero_inds[k]
                A[i,j] = A[j,i] = 0.0001
                idx = u_match[k]
                add_to_expression!(A[i,j], u[idx])
                add_to_expression!(A[j,i], u[idx])
        end
        return A
end

model = Model(Mosek.Optimizer)
@variable(model, u[1:m])
A = build_A(nonzero_inds, u)
@constraint(model, A - C ∈ PSDCone())


## Solve problem
@objective(model, Min, b'*u)
optimize!(model)
@info "Terminated with status $(termination_status(model))"
pstar = value(u[m])
ustar = value.(u)
@info "p⋆ = $pstar, u⋆ = $ustar"


## Setup Chordal SDP (directly in form of paper)
sp_pattern = sparse(unzip(nonzero_inds)..., ones(length(nonzero_inds)))
sp_pattern = sp_pattern + sp_pattern'
P, Cℓs, Tℓs = get_selectors(sp_pattern)

model = Model(Mosek.Optimizer)
@variable(model, u[1:m])
A = build_A(nonzero_inds, u, u_match)
S = Vector{AbstractMatrix}(undef, length(Cℓs))
for p=1:length(Cℓs)
    len_p = length(Cℓs[p])
    S[p] = @variable(model, [1:len_p, 1:len_p] in PSDCone(), base_name="Cℓ_$p")
    global A -= Tℓs[p]'*S[p]*Tℓs[p]
end

# Exploits the sparsity in the constraints
for p = 1:length(Cℓs)
        @constraint(model, Tℓs[p]*A*Tℓs[p]' .== Tℓs[p]*C*Tℓs[p]')
end


## Solve problem
@objective(model, Min, b'*u)
optimize!(model)
@info "Terminated with status $(termination_status(model))"
pstar_c = value(u[m])
ustar_c = value.(u)
@info "p⋆ = $pstar_c, u⋆ = $ustar_c"
