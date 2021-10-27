using Chordal
using JuMP, Hypatia
using SparseArrays, LinearAlgebra, Random
using GLMakie

include("utils.jl")

# Construct Dataset
n, k = 50, 2
x, theta = swissroll(n)
show_data(x, theta)


# Construct problem data
V, VA, b = get_constraint_matrices(x, k)
C = V'*V
perm, iperm, Tls, cliques = get_selectors(sparsity_pattern([VA..., C]); verbose=false, ret_cliques=true, perm=nothing)
sp = Chordal.get_sparsity_pattern_from_cliques(cliques)
sp_inds = zip(findnz(sp)[1:2]...)

# Construct Model
model_c = Model(Hypatia.Optimizer)
Z = @variable(model_c, Z[i=1:n-1, j=1:n-1; (i,j) in sp_inds])

# Constraints
@objective(model_c, Max, Chordal.tr(C, Z))
for i in 1:length(b)
    @constraint(model_c, Chordal.tr(VA[i], Z) == b[i])
end
for (p, cl) in enumerate(cliques)
    len_cl = length(cl)
    Zp = @variable(model_c, [1:len_cl, 1:len_cl], base_name="Z$p")
    for i in 1:len_cl, j in 1:len_cl
        Zp[i,j] = Z[cl[i], cl[j]]
    end
    @constraint(model_c, Zp in PSDCone())
end
optimize!(model_c)
display(solution_summary(model_c))

# Reconstruct solution
Z_ = Chordal.reconstruct_from_sparse_varref(Z, n-1)
Z_ = 0.5*(Z_ + Z_') # Deals with numerical error

Z_factor = Chordal.minrank_completion(Z_)
W, sv, _ = svd(V*Z_factor)
yc = Diagonal(sv[1:3])*W[:, 1:3]'
show_data(yc, theta)
