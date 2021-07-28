cd(@__DIR__)
Pkg.activate("..")
using Chordal
using JuMP, Hypatia, COSMO
using SparseArrays, LinearAlgebra, Random

const CD = Chordal

rand_seed = 1234
n = 200
# optimizer = Hypatia.Optimizer
optimizer = optimizer_with_attributes(COSMO.Optimizer, "decompose" => false)

# min c^Tx st ΣF_ix_i + G ∈ PSDCONE()
# D is dual variable
c, F, G, ystar, D = CD.generate_random_sdp(n, rand_seed=rand_seed)

# ------------------------------------------------------
# ------------------ Standard Form SDP -----------------
# ------------------------------------------------------
## Without decomposition
m3 = Model(optimizer)
sparse(@variable(m3, X[1:n, 1:n] in PSDCone()))
for i in 1:n
    @constraint(m3, sum(F[i] .* X) - c[i] == 0)
end
@objective(m3, Max, -sum(G .* X))
time_non_chordal = @elapsed optimize!(m3)

Xv = value.(X)
@info "termination status: $(termination_status(m3))"
@info "solution status: $(primal_status(m3))"
@info "Optimal value: $(round(objective_value(m3), digits=3))"


## Chordal Setup
m4 = Model(optimizer)
perm, iperm, Tls, cliques = get_selectors(sparsity_pattern([F..., G]); verbose=false, ret_cliques=true)
sp = CD.get_sparsity_pattern_from_cliques(cliques)
sp_inds = zip(findnz(sp)[1:2]...)
X = @variable(m4, X[i=1:n, j=1:n; (i,j) in sp_inds])

# TODO: Package this up into a function?
for i in 1:length(c)
    @constraint(m4, CD.tr(F[i][perm, perm], X) == c[i])
end

for (p, cl) in enumerate(cliques)
    len_cl = length(cl)
    Xp = @variable(m4, [1:len_cl, 1:len_cl], base_name="X$p")
    for i in 1:len_cl, j in 1:len_cl
        Xp[i,j] = X[cl[i], cl[j]]
    end
    @constraint(m4, Xp in PSDCone())
end

@objective(m4, Max, -CD.tr(G[perm, perm], X))

time_chordal = @elapsed optimize!(m4)
@info "termination status: $(termination_status(m4))"
@info "solution status: $(primal_status(m4))"
match = objective_value(m4) ≈ objective_value(m3)
@info "Optimal value: $(round(objective_value(m4), digits=3)), equals: $(match)"
@info "Time with chordal is $(round(time_chordal, digits=3))s vs $(round(time_non_chordal, digits=3))s"


X_ = CD.reconstruct_from_sparse_varref(X, n)
@assert maximum(X_ - X_') < 1e-8
X_ = 0.5*(X_ + X_')
Xcomp = maxdet_completion(X_)[iperm, iperm]
@assert sum(-G .* Xcomp) ≈ objective_value(m4)
