using Random
using LightGraphs: SimpleGraph, is_connected

# TODO: this should probably be a concrete type <: AbstractMatrix that operates
# as an indexer
function build_perm_matrix(p)
    isnothing(p) && return I
    n = length(p)
    P = spzeros(n, n)
    for i in 1:n
        P[i, p[i]] = 1.0
    end
    return P
end


# Generates dual form:
#   min c^Tx st ΣF_ix_i + G ∈ PSDCONE()
#   Build from KKT conditions:
#       Fx = sum(F[i]*xstar[i] for i in 1:n) + G
#       all(eigvals(Matrix(Fx)) .>= 0)
#       all(eigvals(D) .>= 0)
#       all([0 .== c[i] - tr(F[i]*D) for i in 1:n])
#       tr(Fx*D) .== 0
function generate_random_sdp(n; rand_seed=0)
    Random.seed!(rand_seed)

    D = diagm(1 .+ rand(n))
    F = Vector{SparseMatrixCSC{Float64, Int}}(undef, n)
    c = Vector{Float64}(undef, n)
    for i in 1:n
        F[i] = spzeros(n, n)
        block_size = randn() < 1.5 ? 2 : 10
        F[i][i:min(i+block_size,n), i:min(i+block_size,n)] .= 1
        c[i] = tr(D*F[i])
    end
    xstar = rand(n)
    Fx = sum(F[i]*xstar[i] for i in 1:n)
    G = -Fx

    return c, F, G, xstar, D
end


unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

is_separable(sp::SparseMatrixCSC) = !is_connected(SimpleGraph(sp))
