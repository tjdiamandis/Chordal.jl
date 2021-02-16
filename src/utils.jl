using Random

function build_perm_matrix(p)
    isnothing(p) && return I
    n = length(p)
    P = spzeros(n, n)
    for i in 1:n
        P[i, p[i]] = 1.0
    end
    return P
end


# Build from KKT conditions:
# # KKT Conditions
# Fx = sum(F[i]*xstar[i] for i in 1:n) + G
# all(eigvals(Matrix(Fx)) .>= 0)
# all(eigvals(D) .>= 0)
# all([tr(F[i]*D) .== c[i] for i in 1:n])
# tr(Fx*D) .== 0
function generate_random_sdp(n; rand_seed=0)
    Random.seed!(rand_seed)

    # Can assume diagonal WLOG (move around ortho matrix in trace)
    D = diagm(1 .+ rand(n))
    F = Vector{SparseMatrixCSC{Float64}}(undef, n)
    c = Vector{Float64}(undef, n)
    for i in 1:n
        F[i] = spzeros(n, n)
        block_size = randn() > 1 ? 2 : 10
        F[i][i:min(i+block_size,n), i:min(i+block_size,n)] .= 1
        c[i] = tr(D*F[i])
    end
    xstar = rand(n)
    Fx = sum(F[i]*xstar[i] for i in 1:n)
    G = -Fx

    return c, F, G, xstar, D
end
