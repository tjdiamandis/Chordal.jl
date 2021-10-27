# Generate swissroll dataset
function swissroll(n)
    # NOTE: these could be randomized too
    theta = 1.5 .* pi .+ 2.5 .* pi .* LinRange(0, 1, n);
    z = [(i % 2) for i in 1:n]
    return [theta' .* [cos.(theta'); sin.(theta')]; z'], theta;
end

function show_data(x, theta)
    fig = Figure(resolution=(1000,1000))
    ax = LScene(fig)
    GLMakie.scatter!(Point3f0.(eachcol(x)), color = theta, markersize=500)
    # GLMakie.scatter!(Point3f0.(eachcol(x)), color = theta)
    fig[1, 1] = ax
    display(fig)
end


# Constructs Euclidean Distance Matrix given vectors x = [x_1 ... x_n]
function edm(x)
    G = x' * x
    dg = diag(G)
    return -2*G + dg*ones(size(G,1))' + ones(size(G,1))*dg'
end


# Given edge weights E, Eij = d(xi, xj), finds k nearest neighbors
function knn(E, k)
    nn = [Vector{Int}(undef,0) for _ in 1:size(E, 1)]
    for i in 1:size(E, 1)
        nn[i] = partialsortperm(@view(E[:,i]), 1:k+1)
    end
    return nn
end


# from nearest neighbors, constructs the MVU edge list
function knn_edges(nn, clique=true)
    k = length(nn) - 1
    E = CartesianIndex{2}[]

    for i in 1:k+1
        for j in i+1:k+1
            i == j && continue
            !clique && i > 1 && continue

            push!(E, CartesianIndex(nn[i], nn[j]))
        end
    end
    return E
end


# Returns edge list for k nearest neighbors given vectors X
function edges(x, k, clique=true)
    nn = knn(edm(x), k)
    return reduce(vcat, [knn_edges(nn[i], clique) for i in 1:length(nn)])
end


# Constructs V matrix (basis)
function get_V(n)
    return spdiagm(n, n-1, 0 => ones(n-1), -1 => -ones(n-1))
end


# Constructs constrain matrices as described in section 3.1
function get_constraint_matrices(x, k)
    n = size(x, 2)
    D = edm(x)

    E = edges(x, k)
    Ei = [CartesianIndex(e[1], e[1]) for e in E]
    Ej = [CartesianIndex(e[2], e[2]) for e in E]

    b = [D[e] for e in E]
    A = [sparse(
            [e[1], e[2], e[1], e[2]],
            [e[1], e[2], e[2], e[1]],
            [1, 1, -1, -1], n, n
         ) for e in E]
    V = get_V(n)
    VA = [V'*Ai*V for Ai in A]

    return V, VA, b
end


# Returns a SparseMatrixCSC with values of Zref
#   Zref is a SparseAxisArray, which functions like a dictionary
#   - Requires optimizer to have been run
function reconstruct_from_sparse_varref(Zref, n)
    Zcv = value.(Zref)
    Z_uncomp = spzeros(Float64, n, n)
    for ind in eachindex(Zcv)
        Z_uncomp[ind[1], ind[2]] = Zcv[ind]
    end
    return Z_uncomp
end