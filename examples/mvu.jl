#=
# Maximum Variance Unfolding (MVU)
This example uses `Chordal.jl` to perfomr the MVU dimensionality reduction 
technique on a toy dataset.
=#

#-
using Chordal
using JuMP, Hypatia
using SparseArrays, LinearAlgebra, Random
using Plots

#=
Some helper functions are below. This section can be skipped but is included
for completeness.
=#
#-
## Generate swissroll dataset
function swissroll(n)
    ## NOTE: these could be randomized too
    theta = 1.5 .* pi .+ 2.5 .* pi .* LinRange(0, 1, n)
    z = [(i % 2)-1/2 for i in 1:n]
    return [theta' .* [cos.(theta'); sin.(theta')]; z'], theta
end

function show_data(x, theta; c=(70,70))
    Mx = maximum(abs.(x[1,:]))
    My = maximum(abs.(x[2,:]))
    Mz = maximum(abs.(x[3,:]))
    MM = max(Mx, My, Mz)
    Mx = max(Mx*1.1, MM*.9)
    My = max(My*1.1, MM*.3)
    Mz = max(Mz*1.1, MM*.3)
    return scatter(x[1,:], x[2,:], x[3,:],
        legend=false,
        zcolor=theta,
        xlims=(-Mx,Mx),
        ylims=(-My, My),
        zlims=(-Mz,Mz),
        markersize=7,
        camera=c,
    )
end

## Constructs Euclidean Distance Matrix given vectors x = [x_1 ... x_n]
function edm(x)
    G = x' * x
    dg = diag(G)
    return -2*G + dg*ones(size(G,1))' + ones(size(G,1))*dg'
end

## Given edge weights E, Eij = d(xi, xj), finds k nearest neighbors
function knn(E, k)
    nn = [Vector{Int}(undef,0) for _ in 1:size(E, 1)]
    for i in 1:size(E, 1)
        nn[i] = partialsortperm(@view(E[:,i]), 1:k+1)
    end
    return nn
end

## from nearest neighbors, constructs the MVU edge list
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

## Returns edge list for k nearest neighbors given vectors X
function edges(x, k, clique=true)
    nn = knn(edm(x), k)
    return reduce(vcat, [knn_edges(nn[i], clique) for i in 1:length(nn)])
end

## Constructs V matrix (basis)
function get_V(n)
    return spdiagm(n, n-1, 0 => ones(n-1), -1 => -ones(n-1))
end

## Constructs constrain matrices (as in Vandenberghe and Andersen pg 401)
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

## Returns a SparseMatrixCSC with values of Zref
##   Zref is a SparseAxisArray, which functions like a dictionary
##   - Requires optimizer to have been run
function reconstruct_from_sparse_varref(Zref, n)
    Zcv = value.(Zref)
    Z_uncomp = spzeros(Float64, n, n)
    for ind in eachindex(Zcv)
        Z_uncomp[ind[1], ind[2]] = Zcv[ind]
    end
    return Z_uncomp
end

#=
## Problem Setup
First, we construct a "swiss roll" test dataset, which is 3d but has intrinsic 
dimension 2
=#
#-
## Construct Dataset
n, k = 80, 2
x, theta = swissroll(n)
show_data(x, theta)

#=
Next, we setup the MVU SDP
=#

#-
## Construct problem data
V, VA, b = get_constraint_matrices(x, k)
C = V'*V;


#=
## Chordal Decomposition and Solve
=#
## Chordal setup
_, _, _, cliques = get_selectors(sparsity_pattern([VA..., C]); verbose=false, ret_cliques=true, perm=nothing)
sp = Chordal.get_sparsity_pattern_from_cliques(cliques)
sp_inds = zip(findnz(sp)[1:2]...)

## -- Construct Model --
model = Model(Hypatia.Optimizer)
Z = @variable(model, Z[i=1:n-1, j=1:n-1; (i,j) in sp_inds])

## Objective
@objective(model, Max, 1/length(sp_inds) * Chordal.tr(C, Z))

## Distrance-preserving equality constraints 
for i in 1:length(b)
    @constraint(model, Chordal.tr(VA[i], Z) == b[i])
end

## Chordal decomposition constraints (replace big PSD constraints with many small ones)
for cl in cliques
    len_cl = length(cl)
    Zp = Matrix{VariableRef}(undef, len_cl, len_cl)
    for i in 1:len_cl, j in 1:len_cl
        Zp[i,j] = Z[cl[i], cl[j]]
    end
    @constraint(model, Zp in PSDCone())
end

## Optimize
optimize!(model)
solution_summary(model)
#=
## Solution Reconstruction
=#
## Reconstruct solution with min rank PSD completion
Z_ = Chordal.reconstruct_from_sparse_varref(Z, n-1)
Z_ = 0.5*(Z_ + Z_') # Deals with numerical error
Z_factor = Chordal.minrank_completion(Z_)
W, sv, _ = svd(V*Z_factor)
yc = Diagonal(sv[1:3])*W[:, 1:3]'
show_data(yc, theta; c =(40,70))