#=
# Euclidean Distance Matrix Completion
This example illustrates how to use Chordal.jl to complete a partially-specified
Euclidean Distance Matrix (or determine no such completion exists).
=#

#=
## Euclidean Distance Matrices
Given $N$ vectors $x_1,\dots,x_n$ in $\mathbb{R}^n$, the associated Euclidean
Distance Matrix (EDM) $D \in \mathbb{R}^{N \times N}$ records the squared distance
between each pair $x_i, x_j$. Specifically,
$$
D_{ij} = \|x_i - x_j\|_2^2.
$$

Given a partially-specified EDM $D$ which has a chordal sparsity pattern, we would
like to find the missing entries.
=#

using Chordal
using LinearAlgebra, SparseArrays


#=
## Small Example
This small example is taken from Jon Dattorro's [Convex Optimization & Euclidean Distance Geometry](https://ccrma.stanford.edu/~dattorro/EDM.pdf)
eq (1041). First, we construct the partially-specified EDM.
=#
D = sparse([
    0.0 1 5 0;
    1 0 4 1;
    5 4 0 1;
    0 1 1 0
])
@show D

## Next, we complete the matrix
D_complete = sparse(edm_completion(D))
println("\nCompletion:")
@show D_complete


#=
## Larger Random Example

First, we generate an EDM from 17 vectors in $\mathbb{R}^3}.
=#
n = 17
r = 3
V = randn(n, r)
VV = V*V'
D_full = diag(VV)*ones(n)' + ones(n)*diag(VV)' - 2VV

#=
Next, we remove entries. and complete the
The sparsity pattern is from Figure 4.2 in Lieven Vandenberghe and Martin S. Andersen's
[Chordal Graphs and Semidefinite Optimization](https://www.seas.ucla.edu/~vandenbe/publications/chordalsdp.pdf).
=#
ijs = [(3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,3), (5,4), (8,7),
       (9,5), (9,6), (9,7), (9,8), (11,10), (13,10), (13,11), (13,12),
       (14,10), (14,11), (14,12), (14,13),
       (15, 1), (15,3), (15,4), (15,5), (15,7), (15,8), (15,9),
       (16,5), (16,6), (16,9), (16,12), (16,13), (16,14), (16,15),
       (17,10), (17,11), (17,12), (17,13), (17,14), (17,15), (17,16)]
append!(ijs, [(i,i) for i in 1:n])
II, JJ = Chordal.unzip(ijs)
sp = sparse(II, JJ, ones(length(II)))
sp = sp + tril(sp)'

## Remove entries s.t. the remaining entries have a chordal sparsity pattern
D = sp .* D_full

## Complete the matrix
@show Chordal.is_edm(edm_completion(D))
