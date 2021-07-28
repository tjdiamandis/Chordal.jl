cd(@__DIR__)
Pkg.activate("..")
using Chordal
using LinearAlgebra, SparseArrays

const CD = Chordal

## Small Example
# from eq (1041) https://ccrma.stanford.edu/~dattorro/EDM.pdf
D = sparse([
    0.0 1 5 0;
    1 0 4 1;
    5 4 0 1;
    0 1 1 0
])
@show D

D_complete = sparse(edm_completion(D))
println("\nCompletion:")
@show D_complete


## Larger Random Example
# Pattern from Figure 4.2 in [VA15]
n = 17
r = 3
ijs = [(3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,3), (5,4), (8,7),
       (9,5), (9,6), (9,7), (9,8), (11,10), (13,10), (13,11), (13,12),
       (14,10), (14,11), (14,12), (14,13),
       (15, 1), (15,3), (15,4), (15,5), (15,7), (15,8), (15,9),
       (16,5), (16,6), (16,9), (16,12), (16,13), (16,14), (16,15),
       (17,10), (17,11), (17,12), (17,13), (17,14), (17,15), (17,16)]
append!(ijs, [(i,i) for i in 1:n])
II, JJ = CD.unzip(ijs)
sp = sparse(II, JJ, ones(length(II)))
sp = sp + tril(sp)'

V = randn(n, r)
VV = V*V'
Afull = diag(VV)*ones(n)' + ones(n)*diag(VV)' - 2VV
A = sp .* Afull

@show CD.is_edm(Afull)
@show CD.is_edm(edm_completion(A))
