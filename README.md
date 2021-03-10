# ChordalDecomp

[![Build Status](https://github.com/tjdiamandis/ChordalDecomp.jl/workflows/CI/badge.svg)](https://github.com/tjdiamandis/ChordalDecomp.jl/actions)
[![Coverage](https://codecov.io/gh/tjdiamandis/ChordalDecomp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tjdiamandis/ChordalDecomp.jl)


## TODOs
- [ ] Performance and type stability
- [ ] Unit tests -- ensure data structure change worked
- [ ] PSD matrix completion
- [ ] Implement maximal cliques algorithm to avoid LightGraphs dep?
- [ ] Implement different pre-cholesky reorderings
- [ ] Implement clique tree algorithms
- [ ] Allow user-specified weight algorithm for clique graph merging
- [ ] Test with larger problems



## Chordal Decomposition Steps
1. Reorder the vertex set `V = {1, 2, ..., n}` via some heuristics, e.g.,
    - multi-level nested dissection
    - minimum degree
        - P. R. Amestoy, T. A. Davis, and I. S Duﬀ. An approximate minimum degree ordering algorithm. SIAM Journal on Matrix Analysis and Applications, 17(4):886905, 1996.
    - Can be handled by step 2 below with certain Cholesky implementations
2. Symbolic Cholesky factorization
    - This provides the **chordal extension** (get edge set `E` from `L`)
3. Merge cliques into larger blocks
    - form clique graph or tree
    - choose heuristic merging algorithm
4. Index selector matrices

### Clique Merging Strategies
#### Clique Tree Based (Vandenberghe paper)
1. Form the clique tree
    - **Clique intersection property (CIP):** For each pair of cliques, their intersection is contained in every vertex on the unique path connection them in the tree
    - Can be efficiently computed from the chordal graph, see:
        - Blair, J.R.S., Peyton, B. (1993):An introduction to chordal graphs and clique trees. In: George,A., Gilbert, J.R., Liu, J.W.H. (eds.), Graph Theory and Sparse Matrix Computation. Springer-Verlag, pp 1–29
        - Tarjan, R.E., Yannakakis, M. (1984): Simple linear-time algorithms to test chordality of graphs, test acyclicity of hypergraphs, and selectively reduce acyclic hypergraphs. SIAM J. Comput. 13(3), 566–579
        - Lewis, J.G., Peyton, B.W., Pothen, A. (1989): A fast algorithm for reordering sparse matrices for parallel factorization. SIAM J. Sci. Statist. Comput. 10(6), 1146–1173
        - https://www.cs.cornell.edu/info/people/csun/papers/cct.ps
2. Pick an arbitrary clique as the root. Form a topological ordering (number ea. child before parent) that satisfies the running intersection property to get a perfect elimination ordering

#### Clique Graph based ([Garstka et al.](https://arxiv.org/abs/1911.05615), implemented in [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl))
1. Create clique graph
    - Compute weight function e(C_i, C_j) = w_{ij}
    - If eigenvalue is dominant step, can use `e(Ci, Cj) = |Ci|^3 + |Cj|^3 − |Ci ∪ Cj|^3`.
    - See Section 3 of paper
2. Algorithm (greedy):
    1. Merge cliques connected by edge with highest positive weight
    2. Recompute weights
    3. If there are more positive weights, repeat. Else, break.


## References
- Lieven Vandenberghe and Martin Andersen's [Chordal Graphs and Semidefinite Optimization](https://www.seas.ucla.edu/~vandenbe/publications/chordalsdp.pdf) provides an excellent overview of this area.
- The clique graph implementation follows [Garstka et al.](https://arxiv.org/abs/1911.05615). The accompanying solver [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl) has chordal decomposition built-in and several great customization features.
