var documenterSearchIndex = {"docs":
[{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"The source files for all examples can be found in /examples.","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"EditURL = \"https://github.com/tjdiamandis/Chordal.jl/blob/master/examples/sdp_lmi.jl\"","category":"page"},{"location":"examples/sdp_lmi/#Semidefinite-Program-Decomposition","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"","category":"section"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"This example illustrates how to decompose a positive semidefinite constraint with Chordal.jl. The decomposed constraint can then be plugged into your favorite semidefinite program (SDP) solver.","category":"page"},{"location":"examples/sdp_lmi/#Problem-Setup","page":"Semidefinite Program Decomposition","title":"Problem Setup","text":"","category":"section"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"We first generate a random SDP using generate_random_sdp(n), which constructs problems of the form","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"beginarrayll\ntextminimize c^Tx \ntextsubject to  sum_i=1^m F_ix_i + G succeq 0\nendarray","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"using Chordal\nusing JuMP, SCS\nusing SparseArrays, LinearAlgebra, Random\n\nrand_seed = 1234\nn = 500\noptimizer = SCS.Optimizer\n\n# D is the dual variable and xstar is an optimal solution\nc, F, G, xstar, D = Chordal.generate_random_sdp(n, rand_seed=rand_seed);\nnothing #hide","category":"page"},{"location":"examples/sdp_lmi/#Standard-Solve","page":"Semidefinite Program Decomposition","title":"Standard Solve","text":"","category":"section"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"First, we build the SDP using JuMP and solve it using SCS without chordal decomposition.","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"# Without decomposition\nm = Model(optimizer)\nJuMP.set_silent(m)\n@variable(m, x[1:n])\n@constraint(m, sum(F[i]*x[i] for i in 1:n) + G ∈ PSDCone())\n@objective(m, Min, dot(c, x))\ntime_non_chordal = @elapsed optimize!(m)\n\nxv = value.(x)\npstar = dot(c,xv)\n@info \"termination status: $(termination_status(m))\"\n@info \"solution status: $(primal_status(m))\"\n@info \"Difference with true optimal: $(norm(xv - xstar))\"\n@info \"Optimal value: $pstar\"","category":"page"},{"location":"examples/sdp_lmi/#Solving-after-Chordal-Decomposition","page":"Semidefinite Program Decomposition","title":"Solving after Chordal Decomposition","text":"","category":"section"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"Next, we use chordal decomposition tools from Chordal.jl to decompose the PSD constraint into smaller blocks.","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"# Chordal Decomposition Setup\nm2 = Model(optimizer)\nJuMP.set_silent(m2)\n@variable(m2, y[1:n])\n\n# build_A is a helper function to construct the SDP\n# A = sum(F[i]*y[i] for i in 1:n) + G\n# A + S ∈ PSDCone(), where S ⪰ 0\nA = Chordal.build_A(y, F, G)\nnnzA = nnz(A)\n@info \"There are $nnzA nonzeros; density = $(round(nnzA/n^2, digits=3))\"\ntime_constraints = @elapsed build_constraints_lmi!(m2, A)\n# NOTE: Can also call build_constraints_lmi!(m2, y, F, G)\n@info \"Finished adding constraints, time = $(round(time_constraints, digits=3))s\"\n@objective(m2, Min, dot(c,y))\n@info \"Finished setup\"","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"We can solve the decomposed model with any SDP solver (here we use SCS).","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"# Chordal Decomposition Solve\ntime_chordal = @elapsed optimize!(m2)\nyv = value.(y)\npstar_chordal = dot(c,yv)\n\n@info \"termination status: $(termination_status(m2))\"\n@info \"solution status: $(primal_status(m2))\"\n@info \"difference with non-chordal: $(norm(yv - xv))\"\n@info \"difference with true variable: $(norm(yv - xstar))\"\n@info \"Time with chordal is $(round(time_chordal, digits=3))s vs $(round(time_non_chordal, digits=3))s\"","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"","category":"page"},{"location":"examples/sdp_lmi/","page":"Semidefinite Program Decomposition","title":"Semidefinite Program Decomposition","text":"This page was generated using Literate.jl.","category":"page"},{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"","category":"page"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [Chordal]","category":"page"},{"location":"api/#Chordal.build_A-Union{Tuple{S}, Tuple{T}, Tuple{Vector{JuMP.VariableRef}, AbstractArray{SparseArrays.SparseMatrixCSC{T, S}, 1}, SparseArrays.SparseMatrixCSC{T, S}}} where {T<:Number, S<:Integer}","page":"API Reference","title":"Chordal.build_A","text":"build_A(y::Vector{JuMP.VariableRef}, F::AbstractVector{SparseMatrixCSC{T, S}}, G::SparseMatrixCSC{T, S}) where {T <: Number, S <: Integer}\n\nBuilds a JuMP.GenericAffExpr A = F_1 y_1 + F_2 y_2 + ... + F_n y_n + G.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.build_constraints_lmi!-Tuple{JuMP.Model, SparseArrays.SparseMatrixCSC{JuMP.AffExpr, Int64}}","page":"API Reference","title":"Chordal.build_constraints_lmi!","text":"build_constraints_lmi!(m::JuMP.Model, A::SparseMatrixCSC{JuMP.AffExpr, Int}; verbose=false)\n\nAdds the constraint A + S ∈ PSDCone(), where S ⪰ 0 to JuMP model m\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.build_constraints_lmi!-Union{Tuple{S}, Tuple{T}, Tuple{JuMP.Model, Vector{JuMP.VariableRef}, AbstractArray{SparseArrays.SparseMatrixCSC{T, S}, 1}, SparseArrays.SparseMatrixCSC{T, S}}} where {T<:Number, S<:Integer}","page":"API Reference","title":"Chordal.build_constraints_lmi!","text":"build_constraints_lmi!(m::JuMP.Model, A::SparseMatrixCSC{JuMP.AffExpr, Int}; verbose=false)\n\nAdds the constraint F_1 y_1 + F_2 y_2 + ... + F_n y_n + G ∈ PSDCone() to JuMP model m\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.build_perm_matrix-Tuple{Any}","page":"API Reference","title":"Chordal.build_perm_matrix","text":"build_perm_matrix(p)\n\nBuilds a sparse permutation matrix form permutation p such that P*x == x[p]\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.edm_completion-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.edm_completion","text":"edm_completion(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}\n\nReturns a Euclidean distance matrix completion of A if it is EDM-completable using Algorithm 11.1 in [VA15].\n\nReference\n\n[VA15] Chordal Graphs and Semidefinite Optimization by Lieven Vandenberghe and Martin Andersen\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.generate_clique_graph-Tuple{Any, Integer}","page":"API Reference","title":"Chordal.generate_clique_graph","text":"generate_clique_graph(cliques, n::Integer)\n\nGenerates CliqueGraph datastructure for an undirected graph from cliques and the number of nodes n.\n\nReference\n\nA clique graph based merging strategy for decomposable SDPs by Michael Garstka, Mark Cannon, Paul Goulart\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.generate_clique_tree-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.generate_clique_tree","text":"generate_clique_tree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}\n\nGenerates a CliqueTree from chordally sparse matrix A.\n\nReference\n\n[VA15] Chordal Graphs and Semidefinite Optimization by Lieven Vandenberghe and Martin Andersen\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.generate_random_sdp-Tuple{Any}","page":"API Reference","title":"Chordal.generate_random_sdp","text":"generate_random_sdp(n; rand_seed=0)\n\nGenerates a random dual form SDP with side dimension n: min c'*x s.t. sum(F[i]*x[i]) + G ⪰ 0\n\nReturns c, F, G, xstar, D, where xstar and D are optimal primal and dual variables respectively\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_children_from_par-Tuple{Vector{Int64}}","page":"API Reference","title":"Chordal.get_children_from_par","text":"get_children_from_par(par)\n\nCompute the children of each vertex v in a tree specified by par\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_chordal_extension-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.get_chordal_extension","text":"get_chordal_extension(sp_pattern::SparseMatrixCSC; perm=\"amd\", verbose=false)\n\nReturns p, ip, L, where L is a lower-triangular matrix representing the chordal extension of the undirected graph sp_pattern after fill-reducing permutation p (with inverse permutation ip).\n\nOptions for perm include \t- \"nothing\" (no permutation) \t- amd (approximate minimum degree) \t- metis (uses METIS) \t- cuthill-mckee (Cuthill-McKee bandwidth reducing algorithm)\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_cliques-Tuple{Chordal.CliqueGraph}","page":"API Reference","title":"Chordal.get_cliques","text":"get_cliques(cg::CliqueGraph)\n\nReturns a list of the cliques in CliqueGraph cg.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_cliques-Tuple{Chordal.CliqueTree}","page":"API Reference","title":"Chordal.get_cliques","text":"get_cliques(ct::CliqueTree)\n\nReturns the cliques in the graph represented by CliqueTree ct.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_cliques-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.get_cliques","text":"get_cliques(L::SparseMatrixCSC)\n\nReturns the (maximal) cliques of the undirected graph represented by lower triangular matrix L.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_etree-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.get_etree","text":"get_etree(A)\n\nCompute the parent function of the elimination tree of a symmetric sparse matrix A\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_higher_deg-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.get_higher_deg","text":"get_higher_deg(L::SparseMatrixCSC)\n\nReturns the higher degrees of each node in the graph represented by lower triangular matrix L.\n\nNOTE: the algorithm assums that L has order σ = 1:n and zeros on the diagonal.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_postordering-Tuple{Any, Any}","page":"API Reference","title":"Chordal.get_postordering","text":"get_postordering(par, child)\n\nGets a postordering of the tree represented by vector of parents par and vector of children child.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_selectors-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.get_selectors","text":"get_selectors(input_mat::SparseMatrixCSC; verbose=true, ret_cliques=true)\n\nReturns the (merged) cliques of the graph corresponding to the sparsity pattern of input_mat (after a permutation to reduce fill-in) and optionally the cliques. Also returns the fill-reducing permutation and and inverse permutation used.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.get_sparsity_pattern_from_cliques-Tuple{Any}","page":"API Reference","title":"Chordal.get_sparsity_pattern_from_cliques","text":"get_sparsity_pattern_from_cliques(cliques)\n\nReturns a sparse matrix sp corresponding to the graph with maximal cliques cliques\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.is_chordal-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.is_chordal","text":"is_chordal(A)\n\nTests if the graph represented by symmetric sparse matrix A is chordal. This function can also be used with AbstractGraph objects from LightGraphs.jl.\n\nReferences\n\nSimple linear-time algorithms to test chordality of graphs, test acyclicity of hypergraphs, and selectively reduce acyclic hypergraphs by Robert Tarjan and Mihalis Yannakakis\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.is_edm-Tuple{AbstractMatrix{T} where T}","page":"API Reference","title":"Chordal.is_edm","text":"is_edm(M::AbstractMatrix)\n\nChecks if M is a Euclidean Distance Matrix. Note: a matrix is a EDM iff it is negative semidefinite on the subspace on the subspace orthogonal to the all ones vector.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.is_edm_completable-Tuple{SparseArrays.SparseMatrixCSC, Chordal.CliqueTree}","page":"API Reference","title":"Chordal.is_edm_completable","text":"is_edm_completable(A::SparseMatrixCSC, ct::CliqueTree)\n\nChecks if A is EDM completable. Requires the associated CliqueTree ct.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.is_separable-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.is_separable","text":"is_separable(sp::SparseMatrixCSC)\n\nReturns true if the sparsity pattern is separable (i.e., block diagonal).\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.make_selectors_from_clique_graph-Tuple{Chordal.CliqueGraph, Any}","page":"API Reference","title":"Chordal.make_selectors_from_clique_graph","text":"make_selectors_from_clique_graph(cg::CliqueGraph, n)\n\nBuilds selector matrices Tℓ from cg, a CliqueGraph with n vertices.\n\nTℓ*X*Tℓ' is the submatrix of X corresponding to clique ℓ.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.make_selectors_from_cliques-Tuple{Any, Any}","page":"API Reference","title":"Chordal.make_selectors_from_cliques","text":"make_selectors_from_cliques(cliques, n)\n\nBuilds selector matrices Tℓ from cliques, a list of the cliques in a graph with n vertices.\n\nTℓ*X*Tℓ' is the submatrix of X corresponding to clique ℓ.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.max_supernode_etree-Tuple{SparseArrays.SparseMatrixCSC, Vector{Int64}}","page":"API Reference","title":"Chordal.max_supernode_etree","text":"max_supernode_etree(L::SparseMatrixCSC, etree_par::Vector{Int})\n\nConstructs the supernodal elimination tree of the chordal graph represented by lower triangular matrix L and with elimination tree etree_par. Returns the reprsentation vertices, the supernodal elimination tree, and a supernode memership indicator.\n\nImplemented using [VA15, Algorithm 4.1], which was originally formulated by [PS89].\n\nReferences\n\n[PS89] Compact clique tree data structures in sparse matrix factorizations by Alex Pothen and Chunguang Sun\n[VA15] Chordal Graphs and Semidefinite Optimization by Lieven Vandenberghe and Martin Andersen\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.maxdet_completion-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.maxdet_completion","text":"maxdet_completion(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}\n\nReturns the maximum determinant completion of chordal sparse matrix A using Algorithm 10.2 in [VA15].\n\nReference\n\n[VA15] Chordal Graphs and Semidefinite Optimization by Lieven Vandenberghe and Martin Andersen\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.maxdet_completion_etree-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.maxdet_completion_etree","text":"maxdet_completion_etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}\n\nReturns the cholesky factors of the inverse of the maximum determinant completion of chordal sparse matrix A using Algorithm 4.2 in [ADV14]. This algorithm uses the elimination tree of A and, therefore, BLAS level 2 operations.\n\nReference\n\n[ADV12] Logarithmic barriers for sparse matrix cones by Martin S. Andersen, Joachim Dahl, Lieven Vandenberghe\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.maxdet_completion_factors-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.maxdet_completion_factors","text":"maxdet_completion_etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}\n\nReturns the cholesky factors of the inverse of the maximum determinant completion of chordal sparse matrix A using Algorithm 7.3 in [ADV14]. This algorithm uses the supernodal elimination tree of A and, therefore, BLAS level 3 operations.\n\nReference\n\n[ADV12] Logarithmic barriers for sparse matrix cones by Martin S. Andersen, Joachim Dahl, Lieven Vandenberghe\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.maximum_cardinality_search-Tuple{SparseArrays.SparseMatrixCSC}","page":"API Reference","title":"Chordal.maximum_cardinality_search","text":"maximum_cardinality_search(A)\n\nCompute the perfect elimination ordering for a chordal graph represented by a symmetric sparse matrix A. This function can also be used with AbstractGraph objects from LightGraphs.jl.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.merge_cliques!-Tuple{Chordal.CliqueGraph}","page":"API Reference","title":"Chordal.merge_cliques!","text":"merge_cliques!(cg::CliqueGraph; verbose=false)\n\nMerges cliques in cg in a greedy fashion starting with the pair with the largest weight_function. Stops when all pairs of cliques have a non-positive weight.\n\nReference\n\nA clique graph based merging strategy for decomposable SDPs by Michael Garstka, Mark Cannon, Paul Goulart\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.minrank_completion-Union{Tuple{SparseArrays.SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}} where {Tv<:AbstractFloat, Ti<:Integer}","page":"API Reference","title":"Chordal.minrank_completion","text":"minrank_completion(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}\n\nReturns the cholesky factor of the minimum rank completion (A_complete = Y*Yᵀ) of chordal sparse matrix A using Algorithm 2 in [Sun15]. This algorithm uses the supernodal elimination tree of A and, therefore, BLAS level 3 operations.\n\nReference\n\n[Sun15] Decomposition Methods for Semidefinite Optimization (PhD Thesis) by Yifan Sun\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.order_snds!-Tuple{Chordal.CliqueTree}","page":"API Reference","title":"Chordal.order_snds!","text":"order_snds!(ct::CliqueTree)\n\nOrders supernodes such that elements of each supernode are numbered consecutively and the order is a topological ordering of the reprsentative vertices in the supernodal elimination tree. (See [VA15, 4.6])\n\nReference\n\n[VA15] Chordal Graphs and Semidefinite Optimization by Lieven Vandenberghe and Martin Andersen\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.preprocess!-Union{Tuple{SparseArrays.SparseMatrixCSC{T, var\"#s2\"} where var\"#s2\"<:Integer}, Tuple{T}} where T","page":"API Reference","title":"Chordal.preprocess!","text":"preprocess!(mat::SparseMatrixCSC{T, <:Integer}) where {T}\n\nChecks that mat is symmetric and drops numerical zeros.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.reconstruct_from_sparse_varref-Tuple{Any, Any}","page":"API Reference","title":"Chordal.reconstruct_from_sparse_varref","text":"reconstruct_from_sparse_varref(Zref, n)\n\nGets the value of JuMP.Containers.SparseAxisArray Zref and returns the result as a SparseMatrixCSC\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.sparsity_pattern-Union{Tuple{AbstractArray{SparseArrays.SparseMatrixCSC{T, S}, 1}}, Tuple{S}, Tuple{T}} where {T, S<:Integer}","page":"API Reference","title":"Chordal.sparsity_pattern","text":"sparsity_pattern(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T, S <: Integer}\n\nReturns the 0/1 aggregate sparsity pattern of the matrices in mats as a sparse matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.unzip-Tuple{Any}","page":"API Reference","title":"Chordal.unzip","text":"unzip(a)\n\nUnzips a list of tuples a.\n\nExample\n\njulia> unzip([(1,2), (3,4), (5,6)])\n([1, 3, 5], [2, 4, 6])\n\n\n\n\n\n","category":"method"},{"location":"api/#Chordal.weight_function-Tuple{Any, Any}","page":"API Reference","title":"Chordal.weight_function","text":"weight_function(ci, cj)\n\nDefines the weight function used to determine if cliques ci and cj should be merged. Cliques are merged if this is positive.\n\nReference\n\nA clique graph based merging strategy for decomposable SDPs by Michael Garstka, Mark Cannon, Paul Goulart\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Chordal","category":"page"},{"location":"#Chordal","page":"Home","title":"Chordal","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Chordal.jl is a package for working with sparse matrices that have a chordal sparsity pattern. Check out the examples.","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"The source files for all examples can be found in /examples.","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"EditURL = \"https://github.com/tjdiamandis/Chordal.jl/blob/master/examples/edm.jl\"","category":"page"},{"location":"examples/edm/#Euclidean-Distance-Matrix-Completion","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"","category":"section"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"This example illustrates how to use Chordal.jl to complete a partially-specified Euclidean Distance Matrix (or determine no such completion exists).","category":"page"},{"location":"examples/edm/#Euclidean-Distance-Matrices","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrices","text":"","category":"section"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"Given N vectors x_1dotsx_n in mathbbR^n, the associated Euclidean Distance Matrix (EDM) D in mathbbR^N times N records the squared distance between each pair x_i x_j. Specifically,","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"D_ij = x_i - x_j_2^2","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"Given a partially-specified EDM D which has a chordal sparsity pattern, we would like to find the missing entries.","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"using Chordal\nusing LinearAlgebra, SparseArrays","category":"page"},{"location":"examples/edm/#Small-Example","page":"Euclidean Distance Matrix Completion","title":"Small Example","text":"","category":"section"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"This small example is taken from Jon Dattorro's Convex Optimization & Euclidean Distance Geometry eq (1041). First, we construct the partially-specified EDM.","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"D = sparse([\n    0.0 1 5 0;\n    1 0 4 1;\n    5 4 0 1;\n    0 1 1 0\n])\n@show D\n\n# Next, we complete the matrix\nD_complete = sparse(edm_completion(D))\nprintln(\"\\nCompletion:\")\n@show D_complete","category":"page"},{"location":"examples/edm/#Larger-Random-Example","page":"Euclidean Distance Matrix Completion","title":"Larger Random Example","text":"","category":"section"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"First, we generate an EDM from 17 vectors in mathbbR^3.","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"n = 17\nr = 3\nV = randn(n, r)\nVV = V*V'\nD_full = diag(VV)*ones(n)' + ones(n)*diag(VV)' - 2VV;\nnothing #hide","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"Next, we remove entries. and complete the The sparsity pattern is from Figure 4.2 in Lieven Vandenberghe and Martin S. Andersen's Chordal Graphs and Semidefinite Optimization.","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"ijs = [(3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,3), (5,4), (8,7),\n       (9,5), (9,6), (9,7), (9,8), (11,10), (13,10), (13,11), (13,12),\n       (14,10), (14,11), (14,12), (14,13),\n       (15, 1), (15,3), (15,4), (15,5), (15,7), (15,8), (15,9),\n       (16,5), (16,6), (16,9), (16,12), (16,13), (16,14), (16,15),\n       (17,10), (17,11), (17,12), (17,13), (17,14), (17,15), (17,16)]\nappend!(ijs, [(i,i) for i in 1:n])\nII, JJ = Chordal.unzip(ijs)\nsp = sparse(II, JJ, ones(length(II)))\nsp = sp + tril(sp)'\n\n# Remove entries s.t. the remaining entries have a chordal sparsity pattern\nD = sp .* D_full\n\n# Complete the matrix\n@show Chordal.is_edm(edm_completion(D))","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"","category":"page"},{"location":"examples/edm/","page":"Euclidean Distance Matrix Completion","title":"Euclidean Distance Matrix Completion","text":"This page was generated using Literate.jl.","category":"page"}]
}
