import QDLDL
using LightGraphs: SimpleGraph, maximal_cliques
using CuthillMcKee: symrcm
import AMD, Metis


function preprocess!(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T <: Number, S <: Integer}
    for mat in mats
        preprocess!(mat)
    end
end


"""
	preprocess!(mat::SparseMatrixCSC{T, <:Integer}) where {T}

Checks that `mat` is symmetric and drops numerical zeros.
"""
function preprocess!(mat::SparseMatrixCSC{T, <:Integer}) where {T}
    n = size(mat, 1)
    @assert size(mat) == (n, n)
	for i in 1:n, j in i+1:n
		@assert mat[i,j] .== mat[j,i]
	end
    dropzeros!(mat)
end


# Returns 0/1 sparsity pattern of a set of matrices
"""
	sparsity_pattern(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T, S <: Integer}

Returns the 0/1 aggregate sparsity pattern of the matrices in `mats` as a sparse matrix.
"""
function sparsity_pattern(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T, S <: Integer}
	n = size(mats[1], 1)
	for mat in mats
		!(n == size(mat, 1) && issymmetric(mat)) && error(ArgumentError(
			"Matrices must be symmetric and have the same dimensions."
		))
	end

	# Build 0/1 sparsity pattern matrix
	ret = spzeros(n,n)
	for mat in mats
	    @inbounds for col in 1:n, k in nzrange(mat, col)
	        ret[rowvals(mat)[k], col] = 1.0
		end
	end

	# Always assume diagonal can be nonzero
	ret[diagind(ret)] .= 1

    return ret
end
sparsity_pattern(mat::SparseMatrixCSC) = sparsity_pattern([mat])


# Gets chordal extension: a reordering + associated graph from cholesky
"""
	get_chordal_extension(sp_pattern::SparseMatrixCSC; perm="amd", verbose=false)

Returns `p, ip, L`, where `L` is a lower-triangular matrix representing the
chordal extension of the undirected graph `sp_pattern` after fill-reducing permutation
`p` (with inverse permutation `ip`).
"""
function get_chordal_extension(sp_pattern::SparseMatrixCSC; perm="amd", cost=nothing, verbose=false)
	!issymmetric(sp_pattern) && error(ArgumentError("Matrix must be symmetric"))
	n = size(sp_pattern)[1]

	# Permutations: nothing, amd, vertex separators, TODO: add more options?
	# See survey Minimal triangulations of graphs: A survey
	# https://www.sciencedirect.com/science/article/pii/S0012365X05006060
	if perm == "none" || perm == "nothing"
		p = nothing
	elseif perm == "amd"
		p = AMD.amd(sp_pattern)
	elseif perm == "vsep"
		# Convert to Int64 for QDLDL (ret as Int32)
		# first arg = perm (second = iperm)
		p = Int.(Metis.permutation(sp_pattern)[1])
    elseif perm == "cuthill-mckee"
		p = cuthill_mckee(sp_pattern)
	elseif perm == "greedy" && !isnothing(cost)
		# Choose v according to criterion X
		# Let v be the eith vertex in ordering p
		# Let H be the graph obstained by eliminating v from H
		#  - make neighborhood of v a clique and then remove v
		# Return H, p
		# Criteria:
		# Greedy degree: 				v = argminᵤ d(u)
		# Greedy fill in				argminᵤ f(u) = num pairs of non-adj neighbors
		# Greedy degree + fill in		d(u) + f(u)
		# Greedy sparsest subgraph		f(u) - d(u)
		# Greedy fill in degree			d(u) + 1/n^2 * f(u)
		# Greedy degree fill in			f(u) + 1/n * d(u)
		n = size(sp_pattern, 1)
		p = Vector{Int}(undef, n)
		active_inds = Vector(1:n)
		for i in 1:n
			H = @view(sp_pattern[active_inds, active_inds])
			min_cost_vert = get_min_cost_vert(H, cost)
			p[i] = active_inds[min_cost_vert]
			deleteat!(active_inds, [min_cost_vert])
		end

	else
		error(ArgumentError("Invalid permutation"))
	end
	# TODO: add greedy
	# Algorithm: (summarized in Bodlaender & Koster, Treewith computations I. Upper bounds)
	# H = copy(G)
	# for i in 1:n

	F = QDLDL.qdldl(sp_pattern, perm=p, logical=true)

	if verbose
		num_nonzero = 2*nnz(F.L) + n
		num_nonzero_added = (num_nonzero - nnz(sp_pattern)) ÷ 2

		@info "Chordal Extension added $num_nonzero_added nonzeros."
		@info "Density: $(round(num_nonzero/n^2; digits=5))"
	end

    return F.perm, F.iperm, F.L
end


"""
	get_cliques(L::SparseMatrixCSC)

Returns the (maximal) cliques of the undirected graph represented by lower
triangular matrix `L`.
"""
function get_cliques(L::SparseMatrixCSC)
	G = SimpleGraph(L + L')
	return maximal_cliques(G)
end


function cuthill_mckee(sp::SparseMatrixCSC)
	# TODO: Custom implementation for efficiency?
	# Cuthill-McKee (https://en.wikipedia.org/wiki/Cuthill–McKee_algorithm)
	return symrcm(sp)
end


function get_min_cost_vert(H, cost::Function)
	# NOTE: Removed check for performance
	# !issymmetric(H) && error("adj matrix must be symmetric")
	n = size(H, 1)

	min_cost = Inf
	amin_cost = 0
	for i in 1:n
		v_cost = cost(H, i)
		if v_cost < min_cost
			min_cost = v_cost
			amin_cost = i
		end
	end

	return amin_cost
end


# Some cost functions
# H is a (sparse) adj matrix of a graph (CSC format)
# i is the vertex to eval cost on
# https://www.sciencedirect.com/science/article/pii/S0890540109000947
function cost_deg(H::SparseMatrixCSC, i)
	# NOTE: assumes 0/1 matrix
	return length(H.colptr[i]:(H.colptr[i+1]-1))
end

function cost_fill_in(H, i)
	nna_neighbors = 0
	display(fieldnames(H))
	neighbors = @view(H.rowval[H.colptr[i]:(H.colptr[i+1]-1)])
	for i in neighbors, j in neighbors
		i == j && continue
		if H[i,j] == 0
			nna_neighbors += 1
		end
	end
	return nna_neighbors
end

cost_deg_plus_fill_in(H, i) = cost_deg(H, i) + cost_fill_in(H, i)
cost_sparsest_subset(H, i) 	= cost_fill_in(H, i) - cost_deg(H, i)
cost_fill_in_degree(H, i) 	= cost_deg(H, i) + 1/size(H,1)^2 * cost_fill_in(H, i)
cost_deg_fill_in(H, i) 		= 1/size(H,1) * cost_deg(H, i) + cost_fill_in(H, i)
