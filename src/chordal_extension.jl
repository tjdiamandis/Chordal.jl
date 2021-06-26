import QDLDL
using LightGraphs: SimpleGraph, maximal_cliques
import AMD, Metis

function preprocess!(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T <: Number, S <: Integer}
    for mat in mats
        preprocess!(mat)
    end
end

function preprocess!(mat::SparseMatrixCSC{T, <:Integer}) where {T}
    n = size(mat, 1)
    @assert size(mat) == (n, n)
	for i in 1:n, j in i+1:n
		@assert mat[i,j] .== mat[j,i]
	end
    dropzeros!(mat)
end


# Returns 0/1 sparsity pattern of a set of matrices
function sparsity_pattern(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T <: Number, S <: Integer}
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

    return ret
end
sparsity_pattern(mat::SparseMatrixCSC) = sparsity_pattern([mat])


# Gets chordal extension: a reordering + associated graph from cholesky
function get_chordal_extension(sp_pattern::SparseMatrixCSC; perm="amd", verbose=false)
	!issymmetric(sp_pattern) && error(ArgumentError("Matrix must be symmetric"))
	n = size(sp_pattern)[1]

	# Permutations: nothing, amd, vertex separators, TODO: add more options?
	if isnothing(perm)
		p = nothing
	elseif perm == "amd"
		p = AMD.amd(sp_pattern)
	elseif perm == "vsep"
		# Convert to Int64 for QDLDL (ret as Int32)
		# first arg = perm (second = iperm)
		p = Int.(Metis.permutation(sp_pattern)[1])
	else
		error(ArgumentError("Invalid permutation"))
	end

	F = QDLDL.qdldl(sp_pattern, perm=p, logical=true)
	num_nonzero = 2*nnz(F.L) + n
	num_nonzero_added = (num_nonzero - nnz(sp_pattern)) ÷ 2

	verbose && @info "Chordal Extension added $num_nonzero_added nonzeros."
	verbose && @info "Density: $(round(num_nonzero/n^2; digits=5))"
    return F.perm, F.iperm, F.L
end


function get_cliques(L::SparseMatrixCSC)
	G = SimpleGraph(L + L')
	return maximal_cliques(G)
end
