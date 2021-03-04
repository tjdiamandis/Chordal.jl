import QDLDL
using LightGraphs: SimpleGraph, maximal_cliques

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


# Gets sparsity pattern of the sum of mats
function sparsity_pattern(mats::AbstractVector{SparseMatrixCSC{T,S}}) where {T <: Number, S <: Integer}
	n = size(mats[1])[1]
    # TODO: For very large matrices, can just record the CartesianIndex's or I, J's
    pattern = falses(n, n)
    for mat in mats
        @inbounds for col = 1:n, k = mat.colptr[col]:(mat.colptr[col+1]-1)
            pattern[rowvals(mat)[k], col] = true
        end
    end

    # TODO: should this return as CartesianIndex?
	ret = spzeros(n,n)
	for c in findall(pattern)
		i, j = c[1], c[2]
		ret[i,j] = ret[j,i] = 1.0
	end
    return ret
end


function sparsity_pattern(mat::SparseMatrixCSC{T, <:Integer}) where {T}
    n = size(mat, 1)
	ret = spzeros(n,n)

    @inbounds for col = 1:n, k = mat.colptr[col]:(mat.colptr[col+1]-1)
        ret[rowvals(mat)[k], col] = 1.0
	end

    return ret
end


# TODO: Add ref for attribution to COSMO
# NOTE: this assumes a sparse lower triangular matrix L
function _connect_graph!(L::SparseMatrixCSC; verbose=true)
	# unconnected blocks don't have any entries below the diagonal in their right-most column
	count = 0
	m = size(L, 1)
	row_val = L.rowval
	col_ptr = L.colptr
	for j = 1:m-1
		connected = false
		for k in col_ptr[j]:col_ptr[j+1]-1
			if row_val[k] > j
				connected  = true
				break
			end
		end
		if !connected
			L[j+1, j] = 1
			count += 1
		end
	end
	verbose && @info "Connected graph; added $count edges"
end


# Gets chordal extension: a reordering + associated graph from cholesky
function get_chordal_extension(sp_pattern::SparseMatrixCSC; verbose=false)
	# TODO: Min degree preordering for sp_pattern? Or just let Cholesky handle?
    F = QDLDL.qdldl(sp_pattern, logical=true)
	_connect_graph!(F.L; verbose=verbose)
	num_nonzero_added = 2*nnz(F.L) + size(sp_pattern)[1] - nnz(sp_pattern)
	verbose && @info "Chordal Extension added $num_nonzero_added nonzeros."
    return F.perm, F.L
end


function get_cliques(L::SparseMatrixCSC)
	G = SimpleGraph(L + L')
	return maximal_cliques(G)
end
