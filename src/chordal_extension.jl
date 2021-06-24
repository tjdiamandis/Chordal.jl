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


# Gets chordal extension: a reordering + associated graph from cholesky
function get_chordal_extension(sp_pattern::SparseMatrixCSC; perm=nothing, verbose=false)
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
