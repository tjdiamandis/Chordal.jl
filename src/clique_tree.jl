"""
    etree(A)
Compute the elimination tree of a symmetric sparse matrix `A` from `triu(A)`.
"""
function etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && throw(ArgumentError("Matrix must be symmetric"))
    m, n = size(A)
    par = zeros(Ti, n)
    ancestor = zeros(Ti, n)

    _etree(A, par, ancestor)
    return par
end

# (essentially) nonallocating verion of above
function _etree(A::SparseMatrixCSC{Tv, Ti}, par::Vector{Ti}, ancestor::Vector{Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    for col::Ti in 1:n, p in nzrange(A, col) #A.colptr[col]:(A.colptr[col+1]-1)
        i = rowvals(A)[p]
        while !iszero(i) && i < col
            i_next = ancestor[i]
            ancestor[i] = col
            if iszero(i_next)
                par[i] = col
            end
            i = i_next
        end
    end
    return nothing
end
