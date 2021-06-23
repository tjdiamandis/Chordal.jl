
function etree(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    m, n = size(A)
    m != n && throw(DimensionMismatch("Matrix must be square"))

    parent = zeros(Ti, n)
    ancestor = zeros(Ti, n)

    for col in 1:n, p in nzrange(A, col)
        i = rowvals(A)[p]
        while i != 0 && i < col
            i_next = ancestor[i]
            ancestor[i] = col
            if i_next == 0
                parent[i] = col
            end
            i = i_next
        end
    end

    return parent
end
