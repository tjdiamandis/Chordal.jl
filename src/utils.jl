function build_perm_matrix(p)
    n = length(p)
    P = spzeros(n, n)
    for i in 1:n
        P[i, p[i]] = 1.0
    end
    return P
end
