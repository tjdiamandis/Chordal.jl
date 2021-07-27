

function edm_completion(A::SparseMatrixCSC{Tv, Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    !issymmetric(A) && error(ArgumentError("A must be symmetric"))
    n = size(A, 1)

    sp = sparsity_pattern(A)
    !is_chordal(A; peo=1:n) && error(ArgumentError("A must have a chordal sparsity pattern"))
    perm, iperm, L = get_chordal_extension(sp; perm=nothing, verbose=false)
    ct = CliqueTree(L)
    order_snds!(ct)
    iperm = invperm(ct.perm)

    !is_edm_completable(A[ct.perm, ct.perm], ct) && error("A is not EDM completable")

    W = Matrix(A[ct.perm, ct.perm])

    for j in (length(ct.snds)-1):-1:1
        vrep_ind = ct.postordering[j]
        vrep = ct.vreps[vrep_ind]

        ν = ct.snds[vrep_ind]
        α = ct.seps[vrep_ind]
        η = filter(x->(!(x in ν) && !(x in α)), vrep+1:n)

        aj = first(α)

        # See eq (11.4)
        # (I - 1*eₖ') * A * (I - eₖ*1') = A - A[k,k]*1*1' - 1*ãₖ' - aₖ*1'
        Yνν = -0.5*W[ν, ν] .- 0.5*W[aj, aj]
        Yνν .+= 0.5 .* W[ν, aj]
        Yνν .+= 0.5 .* W[aj, ν]'

        Yαα = -0.5*W[α, α] .- 0.5*W[aj, aj]
        Yαα .+= 0.5 .* W[α, aj]
        Yαα .+= 0.5 .* W[aj, α]'

        Yηη = -0.5*W[η, η] .- 0.5*W[aj, aj]
        Yηη .+= 0.5 .* W[η, aj]
        Yηη .+= 0.5 .* W[aj, η]'

        Yηα = -0.5*W[η, α] .- 0.5*W[aj, aj]
        @. Yηα += 0.5 * W[η, aj]
        @. Yηα += 0.5 * W[aj, α]'

        Yαν = -0.5*W[α, ν] .- 0.5*W[aj, aj]
        @. Yαν += 0.5 * W[α, aj]
        @. Yαν += 0.5 * W[aj, ν]'


        d, V = LAPACK.syev!('V', 'L', Yαα)
        λmax = maximum(d)
        Yαα_pinv = V * Diagonal(map(x->(x > λmax*eps() ? 1.0/x : 0.0), d)) * V'
        cache = Yαα_pinv * Yαν
        @views mul!(W[η, ν], Yηα, cache)
        W[η, ν] .*= -2
        W[η, ν] .+= diag(Yηη)
        W[η, ν] .+= diag(Yνν)'

        # W[η, ν] = -2*Yηα*Yαα_pinv*Yαν + ones(length(η))*diag(Yνν)' + diag(Yηη)*ones(length(ν))'
        # BLAS.ger!(1.0, ones(length(η)), diag(Yνν), W[η, ν])
        # BLAS.ger!(1.0, diag(Yηη), ones(length(ν)), W[η, ν])

        # μ = sparsevec([aj], [1.0], n)
        # tmp  = -0.5*(I - ones(n)*μ')*W*(I - μ*ones(n)')
        # rtol = sqrt(eps(real(float(one(eltype(tmp[α, α]))))))
        # W[η, ν] .= -2*tmp[η, α]*pinv(tmp[α, α]; rtol=rtol)*tmp[α, ν] +
        #           ones(length(η))*diag(tmp[ν,ν])' + diag(tmp[η,η])*ones(length(ν))'


        W[ν, η] .= W[η, ν]'
        # display(Matrix(W[iperm, iperm] .* .!(sp .> 0)))
        # break
    end

    return W[iperm, iperm]

end


function is_edm_completable(A::SparseMatrixCSC, ct::CliqueTree)
    cliques = get_cliques(ct)
    for c in cliques
        !is_edm(A[c,c]) && return false
    end
    return true
end


function is_edm(M::AbstractMatrix)
    (!iszero(diag(M)) || !issymmetric(M)) && return false
    n = size(M, 1)
    Q = spdiagm(n, n-1, 0 => ones(n-1), -1 => -ones(n-1))

    return isposdef(-Symmetric(Q'*M*Q) + 1e-10*I)
end
