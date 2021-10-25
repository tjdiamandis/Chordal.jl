@testset "utils" begin

    # Test permutation matrix
    p = [5, 3, 4, 1, 2, 6]
    x = randn(length(p))
    @test x[p] == Chordal.build_perm_matrix(p)*x

    # Test random SDP
    n = 100
    c, F, G, xstar, D = Chordal.generate_random_sdp(n)
    Fx = sum(F[i]*xstar[i] for i in 1:n) + G
    @test (
        all(eigvals(Matrix(Fx)) .>= -sqrt(eps())) &&
        all(eigvals(D) .>= -sqrt(eps())) &&
        all([0 ≈ c[i] - tr(F[i]*D) for i in 1:n]) &&
        tr(Fx*D) ≈ 0
    )
end