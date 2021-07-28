using LightGraphs: SimpleGraph

@testset "Chordal Graph tests" begin
    # From Garstka paper
    n1 = 9
    nonzero_inds = vcat([
            (3,1), (3,2), (5,4), (6,1), (6,3), (7,3), (7,6), (8,3), (8,4),
            (8,5), (8,6), (8,7), (9,6), (9,7), (9,8)
        ],
        [(i,i) for i in 1:n1]
    )
    sp = sparse(CD.unzip(nonzero_inds)..., [1.0i+j for (i, j) in nonzero_inds])
    A1 = sp + tril(sp, -1)'

    # From VA Figure 4.2
    n2 = 17
    ijs = [(3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,3), (5,4), (8,7),
           (9,5), (9,6), (9,7), (9,8), (11,10), (13,10), (13,11), (13,12),
           (14,10), (14,11), (14,12), (14,13),
           (15, 1), (15,3), (15,4), (15,5), (15,7), (15,8), (15,9),
           (16,5), (16,6), (16,9), (16,12), (16,13), (16,14), (16,15),
           (17,10), (17,11), (17,12), (17,13), (17,14), (17,15), (17,16)]
    append!(ijs, [(i,i) for i in 1:n2])
    II, JJ = Chordal.unzip(ijs)
    A2 = sparse(II, JJ, ones(length(II)))
    A2 = (A2 + A2') / 2

    peo = CD.maximum_cardinality_search(A1)
    # @test peo == CD.maximum_cardinality_search(SimpleGraph(A1))
    @test CD.is_chordal(A1; peo=peo)
    # @test CD.is_chordal(SimpleGraph(A1); peo=1:n1)

    peo = CD.maximum_cardinality_search(A2)
    # @test peo == CD.maximum_cardinality_search(SimpleGraph(A2))
    @test CD.is_chordal(A2; peo=peo)
    # @test CD.is_chordal(SimpleGraph(A2); peo=1:n2)
end
