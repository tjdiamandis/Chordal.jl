# Example from A Clique Graph Based Merging Strategy for Decomposable SDPs
# by Michael Garstka, Mark Cannon, Paul Goulart
#       See Figures 1 and 3
#       Note: This does not include chordal extension (already chordal)

# Setup sparsity pattern
n = 9
nonzero_inds = vcat([
        (3,1),
        (3,2),
        (5,4),
        (6,1),
        (6,3),
        (7,3),
        (7,6),
        (8,3),
        (8,4),
        (8,5),
        (8,6),
        (8,7),
        (9,6),
        (9,7),
        (9,8)
    ],
    [(i,i) for i in 1:n]
)
sp = sparse(CD.unzip(nonzero_inds)..., [10i+j for (i, j) in nonzero_inds])
sp = sp + tril(sp, -1)'


@testset "Clique Graph Construction" begin
    cliques_true = Set([
        Set([6, 7, 8, 9]),
        Set([6, 7, 8, 3]),
        Set([8, 4, 5]),
        Set([3, 6, 1]),
        Set([3, 2])
    ])
    final_cliques_true = Set([
        Set([3, 6, 7, 8, 9]),
        Set([8, 4, 5]),
        Set([3, 6, 1]),
        Set([3, 2])
    ])

    sp_full = sparse(Matrix(sp))
    CD.preprocess!(sp_full)
    @test !(0.0 in nonzeros(sp_full)) && all(sp_full .== sp)

    cliques = Set([Set(c) for c in get_cliques(tril(sp))])
    @test issetequal(cliques, cliques_true)

    cg = generate_clique_graph(cliques, n)
    cliques = Set([Set(c) for c in get_cliques(cg)])
    @test issetequal(cliques, cliques_true)

    merge_cliques!(cg; verbose=false)
    cliques = Set([Set(c) for c in get_cliques(cg)])
    @test issetequal(cliques, final_cliques_true)

    # TODO: test selectors

end
