n = 17
ijs = [(3,1), (3,2), (4,1), (4,2), (4,3), (5,1), (5,3), (5,4), (8,7),
       (9,5), (9,6), (9,7), (9,8), (11,10), (13,10), (13,11), (13,12),
       (14,10), (14,11), (14,12), (14,13),
       (15, 1), (15,3), (15,4), (15,5), (15,7), (15,8), (15,9),
       (16,5), (16,6), (16,9), (16,12), (16,13), (16,14), (16,15),
       (17,10), (17,11), (17,12), (17,13), (17,14), (17,15), (17,16)]
append!(ijs, [(i,i) for i in 1:n])
I, J = ChordalDecomp.unzip(ijs)
A = sparse(I, J, ones(length(I)))
A = A + A'
L = qdldl(A; perm=nothing, logical=true).L

@testset "Test elimination tree" begin
    # VA survey, Figure 4.2
    par_true = vec([3 3 4 5 9 9 8 9 15 11 13 13 14 16 16 17 0])
    ch_true =  [Set{Int}() for _ in 1:n]
    ch_true[3] = Set([1, 2])
    ch_true[4] = Set([3])
    ch_true[5] = Set([4])
    ch_true[8] = Set([7])
    ch_true[9] = Set([5, 6, 8])
    ch_true[11] = Set([10])
    ch_true[13] = Set([11, 12])
    ch_true[14] = Set([13])
    ch_true[15] = Set([9])
    ch_true[16] = Set([15, 14])
    ch_true[17] = Set([16])

    et_par = CD.etree(L)
    @test all(et_par .== par_true)
    et_ch = CD.get_children_from_par(et_par)
    @test all([issetequal(et_ch[i], ch_true[i]) for i in 1:n])
    perm = CD.get_postordering(et_par, et_ch)
    iperm = invperm(perm)
    @test all([iperm[i] .< et_par[iperm[i]] for i in 1:17 if et_par[iperm[i]] ≠ 0])
end

