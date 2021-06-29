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


@testset "Test Clique Tree" begin
    deg⁺_true = vec([4 2 3 2 3 2 3 2 2 4 3 4 3 2 2 1 0])

    et_par = CD.etree(L)
    et_ch = CD.get_children_from_par(et_par)
    deg⁺ = CD.get_higher_deg(L)
    @test all(deg⁺_true .== deg⁺)

    vreps, snd_par, snd_mem = CD.max_supernode_etree(L, et_par)
    # Test vrep (Thm 4.2)
    for v in vreps
        @test all([deg⁺[w] < deg⁺[v] + 1 for w in et_ch[v]])
    end
    # Test supernodes
    for (ind, v) in enumerate(vreps)
        nv = count(x->x==ind, snd_mem) - 1
        pv = v
        for k = 1:nv
            pv = et_par[pv]
            @test deg⁺[v] == deg⁺[pv] + k
            @test snd_mem[pv] == ind
        end
    end
    # Test clique tree (Thm 4.3)
    # Test root
    root_et = findfirst(x->x==0, et_par)
    root_ind = findfirst(x->x==0, snd_par)
    @test snd_mem[root_et] == root_ind
    for (ind, v) in enumerate(vreps)
        ind == root_ind && continue
        col_v = vcat([v], rowvals(L)[nzrange(L, v)])
        qv = vreps[snd_par[ind]]
        col_qv = vcat([qv], rowvals(L)[nzrange(L, qv)])
        col_minus_snd = filter(x->(snd_mem[x] != ind), col_v)
        # println(" --- v = $v --- ")
        # @show col_minus_snd
        # @show col_qv
        @test issubset(col_minus_snd, col_qv) && length(col_minus_snd) < length(col_qv)
    end

end
