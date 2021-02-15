function make_selectors_from_cliques(cliques, n)
    selector_mats = [spzeros(length(cliques[i]), n) for i in 1:length(cliques)]
    for i in 1:length(cliques)
        idx = 1
        for node in cliques[i]
            selector_mats[i][idx, node] = 1.0
            idx += 1
        end
    end
    return selector_mats
end


function get_selectors(input_mat::SparseMatrixCSC)
    n = size(input_mat)[1]
    preprocess!(input_mat)

    sp = sparsity_pattern(input_mat)
    P, L = get_chordal_extension(sp)
    P = build_perm_matrix(P)

    cliques = get_cliques(L)
    cg = generate_clique_graph(cliques)
    merge_cliques!(cg; verbose=true)

    Cℓs = get_cliques(cg)
    Tℓs = make_selectors_from_cliques(Cℓs, n)
    return P, Cℓs, Tℓs
end
