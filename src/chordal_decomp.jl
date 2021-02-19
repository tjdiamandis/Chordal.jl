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


# TODO: make this multithreaded
function make_selectors_from_clique_graph(cg::CliqueGraph, n)
    m = length(cg.active_cliques)
    selector_mats = Vector{SparseMatrixCSC}(undef, m)
    for (i, cnum) in enumerate(cg.active_cliques)
        clique = cg.membership_mat[:,cnum]
        selector_mat = spzeros(length(clique), n)

        idx = 1
        for node in findnz(clique)[1]
            selector_mat[idx, node] = 1.0
            idx += 1
        end
        selector_mats[i] = selector_mat
    end
    return selector_mats
end


function get_selectors(input_mat::SparseMatrixCSC; verbose=true, ret_cliques=true)
    n = size(input_mat)[1]
    preprocess!(input_mat)

    sp = sparsity_pattern(input_mat)
    P, L = get_chordal_extension(sp; verbose=verbose)
    P = build_perm_matrix(P)

    cliques = get_cliques(L)
    cg = generate_clique_graph(cliques, n)
    merge_cliques!(cg; verbose=verbose)

    if ret_cliques
        Cℓs = get_cliques(cg)
        Tℓs = make_selectors_from_cliques(Cℓs, n)
        return P, Cℓs, Tℓs
    end

    Tℓs = make_selectors_from_clique_graph(cg, n)
    return P, Tℓs
end
