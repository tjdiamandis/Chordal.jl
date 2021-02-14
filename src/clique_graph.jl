using LightGraphs: src, dst, get_edge
using MetaGraphs: MetaGraph, set_prop!, get_prop

function generate_clique_graph(cliques)
    N = length(cliques)
    cliques = [Set(c) for c in cliques]
    # G = SimpleWeightedGraph(N)
    G = MetaGraph(N)
    for i in 1:N
        set_prop!(G, i, :nodes, cliques[i])
    end
    for i in 1:N, j in i+1:N
        # TODO: iterate through for performance?
        l_int = length(intersect(cliques[i], cliques[j]))
        if l_int > 0
            w_ij = weight_function(cliques[i], cliques[j], l_int)
            # add_edge!(G, i, j, w_ij)
            add_edge!(G, i, j)
            set_prop!(G, i, j, :weight, w_ij)
        end
    end
    return G
end


function weight_function(c_i, c_j, l_int)
    l_ci = length(c_i)
    l_cj = length(c_j)
    l_union = l_ci + l_cj - l_int
    return l_ci^3 + l_cj^3 - l_union^3
end


function _merge_cliques!(cg, i, j)
    new_clique = union(get_prop(cg, i, :nodes), get_prop(cg, j, :nodes))
    n_rem = degree(cg, i) < degree(cg, j) ? i : j
    n_keep = degree(cg, i) < degree(cg, j) ? j : i

    set_prop!(cg, n_keep, :nodes, new_clique)

    rem_edge!(cg, i, j)
    for neighbor in union(neighbors(cg, n_rem), neighbors(cg, n_keep))
        # Some redundancy here. Idk if it matters much
        rem_edge!(cg, n_rem, neighbor)
        add_edge!(cg, neighbor, n_keep)

        # Update weight
        l_int = length(intersect(new_clique, get_prop(cg, neighbor, :nodes)))
        w_ij = weight_function(new_clique, get_prop(cg, neighbor, :nodes), l_int)
        set_prop!(cg, neighbor, n_keep, :weight, w_ij)
    end
    rem_vertex!(cg, n_rem)
end


function merge_cliques!(cg; verbose=false)
    max_val, max_ind = findmax([get_prop(cg, e, :weight) for e in edges(cg)])

    verbose && @info "Starting merge"
    while max_val > 0
        i, j = get_edge(max_ind, edges(cg))
        if verbose
            @info "Merging cliques $(get_prop(cg, i, :nodes)) and $(get_prop(cg, j, :nodes))\n\t\t\tweight: $max_val"
        end

        _merge_cliques!(cg, i, j)
        max_val, max_ind = findmax([get_prop(cg, e, :weight) for e in edges(cg)])
    end
    verbose && @info "Finished merging; max weight = $max_val"
end


function get_edge(max_ind, edge_iter)
    count = 1
    for e in edge_iter
        count == max_ind && return src(e), dst(e)
        count += 1
    end
end


function get_cliques(cg)
    return [get_prop(cg, i, :nodes) for i in 1:nv(cg)]
end
