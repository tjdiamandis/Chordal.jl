using ChordalDecomp
import Plots
import GraphPlot

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
# Example from A Clique Graph Based Merging Strategy for Decomposable SDPs
# by Michael Garstka, Mark Cannon, Paul Goulart
# Figures 1 and 3

## Setup sparsity pattern
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
sp = sparse(unzip(nonzero_inds)..., ones(length(nonzero_inds)))

# c.f. Figure 1
Plots.spy(sp, ms=10)


## Find cliques
sp = sp + sp'
preprocess!(sp)

sparsity_graph = SimpleGraph(sp)
cliques = maximal_cliques(sparsity_graph)

println("Cliques are:")
[println("C$i: $(cliques[i])") for i in 1:length(cliques)]
GraphPlot.gplot(sparsity_graph, nodelabel=1:nv(sparsity_graph))


## Generate clique graph
cg = generate_clique_graph(cliques)
node_labels = ["$i, $(get_prop(cg, i, :nodes))" for i in 1:nv(cg)]
edge_labels = [get_prop(cg, e, :weight) for e in edges(cg)]
plt = GraphPlot.gplot(
        cg,
        nodelabel=node_labels,
        edgelabel=edge_labels,
)


## Merge cliques
merge_cliques!(cg; verbose=true)
node_labels = ["$(get_prop(cg, i, :nodes))" for i in 1:nv(cg)]
edge_labels = [get_prop(cg, e, :weight) for e in edges(cg)]

# c.f. Figure 3
plt = GraphPlot.gplot(
        cg,
        nodelabel=node_labels,
        edgelabel=edge_labels,
)
