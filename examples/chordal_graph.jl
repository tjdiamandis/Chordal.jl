using ChordalDecomp
using LightGraphs, MetaGraphs
import Plots
import GraphPlot

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
# Example from A Clique Graph Based Merging Strategy for Decomposable SDPs
# by Michael Garstka, Mark Cannon, Paul Goulart
#       See Figures 1 and 3
#       Note: This does not include chordal extension (already chordal)

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
@info join(vcat(["Cliques are:"], ["\n\t\tC$i: $(cliques[i])" for i in 1:length(cliques)]))
GraphPlot.gplot(sparsity_graph, nodelabel=1:nv(sparsity_graph))


## Generate clique graph
cg = generate_clique_graph(cliques)
@info ("Initial cliques: $([string(a) for a in get_cliques(cg)])")
node_labels = ["$i, $(get_prop(cg, i, :nodes))" for i in 1:nv(cg)]
edge_labels = [get_prop(cg, e, :weight) for e in edges(cg)]
plt = GraphPlot.gplot(
        cg,
        nodelabel=node_labels,
        edgelabel=edge_labels,
)
# c.f. Figure 3


## Merge cliques
merge_cliques!(cg; verbose=true)
@info ("Final cliques: $([string(a) for a in get_cliques(cg)])")

node_labels = ["$(get_prop(cg, i, :nodes))" for i in 1:nv(cg)]
edge_labels = [get_prop(cg, e, :weight) for e in edges(cg)]
plt = GraphPlot.gplot(
        cg,
        nodelabel=node_labels,
        edgelabel=edge_labels,
)
# c.f. Figure 3




## Selector Matrices Testing
Cℓs = get_cliques(cg)
Tℓs = make_selectors(Cℓs, n)

Cℓs[1]
Tℓs[1]' * sparse([1 1 ; 1 1]) * Tℓs[1]
