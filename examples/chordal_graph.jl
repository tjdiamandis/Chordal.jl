using Pkg
Pkg.activate(@__DIR__)
using Chordal
using LinearAlgebra, SparseArrays
using Plots: spy
using Graphs, GraphPlot, Cairo, Colors

function show_cliques(mat, cliques)
    nzs = findnz(mat)[1:2]
    for clique in cliques
        sp = sparse(nzs..., [(i in clique && j in clique)  ? 10 : 1 for (i, j) in Chordal.zip(nzs...)])
        plt = spy(sp, ms=10, colorbar=false)
        display(plt)
    end
end

function plot_clique_graph(cg::Chordal.CliqueGraph)
        nc = length(cg.active_cliques)
        ind = sort(collect(cg.active_cliques))
        reduced_edge_mat = cg.edge_mat[ind,ind]
        clique_graph = SimpleGraph(reduced_edge_mat)
        node_labels = ["$c" for (i, c) in enumerate(get_cliques(cg))]
        edge_labels = [string(reduced_edge_mat[src(e),dst(e)]) for e in Graphs.edges(clique_graph)]
        plt = gplot(
                clique_graph,
                NODESIZE=0.3/sqrt(nc),
                nodelabel=node_labels,
                nodelabelc=colorant"blue",
                NODELABELSIZE=3,
                # nodesize=100*ones(length(Graphs.vertices(clique_graph))),
                nodefillc=colorant"white",
                edgelabel=edge_labels,
                edgelabelc=colorant"red",
                EDGELINEWIDTH=0.5/sqrt(nc)
        )
end

## Example from A Clique Graph Based Merging Strategy for Decomposable SDPs
# by Michael Garstka, Mark Cannon, Paul Goulart
#       See Figures 1 and 3
#       Note: This does not include chordal extension (already chordal)
# Setup sparsity pattern
n = 9
nonzero_inds = vcat([
        (3,1), (3,2), (5,4), (6,1), (6,3), (7,3), (7,6), (8,3), (8,4), (8,5),
        (8,6), (8,7), (9,6), (9,7), (9,8)
    ],
    [(i,i) for i in 1:n]
)
sp = sparse(Chordal.unzip(nonzero_inds)..., [1 for (i, j) in nonzero_inds])
spy(sp, ms=10)
sp = sparse(Chordal.unzip(nonzero_inds)..., [10i+j for (i, j) in nonzero_inds])
# c.f. Figure 1


## Find cliques
sp = sp + tril(sp, -1)'
Chordal.preprocess!(sp)

sparsity_graph = SimpleGraph(sp)
gplot(sparsity_graph, nodelabel=1:nv(sparsity_graph))
cliques = maximal_cliques(sparsity_graph)
@info join(vcat(["Cliques are:"], ["\n\t\tC$i: $(cliques[i])" for i in 1:length(cliques)]))
show_cliques(sp, cliques)


## Generate clique graph
cg = generate_clique_graph(cliques, n)
@info ("Initial cliques: $([string(a) for a in get_cliques(cg)])")
plot_clique_graph(cg)
# c.f. Figure 3


## Merge cliques
merge_cliques!(cg; verbose=true)
@info ("Final cliques: $([string(a) for a in get_cliques(cg)])")
plot_clique_graph(cg)
# c.f. Figure 3


## Selector Matrices
mat = Matrix(sp)
Cℓs = get_cliques(cg)
Tℓs = make_selectors_from_cliques(Cℓs, n)
println("Full Matrix:")
display(sp)
for i in 1:length(Cℓs)
        println("\nClique $i: $(string(Cℓs[i]))")
        display((Tℓs[i]*mat*Tℓs[i]'))
end