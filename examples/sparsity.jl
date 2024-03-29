using Pkg
Pkg.activate(@__DIR__)
using Chordal
using LinearAlgebra, SparseArrays, Random
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

## Banded setup
n = 15
w = 3
banded = [(i,j) for i in 1:n for j in 1:i if i - j < w]
sp = sparse(Chordal.unzip(banded)..., ones(length(banded)))
sp = sp + tril(sp, -1)'
spy(sp, ms=5, title="Banded")

## Banded sparsity graph
@show Chordal.is_chordal(sp)
sparsity_graph = SimpleGraph(sp)
gplot(sparsity_graph, nodelabel=1:nv(sparsity_graph))

## Banded cliques
cliques = [sort(x) for x in sort(maximal_cliques(sparsity_graph), by=x->minimum(x))]
@info join(vcat(["Cliques are:"], ["\n\t\tC$i: $(cliques[i])" for i in 1:length(cliques)]))
show_cliques(sp, cliques)


## Arrow setup
w = 1
arrow = [(i,j) for i in 1:n for j in 1:i if i - j <= w || n - i < w]
sp = sparse(Chordal.unzip(arrow)..., ones(length(arrow)))
sp = sp + tril(sp, -1)'
spy(sp, ms=5, title="Arrow")

## Arrow sparsity graph
@show Chordal.is_chordal(sp)
sparsity_graph = SimpleGraph(sp)
gplot(sparsity_graph, layout=spectral_layout, nodelabel=1:nv(sparsity_graph))

## Arrow cliques
cliques = [sort(x) for x in sort(maximal_cliques(sparsity_graph), by=x->minimum(x))]
@info join(vcat(["Cliques are:"], ["\n\t\tC$i: $(cliques[i])" for i in 1:length(cliques)]))
show_cliques(sp, cliques)


## Random
Random.seed!(0)
sp = tril(sprand(Bool, n, n, 0.2))
sp[diagind(sp)] .= 1.0
sp = Float64.(sp)
sp = sp + tril(sp, -1)'

spy(sp, ms=7, title="Random")
@show Chordal.is_chordal(sp)
sp_ext = Chordal.get_chordal_extension(sp; perm=nothing, verbose=true)[3]
sp_ext = sp_ext + sp_ext'
sp_ext[diagind(sp_ext)] .= 1.0
for (i, j, _) in Chordal.zip(findnz(sp_ext)...)
    if iszero(sp[i,j])
        sp_ext[i,j] = 5.0
    end
end
spy(sp_ext, ms=7, title="Random - Chordal Extension", colorbar=false, dpi=300)

## Banded sparsity graph
sparsity_graph = SimpleGraph(sp_ext)
gplot(sparsity_graph, nodelabel=1:nv(sparsity_graph))

## Banded cliques
cliques = [sort(x) for x in sort(maximal_cliques(sparsity_graph), by=x->minimum(x))]
@info join(vcat(["Cliques are:"], ["\n\t\tC$i: $(cliques[i])" for i in 1:length(cliques)]))
show_cliques(sp_ext, cliques)