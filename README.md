# ChordalDecomp

[![Build Status](https://github.com/tjdiamandis/ChordalDecomp.jl/workflows/CI/badge.svg)](https://github.com/tjdiamandis/ChordalDecomp.jl/actions)
[![Coverage](https://codecov.io/gh/tjdiamandis/ChordalDecomp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tjdiamandis/ChordalDecomp.jl)

#### NOTE: This package is very much a work in progress. Feature requests welcome! Documentation forthcoming.

## Overview
This package implements data structures and subroutines for chordal matrices. It is largely a (work in progress) port of [`CHOMPACK.py`](https://chompack.readthedocs.io/en/latest/) to Julia and currently includes routines for
- Chordal Decomposition of PSD matrices
- Maximum Determinant PSD completion
- Chordality tests and perfect elimination ordering

These algorithms are (mostly) based on the multi-frontal algorithms described in [Logarithmic Barriers for Sparse Matrix Cones](https://arxiv.org/abs/1203.2742) and the survey paper [Chordal Graphs and Semidefinite Optimization](https://www.seas.ucla.edu/~vandenbe/publications/chordalsdp.pdf).

Additionally, `Chordal.jl` includes utility functions for the elimination tree and clique tree data structures described in these papers.


## Chordal Decomposition
Chordal decomposition breaks up a large `n x n` PSD matrix constraint into `K` smaller `n_k x n_k` PSD constraints. These smaller constraints correspond to the cliques in the chordal graph associated with the aggregate sparsity pattern of the problem (i.e., the matrices associated with the constraints and objective function involving the PSD variable).

These cliques can be further combined (replacing structural zeros with numeric zeros) to optimize optimization algorithm performance. This package uses the clique graph-based merging scheme introduced by Garstka et al. in [A clique graph based merging strategy for decomposable SDPs](https://arxiv.org/abs/1911.05615) and implemented in [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)).


## Additional functionality to include and TODOs
- [ ] Min rank matrix completion
- [ ] Utility functions for SDP decomposition with JuMP
- [ ] Symbolic and numeric cholesky factorization
- [ ] Allow user-specified weight algorithm for clique graph merging
- [ ] Euclidean distance matrix completion

## References
- Lieven Vandenberghe and Martin Andersen's [Chordal Graphs and Semidefinite Optimization](https://www.seas.ucla.edu/~vandenbe/publications/chordalsdp.pdf)
- Martin Andersen, Joachim Dahl, and Lieven Vandenberghe's [Logarithmic Barriers for Sparse Matrix Cones](https://arxiv.org/abs/1203.2742)
- Yifan Sun's thesis [Decomposition Methods for Semidefinite optimization](https://escholarship.org/content/qt1cv6981p/qt1cv6981p.pdf)
- Yifan Sun, Martin S. Andersen, and Lieven Vandenberghe's [Decomposition in conic optimization with partially separable structure](https://arxiv.org/abs/1306.0057)
- Michael Garstka, Mark Cannon, and Paul Goulart's [A clique graph based merging strategy for decomposable SDPs](https://arxiv.org/abs/1911.05615)


## See also
- [`COSMO.jl`](https://github.com/oxfordcontrol/COSMO.jl) is a conic solver that uses chordal decomposition for large PSD constraints.
- [`CHOMPACK.py`](https://github.com/cvxopt/chompack) includes the same algorithms in Python.
