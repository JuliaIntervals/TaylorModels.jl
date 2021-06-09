# TaylorModels

[![CI][CI_badge]][CI_url]
[![Coverage Status][Coverage_badge]][Coverage_url]

<!-- [![codecov badge][codecov_badge]][codecov_url] -->

[![docs stable][docsbadge_stable]][documenter_stable]
[![docs latest][docsbadge_latest]][documenter_latest]
[![DOI][Zenodo_badge]][Zenodo_url]

[CI_badge]: https://github.com/JuliaIntervals/TaylorModels.jl/actions/workflows/ci.yml/badge.svg
[CI_url]: https://github.com/JuliaIntervals/TaylorModels.jl/actions/workflows/ci.yml

[Coverage_badge]: https://coveralls.io/repos/github/JuliaIntervals/TaylorModels.jl/badge.svg?branch=lb/github_actions
[Coverage_url]: https://coveralls.io/github/JuliaIntervals/TaylorModels.jl?branch=lb/github_actions

[Zenodo_badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.2613102.svg
[Zenodo_url]: https://doi.org/10.5281/zenodo.2613102

<!-- [codecov_badge]: http://codecov.io/github/JuliaIntervals/TaylorModels.jl/coverage.svg?branch=master
[codecov_url]: http://codecov.io/github/JuliaIntervals/TaylorModels.jl?branch=master -->

[docsbadge_stable]: https://img.shields.io/badge/docs-stable-blue.svg
[documenter_stable]: https://juliaintervals.github.io/TaylorModels.jl/stable
[docsbadge_latest]: https://img.shields.io/badge/docs-latest-blue.svg
[documenter_latest]: https://juliaintervals.github.io/TaylorModels.jl/latest


This package combines the `IntervalArithmetic.jl` and `TaylorSeries.jl` packages to provide **Taylor models**, i.e.
Taylor polynomials with guaranteed error bounds to approximate functions.

An introduction is available in [this video from JuliaCon 2018](https://www.youtube.com/watch?v=o1h7BUW04NI).


### Authors
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders), Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)


## Bibliography

- Kyoko Makino, [Rigorous Analysis of Nonlinear Motion in Particle Accelerators](https://bt.pa.msu.edu/pub/papers/makinophd/makinophd.pdf), PhD Thesis, Michigan State University (1991).

- Mioara Joldes, [Rigorous Polynomial Approximations and Applications](https://tel.archives-ouvertes.fr/tel-00657843), PhD Thesis, École Normale Supérieure de Lyon, ENS-Lyon (2011).


## Acknowledgements ##

Financial support is acknowledged from DGAPA-UNAM PAPIIT grants IN-117117, IG-100616 and IG-100819. DPS acknowledges support through a *Cátedra Marcos Moshinsky* (2018).
