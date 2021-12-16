# CryoFaBs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hsgg.github.io/CryoFaBs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hsgg.github.io/CryoFaBs.jl/dev)
[![Build Status](https://github.com/hsgg/CryoFaBs.jl/workflows/CI/badge.svg)](https://github.com/hsgg/CryoFaBs.jl/actions)

This code accompanies a paper on [Harmonic analysis of isotropic fields on the
sphere with arbitrary masks](https://arxiv.org/abs/2109.13352). It is primarily
intended to be used in the context of large galaxy surveys. We use the code to
derive custom eigenfunctions to the Laplacian on the sphere for arbitrary
domains. Also, custom eigenfunctions for 3D analysis separable in the angular
and radial directions.

The module can create a large cache in `~/.julia/scratchspaces/`, where it
caches the Green's matrix.
