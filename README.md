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


## Copyright notice

Copyright (c) 2021 California Institute of Technology (“Caltech”).
U.S. Government sponsorship acknowledged.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

- Neither the name of Caltech nor its operating division, the Jet Propulsion
  Laboratory, nor the names of its contributors may be used to endorse or
  promote products derived from this software without specific prior written
  permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
