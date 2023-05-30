```@meta
CurrentModule = CryoFaBs
```

# CryoFaBs

This is the documentation for [CryoFaBs](https://github.com/hsgg/CryoFaBs.jl).
The paper describing the details of the algorithm will be on arXiv shortly.

This package is similar in spirit to
[Healpix.jl](https://github.com/ziotom78/Healpix.jl) and
[SphericalFourierBesselDecompositions.jl](https://github.com/hsgg/SphericalFourierBesselDecompositions.jl).
It performs a harmonic analysis (usually of galaxy surveys) exploiting the
isotropy on the sky despite non-isotropic window functions and masks. It does
this in both 2D and 3D.

## Installation

To install, add the github repo to your julia project, like so:
```julia
julia>] add https://github.com/hsgg/CryoFaBs.jl.git
```


## Basic usage

The central objects of this package are of the abstract type `CryoFaB`, which
behaves like a matrix. For example, if `δr` is the density contrast in
configuration space and `cfb` is a *CryoFaB*, then the transform into cryospace
is performed via
```julia
δnlm = cfb \ δr
```

`CryoFaB`-objects like `AngularCryoFaB` (for 2D on the sphere) and
`AngRadCryoFaB` (for 3D) are also used to convert between pixels (or voxels in
3D) and a cell index, e.g.,
```julia
pix = coord2cell(cfb, θ, ϕ)
```
This is similar to the function `ang2pix()` from
[Healpix.jl](https://github.com/ziotom78/Healpix.jl) and indeed it is the same
in the 2D case. In 3D the same is done by using a different `cfb` object of
type `AngRadCryoFaB` and additionally passing a radial coordinate.

The `CryoFaB`-objects also encode everything that is needed for the pixel
window.
