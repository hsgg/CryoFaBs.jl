# 3D Tutorial

In this tutorial we will do everything to show how to create your own
*cryofunks* and use `CryoFaB` for a 3D analysis with spherical symmetry, but
arbitrary masks.

First, let's load some packages you will need:
```julia
using Healpix
using CryoFaBs
```

Next, we assume a galaxy survey with distances between $500\,h^{-1} {\rm Mpc}$
and $1500\,h^{-1}{\rm Mpc}$. For the angular resolution of the HEALPix map we
choose $n_{\rm side}=32$. Don't go too high in $n_{\rm side}$, as that will
take a lot of memory! For the radial resolution we choose 50 bins.
```julia
rmin = 500.0
rmax = 1500.0
nside = 32
nbins = 50
npix = nside2npix(nside)
mask = fill(1.0, npix)
mask[(npix ÷ 2):end] .= 0
mask = HealpixMap{Float64, Healpix.RingOrder}(mask)
cfb = AngRadCryoFaB(mask, rmin, rmax, nbins)
```
Here we created a mask that covers just about half the sky. In the last line we
create the `CryoFaB` that contains the *cryofunks* for the harmonic transform.
The function `nside2npix()` is provided by
[Healpix.jl](https://github.com/ziotom78/Healpix.jl).

To demonstrate the analysis of a 3D survey, let's create some $10^5$ galaxies
with random positions:
```julia
Ngals = 10^5
rθϕ = fill(NaN, 3, Ngals)
@. rθϕ[1,:] = rmin + (rmax - rmin) * rand()
@. rθϕ[2,:] = π/2 * rand()
@. rθϕ[3,:] = 2 * π * rand() - π
```
The coordinates are the radial distance, $θ$, and $ϕ$. This is a non-uniform
galaxy catalog. Therefore, we expect some funkyness in the radial direction. It
also isn't isotropic, so we really should find a better example. We have also
ensured that only the Northern cap is populated with galaxies, but if you need,
there is the `isinsurvey()` function to check if a given galaxy is in the
*CryoFaB*'s footprint.

In any case, the array `rθϕ` contains the $r$-coordinates of all the galaxies
in `rθϕ[1,:]`, all the $\theta$-coordinates are in `rθϕ[2,:]`, and all the
$\phi$-coordinates are in `rθϕ[3,:]`.

Next, we assign each galaxy to its nearest grid cell by calling
```julia
n_g = calc_numgals_per_cell(rθϕ, cfb)
δr = numgals2densitycontrast(n_g)
```
where the second line further estimates the density contrast. As if you
wouldn't be able to guess that from the name.

Now we can convert to cryospace,
```julia
δnlm = cfb \ δr
```
The power spectrum is constructed as
```julia
P = @. δnlm * conj(δnlm)
```
We have included the complex conjugate eventhough the *cryofunks* are real,
just to remind ourselves.

One more step is needed: the combination of pixel-window and survey geometry,
aka the *cryo-window*.

...
