# 2D Tutorial

Here we give a brief introduction on how to use `CryoFaB` on the 2D sphere.

```julia
using Healpix
using CryoFaBs
```

Then, create a mask, or rather, in inverse of the mask where the observed
pixels are ones and the masked pixels are zero. For example:

```julia
nside = 4

npix = nside2npix(nside)

mask_Float64 = Float64.(rand(Bool, npix))

mask = HealpixMap{Float64, Healpix.RingOrder}(mask_Float64)
```

That allows us to create an `AngularCryoFaB` that contains all the information
about the transform:
```julia
cfb = AngularCryoFaB(mask)
```

The object `cfb` acts like a matrix. It can be queried with some standard
functions. For example,
```julia
size(cfb)
```
Whether a point on the sky is in the survey or which index in the input field
it has:
```julia
yes = isinsurvey(cfb, theta, phi)
idx = coord2cell(cfb, theta, phi)
```

The angular *cryofab* acts on fields on the incomplete sphere. To get such a
field, we first create a random field on the full sky, then select only those
pixels that are not masked:
```julia
delta_rhat = rand(npix)[mask .> 0]
```
The function `coord2cell()` can also be very helpful in constructing
the field `delta_rhat`, which is an array with length the number of pixels that
are observed.

To transform from configuration space to cryospace,
```julia
delta_lm = cfb \ delta_rhat
```

Compute the power spectrum:
```julia
Cl = @. delta_lm * conj(delta_lm)
```
This will be afflicted by the pixel+geometry window, which we calculate next.


## Pixel and Geometry Window

Calculate the pixel and geometry window:
```julia
ellavg, pixwin = angular_cryo_window(cfb)
```
Then, given a prediction `Cl_theory` evaluated at modes `ellavg`, the
window-convolved power spectrum would be
```julia
Cl_theory_mask = pixwin * Cl_theory
```
to be compared with the measured `Cl` at modes `cfb.ell`.
