using CryoFaBs
using Healpix

nside = 8
npix = nside2npix(nside)

mask = fill(1.0, npix)

cfb = AngularCryoFaB(mask)
