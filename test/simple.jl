using CryoFaBs
using Healpix
using Profile, PProf
using Random


Random.seed!(1234567890)


nside = 8

npix = nside2npix(nside)

mask = rand([0.0,1.0], npix)
#mask = fill(1.0, npix)

@show nside npix sum(mask)
println("Expected size of Greens matrix: ", sum(mask)^2 * 8 / 1024^2, " MB")

@show Sys.maxrss() / 1024^3
AngularCryoFaB(mask)

@show Sys.maxrss() / 1024^3

Profile.Allocs.clear()

Profile.Allocs.@profile sample_rate=1 AngularCryoFaB(mask)

PProf.Allocs.pprof()

@show Sys.maxrss() / 1024^3

#:Profile.clear_alloc_data()


# vim: set sw=4 et sts=4 :
