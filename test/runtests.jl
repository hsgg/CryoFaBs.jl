using CryoFaBs
using Healpix
using BlockDiagonals
using Statistics
using Test


@testset "CryoFaBs.jl" begin
    @testset "AngularCryoFaB" begin
        nside = 8
        npix = nside2npix(nside)
        mask = fill(1.0, npix)
        cfb = AngularCryoFaB(mask; nside_hires=npix2nside(length(mask)))
    end


    @testset "Radial Cryofunks" begin
        ℓ = 3.14
        rmin = 500.0
        rmax = 1500.0
        nbins = 50
        Δr = (rmax - rmin) / nbins
        reff = range(rmin + Δr/2, rmax - Δr/2, length=nbins)
        k, g, ginv = CryoFaBs.get_radial_cryofunks(ℓ, reff, Δr)
    end


    @testset "AngRadCryoFaB" begin
        nside = 4
        rmin = 500.0
        rmax = 1500.0
        nbins = 50

        mask = fill(1.0, nside2npix(nside))

        cfb = AngRadCryoFaB(mask, rmin, rmax, nbins)

        for n=1:1:nbins, iL=1:10
            # the mode to excite:
            NL = (iL - 1) * nbins + n
            @show iL,n,NL

            # excite that mode:
            R = blocks(cfb.RadTrans)[iL]
            fofr = R[:,n]
            fofa = cfb.angfab.Trans[:,iL]
            δr = [fofr[N]*fofa[i] for N=1:nbins for i=1:cfb.angfab.npixsurvey]

            # manual transform
            δrlm = cfb.AngTransInv * δr
            δrlm_sort = δrlm[cfb.sortbyl]
            δnlm_new = cfb.RadTransInv * δrlm_sort

            # transform
            δnlm = cfb \ δr

            @test δnlm == δnlm_new

            idx = (1:length(δnlm))
            idx = idx[idx .!= NL]
            prec = √eps(eltype(δnlm))

            @show prec
            @show δnlm[NL]
            @show extrema(δnlm[idx])

            @test δnlm[NL] ≈ 1
            @test all(abs.(extrema(δnlm[idx])) .< prec)
        end


        # compile test
        @time k, ell, W = CryoFaBs.cryo_window(cfb; do_calc=true, lmax=10)
    end


    @testset "cryobin" begin
        println("No binning:")
        kk = collect(0:0.1:1)
        ll = fill(1, length(kk))
        P = collect(1:length(kk))
        @show kk ll P
        kbin, lbin, Pbin = cryobin(P, kk, ll)
        @show kbin lbin Pbin
        @test kbin == kk
        @test lbin == ll
        @test Pbin == P

        println("k binning, ktol=0.15:")
        kk = collect(0:0.1:1)
        ll = fill(1, length(kk))
        P = collect(1:length(kk))
        @show kk ll P
        kbin, lbin, Pbin = cryobin(P, kk, ll; ktol=0.15)
        kbin_th = [middle.(kk[1:2:end-1],kk[2:2:end]); [kk[end]]]
        lbin_th = fill(1, length(kbin))
        Pbin_th = [middle.(P[1:2:end-1],P[2:2:end]); [P[end]]]
        @show kbin lbin Pbin
        @show kbin_th lbin_th Pbin_th
        @test kbin == kbin_th
        @test lbin == lbin_th
        @test Pbin == Pbin_th

        println("l binning, ltol=1.5:")
        kk = fill(0.1, 10)
        ll = 1:length(kk)
        P = 1:length(kk)
        @show kk ll P
        kbin, lbin, Pbin = cryobin(P, kk, ll; ktol=0, ltol=1.5)
        lbin_th = middle.(ll[1:2:end-1],ll[2:2:end])
        kbin_th = fill(0.1, length(lbin))
        Pbin_th = middle.(P[1:2:end-1],P[2:2:end])
        @show kbin lbin Pbin
        @show kbin_th lbin_th Pbin_th
        @test kbin == kbin_th
        @test lbin == lbin_th
        @test Pbin == Pbin_th

        println("k & l binning:")
        kk = collect(0.1:0.1:1)
        ll = collect(1:length(kk))
        P = collect(1:length(kk))
        @show kk ll P
        kbin, lbin, Pbin = cryobin(P, kk, ll; ktol=0.15, ltol=1.5)
        kbin_th = collect(0.15:0.2:1)
        lbin_th = collect(1.5:2:10)
        Pbin_th = collect(1.5:2:10)
        @show kbin lbin Pbin
        @show kbin_th lbin_th Pbin_th
        @test kbin ≈ kbin_th
        @test lbin ≈ lbin_th
        @test Pbin ≈ Pbin_th
    end
end


# vim: set sw=4 et sts=4 :
