# Copyright 2021, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged. Any commercial
# use must be negotiated with the Office of Technology Transfer at the
# California Institute of Technology.
# 
# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

using CryoFaBs
using Healpix
using BlockDiagonals
using Statistics
using Test
using Random
using LinearAlgebra


@testset "CryoFaBs.jl" begin
    @testset "AngularCryoFaB" begin
        Random.seed!(1234567890)
        nside = 2
        npix = nside2npix(nside)
        #mask = fill(1.0, npix)
        mask = rand([0.0,1.0], npix)
        cfb = AngularCryoFaB(mask; nside_hires=npix2nside(length(mask)))

        fname = (@__DIR__) * "/test_rand_nside$nside.cfb"
        #write(fname, cfb)  # write reference
        cfb0 = AngularCryoFaB(fname)  # read reference

        for f in fieldnames(AngularCryoFaB)
            A = getfield(cfb, f)
            B = getfield(cfb0, f)
            @show f,norm(A .- B),norm(A.-B) / norm(A)
            @test getfield(cfb, f) ≈ getfield(cfb0, f)  rtol=1e-14
        end
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


    @testset "estimate_density_contrast" begin
        nside = 8
        rmin = 500.0
        rmax = 1500.0
        nbins = 50
        Ngals = 10^5

        mask = fill(1.0, nside2npix(nside))

        cfb = AngRadCryoFaB(mask, rmin, rmax, nbins)

        rθϕ = fill(NaN32, 3, Ngals)
        @. rθϕ[1,:] = rmin + (rmax - rmin) * rand()
        @. rθϕ[2,:] = π * rand()
        @. rθϕ[3,:] = 2 * π * rand() - π

        println("First:")
        @time n_g = CryoFaBs.calc_numgals_per_cell(rθϕ, cfb)
        println("Second:")
        @time n_g = CryoFaBs.calc_numgals_per_cell(rθϕ, cfb)

        @time n_g1 = CryoFaBs.calc_numgals_per_cell_v1(rθϕ, cfb)
        @test n_g == n_g1
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
