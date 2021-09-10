#!/usr/bin/env julia

@doc raw"""
This module implements functions to create your own fabulous basis functions,
    primarily for the analysis of galaxy redshift surveys.
"""
module CryoFaBs


export CryoFaB, AngularCryoFaB, AngRadCryoFaB
export isinsurvey
export estimate_density_contrast
export angular_cryo_window
export cryobin


using LinearAlgebra
using BlockDiagonals
using Healpix
using Statistics
using SparseArrays
using SpecialFunctions

include("Splines.jl")
using .Splines

using TimerOutputs

#include("SciPy.jl")
#using .SciPy
using FastTransforms
#using SphericalHarmonics

using Scratch
cachedir_greensmatrix = ""
function __init__()
    global cachedir_greensmatrix = @get_scratch!("greensmatrix")
end


################### CryoFaB basis ######################

@doc raw"""
The abstract type `CryoFaB` defines a common interface for the bases defined in
        this module.

            δnlm = cfb \ δθϕr[:]
"""
abstract type CryoFaB end  # Create Your Own Fabulous Basis
# TODO: CryoFaB should probably be an AbstractMatrix.


struct AngularCryoFaB <: CryoFaB
    # define angular mask
    mask::HealpixMap{Float64}
    nside::Int
    npixsurvey::Int
    hpix_full2survey::Array{Int,1}  # surveypixel = hpix_full2survey[fullskypixel]

    # transform
    ell::Vector{Float64}
    Trans::Matrix{Float64}
    TransInv::Matrix{Float64}
end


struct AngRadCryoFaB <: CryoFaB  # CryoFaB that separates angular and radial basis
    # angular transform
    angfab::AngularCryoFaB
    Angl::Vector{Float64}
    AngTrans::BlockDiagonal{Float64}
    AngTransInv::BlockDiagonal{Float64}

    sortbyl::Array{Int,1}  # to transpose the input vector between angular and radial transforms

    # define radial boundary
    rmin::Float64
    rmax::Float64
    rbounds::AbstractArray{Float64,1}  # bin boundaries
    reff::AbstractArray{Float64,1}  # effective center of bin

    # radial transform
    RadK::Array{Float64,1}
    RadTrans::BlockDiagonal{Float64}
    RadTransInv::BlockDiagonal{Float64}
end


include("io.jl")
include("cryofunks_angular.jl")
include("cryofunks_radial.jl")
include("cryowindow_angular.jl")
include("cryowindow_radial.jl")
include("cryowindow_sfb.jl")


function AngularCryoFaB(mask; nside_hires=npix2nside(length(mask)), method=:cryofunk)
    if !(typeof(mask) <: HealpixMap)
        mask = HealpixMap{Float64,Healpix.RingOrder}(mask)
    end
    nside = npix2nside(length(mask))
    npix = nside2npix(nside)
    npixsurvey = sum(mask .> 0)
    hpix_full2survey = fill(0, npix)
    @show typeof(hpix_full2survey) size(hpix_full2survey) npix
    p = 1
    for i=1:npix
        if mask[i] > 0
            hpix_full2survey[i] = p
            p += 1
        end
    end
    if method == :cryofunk
        ℓ, Y, Yinv = get_angular_cryofunks(mask; nside_hires=nside_hires)
    elseif method == :integerLM
        ℓ, Y, Yinv = get_angular_continuous_cryofunks(mask; nside_hires=nside_hires)
    elseif method == :continuousfunks
        ℓ, = get_angular_cryofunks(mask; nside_hires=nside_hires)
        Y, Yinv = get_angular_continuous_cryofunks(ℓ, mask; nside_hires=nside_hires)
    end
    AngularCryoFaB(mask, nside, npixsurvey, hpix_full2survey, ℓ, Y, Yinv)
end


function block_vector(v, n)
    bv = eltype(v)[]
    for i=1:n
        append!(bv, v)
    end
    return bv
end


function AngRadCryoFaB(mask, rmin, rmax, nbins)
    # angular part
    angfab = AngularCryoFaB(mask)
    Angℓ = block_vector(angfab.ell, nbins)
    AngTrans = BlockDiagonal(fill(angfab.Trans, nbins))
    AngTransInv = BlockDiagonal(fill(angfab.TransInv, nbins))

    # radial part
    Δr = (rmax - rmin) / nbins
    reff = range(rmin + Δr/2, rmax - Δr/2, length=nbins)
    kn = []
    gn = Matrix{Float64}[]
    gninv = Matrix{Float64}[]
    @timeit "radial cryofunk" for l=1:length(angfab.ell)
        k, g, ginv = get_radial_cryofunks(angfab.ell[l], reff, Δr)
        append!(kn, k)
        push!(gn, g)
        push!(gninv, ginv)
    end
    RadK = collect(kn)
    RadTrans = BlockDiagonal(gn)
    RadTransInv = BlockDiagonal(gninv)

    #sortbyl = sortperm(Angℓ)  # No! Need to keep the order of shells for equal ℓ the same.
    sortbyl = Int[]
    for iLM=1:length(angfab.ell)
        for ik=1:nbins
            j = (ik - 1) * angfab.npixsurvey + iLM
            push!(sortbyl, j)
        end
    end

    rbounds = range(rmin, rmax, length=nbins+1)

    AngRadCryoFaB(angfab, Angℓ, AngTrans, AngTransInv, sortbyl,
        rmin, rmax, rbounds, reff, RadK, RadTrans, RadTransInv)
end


Base.size(cfb::AngularCryoFaB) = size(cfb.Trans)
Base.size(cfb::AngularCryoFaB, i) = size(cfb.Trans,i)

Base.size(cfb::AngRadCryoFaB) = size(cfb.RadTrans,1), size(cfb.AngTrans,2)
Base.size(cfb::AngRadCryoFaB, i) = (i==1) ? size(cfb.RadTrans,1) : size(cfb.AngTrans,2)


function isinsurvey(cfb::AngularCryoFaB, θ, ϕ)
    p = ang2pix(cfb.mask, θ, ϕ)
    p = cfb.hpix_full2survey[p]
    if p == 0
        return false
    end
    return true
end


function coord2cell(cfb::AngularCryoFaB, θ, ϕ)
    p = ang2pix(cfb.mask, θ, ϕ)
    p = cfb.hpix_full2survey[p]
    if p == 0
        error("angle outside of mask: (θ,ϕ)=($θ,$ϕ)")
    end
    return p
end

function bin2coord(cfb::AngularCryoFaB, p)
    error("unimplemented")
end

function mode2idx(cfb::AngularCryoFaB, l)
    error("unimplemented")
end

function idx2mode(cfb::AngularCryoFaB, i)
    return cfb.ell[i]
end


function isinsurvey(cfb::AngRadCryoFaB, r, θ, ϕ)
    if !isinsurvey(cfb.angfab, θ, ϕ)
        return false
    end

    rbounds = get_rbounds(cfb)
    i = (r == cfb.rmax) ? length(cfb.reff) : searchsortedlast(rbounds, r)
    if !(1 <= i <= length(cfb.reff))
        return false
    end
    return true
end


function coord2cell(cfb::AngRadCryoFaB, r, θ, ϕ)
    p = coord2cell(cfb.angfab, θ, ϕ)

    rbounds = get_rbounds(cfb)
    i = (r == cfb.rmax) ? length(cfb.reff) : searchsortedlast(rbounds, r)
    if !(1 <= i <= length(cfb.reff))
        error("radial dist outside of survey: (r,θ,ϕ)=($r,$θ,$ϕ), $i, $(length(cfb.reff))")
    end

    v = p + (i - 1) * cfb.angfab.npixsurvey
    if v > size(cfb,2)
        @error r,θ,ϕ p i v size(cfb) cfb.angfab.npixsurvey
    end
    return v
end


function bin2coord(cfb::AngRadCryoFaB, v)
    error("unimplemented")
end


function idx2mode(cfb::AngRadCryoFaB, i)
    return cfb.RadK[i], cfb.Angl[cfb.sortbyl][i]
end


function mode2idx(cfb::AngRadCryoFaB, k, l)
    error("unimplemented")
end


Base.:\(cfb::AngularCryoFaB, v) = cfb.TransInv * v
Base.:*(cfb::AngularCryoFaB, v) = cfb.Trans * v
Base.adjoint(cfb::AngularCryoFaB) = cfb.TransInv'
Base.inv(cfb::AngularCryoFaB) = cfb.TransInv


Base.:\(cfb::AngRadCryoFaB, v::AbstractVector) = begin
    cfb.RadTransInv * (cfb.AngTransInv * v)[cfb.sortbyl]
end

Base.:*(cfb::AngRadCryoFaB, v::AbstractVector) = begin
    error("unimplemented")
end


function get_pixs_at_r(cfb::AngRadCryoFaB, i)
    return (1:cfb.angfab.npixsurvey) .+ (i-1) * cfb.angfab.npixsurvey
end


function get_rbounds(cfb::AngRadCryoFaB)
    range(cfb.rmin, cfb.rmax, length=length(cfb.reff)+1)
end



################### estimate density contrast ##########

@doc raw"""
calc_numgals_per_cell(rθϕ, cfb::CryoFaB)

Calculate the number of galaxies in `rθϕ` in each pixel or voxel in the survey
that `cfb` can act on. The radial coordinates are in `rθϕ[1,:]`, etc.
"""
function calc_numgals_per_cell(rθϕ, cfb::CryoFaB)
    T = Float64  # need to allow floats for easy conversion to a density contrast
    N = fill(T(0), size(cfb,1))
    for i=1:size(rθϕ,2)
        p = coord2cell(cfb, rθϕ[:,i]...)
        #p = coord2cell(cfb, rθϕ[1,i], rθϕ[2,i], rθϕ[3,i])
        N[p] += 1
    end
    return N
end


function numgals2densitycontrast!(δ, cfb::AngRadCryoFaB)
    rbounds = get_rbounds(cfb)
    for i=2:length(rbounds)
        p = get_pixs_at_r(cfb, i-1)
        ntotz = sum(δ[p])
        nofz = ntotz / length(p)  # average number of galaxies per pixel
        @. δ[p] = δ[p] / nofz - 1
    end
    return δ
end

function numgals2densitycontrast(δ, cfb::AngRadCryoFaB)
    return numgals2densitycontrast!(deepcopy(δ), cfb)
end


@doc raw"""
estimate_density_contrast(rθϕ, cfb::CryoFaB)

Estimate the density contrast from the galaxy positions in `rθϕ`.
"""
function estimate_density_contrast(rθϕ, cfb::AngRadCryoFaB)
    δ = calc_numgals_per_cell(rθϕ, cfb)
    numgals2densitycontrast!(δ, cfb)
    return δ
end


function calc_angular_average(nr, cfb::CryoFaB)
    rbounds = get_rbounds(cfb)
    nbar = fill(float(eltype(nr)(NaN)), length(rbounds) - 1)
    for i=1:length(rbounds)-1
        p = get_pixs_at_r(cfb, i)
        nbar[i] = mean(nr[p])
    end
    return nbar
end


function cryobin(P, cfb::AngRadCryoFaB; ltol=0.2, ktol=minimum(cfb.RadK), ignore_l=[])
    kk, ll = idx2mode(cfb, 1:length(P))

    # ignore ℓ=0 and ℓ>2*nside
    s = fill(true, length(ll))
    sl = @. cfb.angfab.ell > 1.5 * cfb.angfab.nside
    ignore_l = [ignore_l; cfb.angfab.ell[sl]]
    for l in ignore_l
        @. s &= (ll != l)
    end

    # ignore large k
    Δr = (cfb.rmax - cfb.rmin) / length(cfb.reff)
    @. s &= (kk < 3 / Δr)

    return cryobin(P[s], kk[s], ll[s]; ltol=ltol, ktol=ktol)
end

function cryobin(P, kk, ll; ltol=0.5, ktol=minimum(abs,kk))
    kbin = float(eltype(kk))[]
    lbin = float(eltype(ll))[]
    Pbin = float(eltype(P))[]
    Nbin = Int[]
    while length(kk) >= 1
        k = kk[1]
        l = ll[1]
        sel = @. (abs(kk - k) <= ktol) & (abs(ll - l) <= ltol)
        #@show length(sel),sum(sel)
        #@show k,l kk kk[sel] ll[sel] P[sel]
        sel[1] = true  # include at least one in case 'sel' is false
        push!(kbin, mean(kk[sel]))
        push!(lbin, mean(ll[sel]))
        push!(Pbin, mean(P[sel]))
        push!(Nbin, sum(sel))
        kk = kk[.!sel]
        ll = ll[.!sel]
        P = P[.!sel]
    end
    return kbin, lbin, Pbin, Nbin
end


end


# vim: set sw=4 et sts=4 :
