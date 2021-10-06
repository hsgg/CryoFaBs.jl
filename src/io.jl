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

struct FormatError <: Exception end

Base.showerror(io::IO, e::FormatError) = print(io, "Error parsing file.")

const formatversion_angfab = 1
const formatversion_angradfab = 3


######## I/O AngularCryoFaB
function Base.write(io::IO, cfb::AngularCryoFaB)
    n = 0
    # assumes RingOrder
    n += write(io, Int64(formatversion_angfab))
    n += write(io, Int64(cfb.nside))
    n += write(io, Array{Float64,1}(cfb.mask))
    n += write(io, Int64(cfb.npixsurvey))
    n += write(io, Array{Int64,1}(cfb.hpix_full2survey))
    n += write(io, Vector{Float64}(cfb.ell))
    n += write(io, Matrix{Float64}(cfb.Trans))
    n += write(io, Matrix{Float64}(cfb.TransInv))
    #@show size(cfb.ell) size(cfb.Trans) size(cfb.TransInv)
    return n
end


function Base.read(io::IO, ::Type{AngularCryoFaB})
    # assumes RingOrder
    diskversion = read(io, Int64)
    if diskversion != formatversion_angfab
        @error "Format version mismatch." diskversion formatversion_angfab io
        throw(FormatError())
    end
    nside = read(io, Int64)
    npix = nside2npix(nside)
    mask = HealpixMap{Float64,Healpix.RingOrder}(read!(io, Array{Float64,1}(undef, npix)))
    npixsurvey = read(io, Int64)
    hpix_full2survey = read!(io, Array{Int64,1}(undef, npix))
    ell = read!(io, Vector{Float64}(undef, npixsurvey))
    Y = read!(io, Matrix{Float64}(undef, npixsurvey, npixsurvey))
    Yinv = read!(io, Matrix{Float64}(undef, npixsurvey, npixsurvey))
    return AngularCryoFaB(mask, nside, npixsurvey, hpix_full2survey, ell, Y, Yinv)
end


function AngularCryoFaB(fname::AbstractString)
    io = open(fname, read=true)
    cfb = read(io, AngularCryoFaB)
    close(io)
    return cfb
end


function AngRadCryoFaB(fname::AbstractString)
    io = open(fname, read=true)
    cfb = read(io, AngRadCryoFaB)
    close(io)
    return cfb
end


######## I/O BlockDiagonals
function Base.write(io::IO, bd::BlockDiagonal{T}) where {T}
    n = 0
    allblocks = blocks(bd)
    numblocks = length(allblocks)
    n += write(io, Int64(numblocks))
    for i=1:numblocks
        block = allblocks[i]
        len1, len2 = size(block)
        n += write(io, Int64(len1), Int64(len2))
        n += write(io, Array{T,2}(block))
    end
    return n
end


function Base.read(io::IO, ::Type{BlockDiagonal{T}}) where {T}
    numblocks = read(io, Int64)
    allblocks = Matrix{T}[]
    for i=1:numblocks
        len1 = read(io, Int64)
        len2 = read(io, Int64)
        block = read!(io, Array{T,2}(undef, len1, len2))
        push!(allblocks, block)
    end
    return BlockDiagonal(allblocks)
end


######## I/O AngRadCryoFaB
function Base.write(io::IO, cfb::AngRadCryoFaB)
    n = 0
    n += write(io, Int64(formatversion_angradfab))

    n += write(io, cfb.angfab)
    n += write(io, Int64(length(cfb.Angl)))
    n += write(io, Array{Float64,1}(cfb.Angl))
    n += write(io, Int64(length(cfb.reff)))  # save 'nbins'
    #n += write(io, cfb.AngTrans)    # not needed, since easily reconstructed
    #n += write(io, cfb.AngTransInv) # not needed, since easily reconstructed

    n += write(io, Int64(length(cfb.sortbyl)))  # save length
    n += write(io, Array{Int64,1}(cfb.sortbyl))

    n += write(io, Float64(cfb.rmin))
    n += write(io, Float64(cfb.rmax))
    n += write(io, Int64(length(cfb.rbounds)))  # save length
    n += write(io, Array{Float64,1}(cfb.rbounds))
    n += write(io, Int64(length(cfb.reff)))  # save length
    n += write(io, Array{Float64,1}(cfb.reff))

    n += write(io, Int64(length(cfb.RadK)))  # save length
    n += write(io, Array{Float64,1}(cfb.RadK))
    n += write(io, cfb.RadTrans)
    n += write(io, cfb.RadTransInv)
    return n
end


function Base.read(io::IO, ::Type{AngRadCryoFaB})
    # assumes RingOrder
    diskversion = read(io, Int64)
    if diskversion != formatversion_angradfab
        @error "Format version mismatch." diskversion formatversion_angradfab io
        throw(FormatError())
    end

    angfab = read(io, AngularCryoFaB)
    Angl_len = read(io, Int64)
    Angl = read!(io, Vector{Float64}(undef, Angl_len))
    nbins = read(io, Int64)
    AngTrans = BlockDiagonal(fill(angfab.Trans, nbins))
    AngTransInv = BlockDiagonal(fill(angfab.TransInv, nbins))

    len = read(io, Int64)
    sortbyl = read!(io, Array{Int64,1}(undef, len))

    rmin = read(io, Float64)
    rmax = read(io, Float64)
    len = read(io, Int64)
    rbounds = read!(io, Array{Float64,1}(undef, len))
    len = read(io, Int64)
    reff = read!(io, Array{Float64,1}(undef, len))

    len = read(io, Int64)
    RadK = read!(io, Array{Float64,1}(undef, len))
    RadTrans = read(io, BlockDiagonal{Float64})
    RadTransInv = read(io, BlockDiagonal{Float64})
    return AngRadCryoFaB(angfab, Angl, AngTrans, AngTransInv,
        sortbyl, rmin, rmax, rbounds, reff,
        RadK, RadTrans, RadTransInv)
end



# vim: set sw=4 et sts=4 :
