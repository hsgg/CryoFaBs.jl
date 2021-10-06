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


##################### Radial cryowindow ################################

function calc_Rb̃llnn(ℓ, rlo, rhi, keff)
    I, E = quadgk(r->sphericalbesselj(ℓ, r*keff), rlo, rhi)
    return I * keff * √(2/π) / (rhi - rlo)
end


function calc_Rb̃llnn(ℓ, rlo, rhi, keff, jl_spl)
    xlo = rlo * keff
    xhi = rhi * keff
    I = √(2/π) / (rhi - rlo) * integrate(jl_spl, xlo, xhi)
    if !isfinite(I) # || !(I ≈ I0)
        I0 = calc_Rb̃llnn(ℓ, rlo, rhi, keff)
        @error "Rb̃llnn error" ℓ rlo rhi keff xlo xhi jl_spl(xlo) jl_spl(xhi) I I0 I≈I0
        error("Rb̃llnn error")
    end
    return I
end


function get_kN(cfb, L)
    kk, ll = CryoFaBs.idx2mode(cfb, 1:size(cfb,2))
    kN = kk[ll .== L]
    # There are ~2L+1 sets of keff that are equal, choose one. (If we don't know which one it is, then we don't care which one it is.)
    rlo = cfb.rbounds[1:end-1]
    kN = kN[1:length(rlo)]
    return kN
end


function get_kj(cfb, L; nmult=10)
    kN = get_kN(cfb, L)
    #kj = kN
    Δr = (cfb.rmax - cfb.rmin) / length(cfb.reff)
    kmax = 2 / Δr
    #kj = kN
    kj = range(0.0, π/2*maximum(kN), length=nmult*length(kN))
    #kj = range(minimum(kN), maximum(kN), length=4*length(kN))
    return kj
end


function radial_cryo_window(cfb::AngRadCryoFaB, iL::Integer, ℓ::Integer, k::AbstractVector{Float64})
    rlo = cfb.rbounds[1:end-1]
    rhi = cfb.rbounds[2:end]
    Rb̃ll = calc_Rb̃llnn.(ℓ, rlo, rhi, k')
    Rinv = blocks(cfb.RadTransInv)[iL]
    b̃ll = Rinv * Rb̃ll
    return b̃ll::Matrix{Float64}
end


function radial_cryo_window2(cfb::AngRadCryoFaB, iL::Integer, ℓ::Integer, k::AbstractVector{Float64})
    rlo = cfb.rbounds[1:end-1]
    rhi = cfb.rbounds[2:end]
    xmin = k[1]*rlo[1]
    xmax = k[end]*rhi[end]
    nsamps = (xmax - xmin) * 20 / π
    x = range(xmin, xmax, length=ceil(Int, nsamps))
    jl_spl = Spline1D(x, sphericalbesselj.(ℓ, x))
    Rb̃ll = calc_Rb̃llnn.(ℓ, rlo, rhi, k', jl_spl)
    Rinv = blocks(cfb.RadTransInv)[iL]
    b̃ll = Rinv * Rb̃ll
    return b̃ll::Matrix{Float64}
end


# vim: set sw=4 et sts=4 :
