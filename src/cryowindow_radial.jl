# Copyright (c) 2021 California Institute of Technology (“Caltech”).
# U.S. Government sponsorship acknowledged.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# - Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# - Neither the name of Caltech nor its operating division, the Jet Propulsion
#   Laboratory, nor the names of its contributors may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
