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


# Here we collect the angular cryo-window functions.

function sphericalharmonicsy(l, m, θ, ϕ)
    #return scipy.special.sph_harm(m, l, ϕ, θ)  # scipy has ℓ,m and θ,ϕ reversed
    rY1 = sphevaluate(θ, ϕ, l, abs(m))
    rY2 = sphevaluate(θ, ϕ, l, -abs(m))
    if m < 0
        return (rY1 - im*rY2) / √2
    elseif m == 0
        return rY1 + 0*im
    else
        return (-1)^m * (rY1 + im*rY2) / √2
    end
end

function realsphericalharmonicsy(l, m, θ, ϕ)
    #Ylm = sphericalharmonicsy(l, m, θ, ϕ)
    #if m < 0
    #    return √2 * (-1)^m * imag(Ylm)
    #elseif m == 0
    #    return real(Ylm)
    #else
    #    return √2 * (-1)^m * real(Ylm)
    #end
    return sphevaluate.(θ, ϕ, l, m)
end


function get_wl_Y_j_lm(bmask, nside_hires=npix2nside(length(bmask)); lmax=Inf)
    npix = length(bmask)
    lmsize = npix
    pix = collect((1:npix)[bmask])
    nside = npix2nside(npix)
    pixinv = fill(0, npix)
    for i=1:length(pix)
        pixinv[pix[i]] = i
    end

    pw = Healpix.pixwin(nside)
    pixwin = Spline1D(0:length(pw)-1, pw)

    for i in pix
        θ, ϕ = Healpix.pix2angRing(nside, [i])
        n = Healpix.ang2pixRing(Resolution(nside), θ[1], ϕ[1])
        @assert n==i
    end

    bmask_hires = udgrade(HealpixMap{Bool,Healpix.RingOrder}(bmask), nside_hires)
    npix_hires = nside2npix(nside_hires)
    pix_hires = collect((1:npix_hires)[bmask_hires])
    θ, ϕ = Healpix.pix2angRing(nside_hires, pix_hires)
    pixmap = [Healpix.ang2pixRing(Resolution(nside), θ[i], ϕ[i]) for i=1:length(pix_hires)]

    wl_Yjlm = fill(0.0, length(pix), lmsize)
    l = 0
    m = 0
    for lm=1:lmsize
        (l > lmax) && continue
        Ylm = CryoFaBs.realsphericalharmonicsy(l, m, θ, ϕ)
        for i=1:length(pix_hires)
            n = pixinv[pixmap[i]]
            n==0 && @show i,pixmap[i],n
            wl_Yjlm[n,lm] += Ylm[i]
            if !isfinite(wl_Yjlm[n,lm])
                @error "Non-finite value in w_Yjlm! Setting to zero." l m n lm wl_Yjlm[n,lm] i Ylm[i] pixmap[i] θ[i] ϕ[i] CryoFaBs.sphericalharmonicsy(l, m, θ[i], ϕ[i])
                wl_Yjlm[n,lm] = 0
            end
        end
        wl_Yjlm[:,lm] .*= npix / npix_hires

        if m < 0
            m = -m
        else # m >= 0
            m = -m - 1
        end

        if m < -l
            l += 1
            m = 0
        end
    end
    return wl_Yjlm
end


function fill_l(n)
    ℓ = Int[]
    l = 0
    m = 0
    for i=1:n
        push!(ℓ, l)
        m += 1
        if m > l
            l += 1
            m = -l
        end
    end
    return ℓ
end


function angular_cryo_window_sum_m(w̃)
    ell = fill_l(size(w̃,2))
    ell_uniq = unique(ell)
    ww = fill(0.0, size(w̃,1), length(ell_uniq))
    for j=1:length(ell_uniq)
        m = (ell .== ell_uniq[j])
        for i=1:size(w̃,1)
            ww[i,j] = sum(abs2.(w̃[i,m]))
        end
    end
    return ell_uniq, ww
end


function angular_cryo_window(cfb::AngularCryoFaB; do_calc=true, lmax=Inf)
    nside = cfb.nside
    npix = nside2npix(nside)
    pwspl = Spline1D(0:4*nside, Healpix.pixwin(nside))  # Don't really need this to be a spline.
    pwlm = pwspl.(fill_l(npix))
    bmask = @. cfb.mask > 0
    if do_calc
        wjlm = pwlm' .* get_wl_Y_j_lm(bmask, nside; lmax=lmax)
        w̃ = cfb.TransInv * wjlm
    else
        w̃ = zeros(size(cfb,2), length(pwlm))
    end
    unique_ell, ww = angular_cryo_window_sum_m(w̃)
    return unique_ell, ww
end

function angular_cryo_window(cfb::AngRadCryoFaB; do_calc=true, lmax=Inf)
    return angular_cryo_window(cfb.angfab; do_calc=do_calc, lmax=lmax)
end




# vim: set sw=4 et sts=4 :
