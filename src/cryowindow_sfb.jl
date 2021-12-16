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


#################### 3D cryowindow #######################33

function find_acceptable_tolerance(Hrow, htol)
    Hsort = sort(abs.(Hrow))  # sort ascending
    Hrowtotal = sum(Hsort)
    s = 0.0
    i = length(Hsort)
    while 1 - s / Hrowtotal > htol && i > 1
        s += Hsort[i]
        #@show s,Hrowtotal,s/Hrowtotal,i
        i -= 1
    end
    tolerance = Hsort[i]
    return tolerance
end


# This version writes one large sparse matrix
function cryo_window(cfb; htol=0.01, do_calc=true, lmax=Inf)
    @time ℓ, H = angular_cryo_window(cfb; do_calc, lmax)  # ℓ contains full-sky (integer) ℓ
    @assert all(isfinite.(H))

    ell = Vector{Int}[]  # will contain continuous ℓ
    k = Vector{Float64}[] # will contain continuous k
    fulllen = 0
    @time for l in ℓ
        # find a closely matching L ≈ l, as that is most likely to give relevant k.
        iL = argmin(@. abs(l - cfb.angfab.ell))
        L = cfb.angfab.ell[iL]

        kj = get_kj(cfb, L)
        push!(ell, fill(l, length(kj)))
        push!(k, kj)
        fulllen += length(kj)
    end
    @show fulllen

    cryolen = 0
    I = Int[]
    J = Int[]
    V = Float64[]
    @time for i=1:size(H,1)
        #L = cfb.angfab.ell[i]  # cryo-ℓ
        #kN = get_kN(cfb, L)  # cryo-k
        kN_len = length(cfb.reff)
        cryolen += kN_len

        (!do_calc) && continue
        (cfb.angfab.ell[i] > lmax) && continue

        tolerance = find_acceptable_tolerance(H[i,:], htol)
        @show i,cfb.angfab.ell[i],sum(abs.(H[i,:])),tolerance
        flen = 0
        @time for j=1:size(H,2)
            flen += length(k[j])
            if abs(H[i,j]) >= tolerance
                ibase = cryolen - kN_len
                jbase = flen - length(k[j])
                Δkj = k[j][2] - k[j][1]
                b̃ll = radial_cryo_window2(cfb, i, ℓ[j], k[j])
                b̃ll2 = sparsify(abs2.(b̃ll), rtol=htol)
                bI, bJ, bV = findnz(b̃ll2)
                @assert all(isfinite.(bV))
                append!(I, ibase .+ bI)
                append!(J, jbase .+ bJ)
                append!(V, H[i,j] * Δkj .* bV)
            end
        end
    end
    W = sparse(I, J, V, cryolen, fulllen)
    k = [kj for kv in k for kj in kv]
    ell = [lj for lv in ell for lj in lv]
    return k, ell, W
end


function sparsify(m; rtol=0.01)
    for i=1:size(m,1)
        row = m[i,:]
        tolerance = find_acceptable_tolerance(row, rtol)
        s = @. abs(row) < tolerance
        m[i,s] .= 0
    end
    return sparse(m)
end



# vim: set sw=4 et sts=4 :
