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
