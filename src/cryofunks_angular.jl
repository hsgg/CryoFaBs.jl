
############### Angular cryofunks ########################

function get_angular_cryofunks(mask; G0=nothing, nside_hires=npix2nside(length(mask)))
    if G0 == nothing
        @time G0 = get_green(npix2nside(length(mask)), nside_hires=nside_hires)
    end
    G = Hermitian(collect(G0[mask.>0,mask.>0]))

    #fsky = sum(mask .> 0) / length(mask)
    #G = Hermitian(G + I / (4*π*fsky*length(mask)))

    G = monopolize_green(G)

    B = get_B(mask)
    @time ℓ, Z = get_cryobasis(G)

    ℓ[1] = 0  # monopole
    Z[:,1] .= mean(Z[:,1])  # monopole

    ℓ, Z = sort_cryobasis(ℓ, Z)
    Y = B \ Z
    Yinv = Z' * B

    return ℓ, Y, Yinv
end


function get_cryobasis(G)
    f = eigen(G)
    λ = f.values  # already sorted smallest to largest
    Z = f.vectors

    @show λ[1:4] λ[end-4:end]
    epsmax = eps(maximum(abs.(λ)))

    if maximum(abs.(imag(λ))) <= epsmax
        # assume this is numerical noise, and imag really vanishes, as for all
        # Hermitian matrices.
        λ = real(λ)
    end

    if abs(λ[1]) <= epsmax
        # assume this is numerical noise, and it should really vanish, as for
        # the monopole.
        λ[1] = 0
    end

    @show λ[1] λ[end]
    #λ = complex(λ)

    ℓ = @. sqrt(1 / λ + 1/4) - 1/2
    #sel = @. imag(ℓ) != 0
    #ℓ[sel] .= abs.(ℓ[sel]) .+ maximum(real(ℓ))
    #ℓ = real(ℓ)

    return ℓ, Z
end


# simple form to fill in the Green's matrix
function get_green1(mask; kwargs...)  # 'kwargs' are ignored
    bmask = collect(mask .> 0)

    nside = npix2nside(length(bmask))
    npix = nside2npix(nside)
    Ωᵢ = 4 * π / npix

    Gii_same = Ωᵢ / (4 * π) * (1 - log(Ωᵢ/(2*π)))
    θ = fill(NaN, npix)
    ϕ = fill(NaN, npix)
    for i=1:npix
        θ[i], ϕ[i] = pix2angRing(Resolution(nside), i)
    end
    sinθ = sin.(θ)

    @show npix^2 * 8/1024^3
    Gij = fill(NaN, npix, npix)
    for j=1:npix, i=j:npix
        (bmask[i] && bmask[j]) || continue
        if i == j
            Gij[i,j] = Gii_same
        else
            Δθ = θ[j] - θ[i]
            Δϕ = ϕ[j] - ϕ[i]
            s²ρ2 = sin(Δθ/2)^2 + sinθ[i]*sinθ[j]*sin(Δϕ/2)^2
            Gij[i,j] = - log(2 * s²ρ2) * Ωᵢ / (4 * π)
            Gij[j,i] = Gij[i,j]'
        end
    end

    return Hermitian(Gij)
end


# integrate over each pixel
function get_green3(mask; nside_hires=1*npix2nside(length(mask)))
    T = Float64
    nside = npix2nside(length(mask))
    npix = nside2npix(nside)
    npix_hires = nside2npix(nside_hires)

    bmask = collect(mask .> 0)

    Ωᵢhires = 4 * π / npix_hires
    Gii_same = Ωᵢhires / (4 * π) * (1 - log(Ωᵢhires/(2*π)))
    θ = fill(NaN, npix_hires)
    ϕ = fill(NaN, npix_hires)
    for i=1:npix_hires
        θ[i], ϕ[i] = pix2angRing(Resolution(nside_hires), i)
    end
    sinθ = sin.(θ)
    cosθ = cos.(θ)

    pixmap = [ang2pix(mask, θ[i], ϕ[i]) for i=1:npix_hires]

    @show nside nside_hires
    @show npix npix_hires
    @show npix^2 * 8/1024^3
    @show npix_hires^2 * 8/1024^3
    Gmn = fill(T(0), npix, npix)
    for j=1:npix_hires
        n = pixmap[j]
        bmask[n] || continue
        for i=1:npix_hires
            m = pixmap[i]
            bmask[m] || continue
            if i == j
                Gmn[m,n] += Gii_same
            else
                Δθ = θ[j] - θ[i]
                Δϕ = ϕ[j] - ϕ[i]

                # arcsin formula
                s²ρ2 = sin(Δθ/2)^2 + sinθ[i]*sinθ[j]*sin(Δϕ/2)^2
                Gij = - log(2 * s²ρ2) * Ωᵢhires / (4 * π)

                # arctan formula
                r² = (sinθ[j] * sin(Δϕ))^2 + (sinθ[i] * cosθ[j] - cosθ[i] * sinθ[j] * cos(Δϕ))^2
                d = cosθ[i] * cosθ[j] + sinθ[i] * sinθ[j] * cos(Δϕ)
                ρ = atan(sqrt(r²), d)  # Don't do a naive atan! atan(x) looses info.
                Gij2 = - Ωᵢhires / (4*π) * log(2 * sin(ρ/2)^2)

                #@show Gij Gij2 θ[i],ϕ[j] Δθ,Δϕ 2asin(√s²ρ2)-π/2 ρ+π/2
                #@show tan(2asin(√s²ρ2)) tan(ρ)
                #tanρ = sqrt(r²) / d
                #@show s²ρ2 tanρ^2 / (1 + sqrt(1 + tanρ^2)) / sqrt(1 + tanρ^2) / 2
                false && if abs(Gij - Gij2) > 1e-10 * (abs(Gij) + abs(Gij2)) / 2
                    θ1 = big(θ[i])
                    θ2 = big(θ[j])
                    ϕ1 = big(ϕ[i])
                    ϕ2 = big(ϕ[j])
                    Δθ = θ2 - θ1
                    Δϕ = ϕ2 - ϕ1

                    s²ρ2 = sin(Δθ/2)^2 + sin(θ1)*sin(θ2)*sin(Δϕ/2)^2
                    Gij = - log(2 * s²ρ2) * Ωᵢhires / (4 * π)

                    # arctan formula
                    r² = (sin(θ2) * sin(Δϕ))^2 + (sin(θ1) * cos(θ2) - cos(θ1) * sin(θ2) * cos(Δϕ))^2
                    d = cos(θ1) * cos(θ2) + sin(θ1) * sin(θ2) * cos(Δϕ)
                    ρ = atan(sqrt(r²) / d)
                    @error "Mismatch!" θ1 ϕ1 θ2 ϕ2 Δθ Δϕ atan(tan(2asin(√s²ρ2))) atan(tan(ρ))
                    error("Wrong!")
                end

                Gmn[m,n] += Gij2
            end
        end
    end
    Ωᵢ = 4 * π / npix
    Gmn .*= Ωᵢhires / Ωᵢ
    return Hermitian(Gmn)
end


function calc_green(mask; nside_hires=npix2nside(length(mask)))
    bmask = collect(mask .> 0)
    mask = HealpixMap{eltype(mask), Healpix.RingOrder}(mask)
    return Hermitian(collect(get_green3(mask, nside_hires=nside_hires)[bmask,bmask]))
end


function get_green(nside::Int, fname=nothing; nside_hires=nside)
    if fname == nothing
        fname = "$cachedir_greensmatrix/Gang_nside$(nside)_$(nside_hires).bin"
    end
    npix = nside2npix(nside)
    if !isfile(fname)
        @timeit "Angular Green's matrix (nside,nside_hires)=($nside,$nside_hires)" G = calc_green(fill(1.0, npix), nside_hires=nside_hires)
        write(fname, G)
    end
    return Hermitian(read!(fname, fill(NaN, npix, npix)))
end


function monopolize_green(G)
    u = fill(1.0, size(G,1))
    u = u / sqrt(u'*u)
    P = Hermitian(I - u*u')
    G′ = P * G * P'
    if typeof(G) <: Hermitian
        # Note: Here we set the type, because that will significantly speed up
        # the eigenvector decomposition. It really should be considered a bug
        # in Julia that P*G*P' does not automatically result in the best type.
        G′ = Hermitian(G′)
        G′.data .= G′
    end
    @show u'*u typeof(u*u')
    @show typeof(G′)
    return G′
end


function sort_cryobasis(ℓ, Z)
    perm = sortperm(real(ℓ))
    ℓ .= ℓ[perm]
    Z[:,:] .= Z[:,perm]
    return collect(ℓ), collect(Z)
end


function get_B(mask)
    nside = npix2nside(length(mask))
    npix = size(mask,1)
    Ωᵢ = 4 * π / npix
    B = √Ωᵢ* I
    return B
end

######################### Angular Functions, integrate directly ######

function Healpix.pix2angRing(nside::Integer, idxs::AbstractArray)
    reso = Resolution(nside)
    θ = fill(NaN, length(idxs))
    ϕ = fill(NaN, length(idxs))
    for i=1:length(idxs)
        θ[i], ϕ[i] = pix2angRing(reso, idxs[i])
    end
    return θ, ϕ
end


function get_angular_continuous_cryofunks(mask; nside_hires=nside_hires)
    npix = length(mask)
    npix_hires = nside2npix(nside_hires)
    nside = npix2nside(npix)
    Ω = 4*π / npix

    ℓ, Y = calc_pixwin(nside, nside_hires=nside_hires)

    Yinv = Y' * Ω  # one of the eigenvalues vanishes, but this is all we need.

    return ℓ, Y, Yinv
end



# vim: set sw=4 et sts=4 :
