
######################### Radial Functions #########################3

function radial_greens_function(L, r1, r2)
    if r1 <= r2
        return r2 / (2*L+1) * (r1/r2)^L
    else
        return r2 / (2*L+1) * (r2/r1)^(L+1)
    end
end


function get_radial_greens_matrix(L, ri, Δr)
    G = fill(NaN, length(ri), length(ri))
    for j=1:length(ri), i=j:length(ri)
        G[i,j] = G[j,i] = (ri[i] / ri[j]) * radial_greens_function(L, ri[i], ri[j])
    end
    G .*= Δr
    return Hermitian(G)
end


function get_radial_cryofunks(L, ri, Δr)
    B = diagm(@. √Δr * ri)
    G = get_radial_greens_matrix(L, ri, Δr)
    e = eigen(G)
    k = @. 1 / √e.values
    p = sortperm(k)
    k = k[p]
    Z = e.vectors[:,p]
    g = inv(B) * Z
    ginv = Z' * B
    return k, g, ginv
end



# vim: set sw=4 et sts=4 :
