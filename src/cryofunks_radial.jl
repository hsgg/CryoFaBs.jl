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
