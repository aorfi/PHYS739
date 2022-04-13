include("TFIM.jl")
using Arpack
using LinearAlgebra

function make_rhoA(rho,A_spins)
    dim = length(rho[1,:])
    A_size = 2^A_spins
    B_size = Int(dim/A_size)
    reshaped = reshape(rho,(B_size,A_size,B_size,A_size))
    rhoA = zeros(A_size,A_size)
    for i in (1:A_size)
        for j in (1:A_size)
            current = reshaped[:,j,:,i]
            rhoA[j,i] = tr(current)
        end
    end
    return rhoA
end

