include("TFIM.jl")
include("XXZ.jl")
using LinearAlgebra
using Arpack
using SparseArrays
# using PyPlot

# This isn't needed
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

function entropy(v,A_spins)
    dim = length(v)
    v = reshape(v,(dim,1))
    A_size = Int(2^A_spins)
    B_size = Int(dim/A_size)
    reshaped = reshape(v,(B_size,A_size))
    U,s,Vt= svd(reshaped)
    lambda  = s.^2
    S = -lambda'*log.(lambda+0.0001*ones(length(lambda)))
    return S
end


# # TFIM 
# N=15
# g=2
# H, m_basis = tf_hamiltonian(N,g)
# e,v  = eigs(H, nev = 1, which=:SR)
# x_all = range(1,N-1)
# S_all = zeros(N-1)
# for i in (1:length(x_all))
#     S_all[i] = entropy(v,x_all[i])
# end
# save_object("Data/Q3/TFIM_g2", S_all)

# XXZ
# N=16
# Delta = 1
# J_perp = 2
# H = XXZ_sz0(N,Delta,J_perp)
# e,v  = eigs(H, nev = 1, which=:SR)
# v_full = convert_full(N,v)
# x_all = range(1,N-1)
# S_all = zeros(N-1)
# for i in (1:length(x_all))
#     println("Working on ",i)
#     S_all[i] = entropy(v_full,x_all[i])
# end
# save_object("Data/Q3/XXY2", S_all)



