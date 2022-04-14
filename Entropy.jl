include("TFIM.jl")
using Arpack
using LinearAlgebra

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
    A_size = 2^A_spins
    B_size = Int(dim/A_size)
    reshaped = reshape(v,(B_size,A_size))
    U,s,Vt= svd(reshaped)
    lambda  = s.^2
    S = -lambda'*log.(lambda+0.0001*ones(length(lambda)))
    return S
end


# g = 0 
# N = 2
# H, m_basis = tf_hamiltonian(N,g)
# e,v  = eigs(H, nev = 1, which=:SR)
# print(v[1:2^N])
psi = [1,1,0,0,1,1,0,0]
dim = length(psi)
A_spins = 1
A_size = 2^A_spins
B_size = Int(dim/A_size)
reshaped = reshape(psi,(A_size,B_size))
# print(reshaped)
U,s,Vt= svd(reshaped)
lambda  = s.^2
println("rho: ")
display(psi*psi')
println("rhoA: ")
display(U*Diagonal(lambda)*transpose(U))

print(lambda)
print(log.(lambda+0.0001*ones(length(lambda))))
S = -lambda'*log.(lambda+0.0001*ones(length(lambda)))
println("Entropy: ", S)

# rho = psi*psi'
# display(make_rhoA(rho,A_spins))
