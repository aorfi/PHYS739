using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function tf_hamiltonian(N,g)
    dim = (2)^N
    H = zeros(dim,dim)
    m_basis = zeros(dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        for SpinIndex in (0:N-1)
            m_basis[ket+1] += 2*((ket>>SpinIndex)&1)-1
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            H[bra+1,ket+1] += -g
            if SpinIndex == N-1
                Si = 2*((ket>>SpinIndex)&1)-1
                Si_next = 2*((ket>>(0))&1)-1
                Bond = Si*Si_next
                H[ket+1,ket+1] += -1 * Bond
                break
            end
            Si = 2*((ket>>SpinIndex)&1)-1
            Si_next = 2*((ket>>(SpinIndex+1))&1)-1
            Bond = Si*Si_next
            Diagonal += Bond
        end
        H[ket+1,ket+1] += -1 * Diagonal
    end
    return H |> sparse, m_basis/N
end


function magn(m_basis, eigvector)
    square(x) = conj(x)*x
    m = abs.(m_basis)'square.(eigvector)
    return m
end


g_beg = range(0.001; stop=0.75, length=10)
g_fine = range(0.76,stop=1.25, length = 50)
g_end = range(1.325; stop=2, length=10)
g = vcat(g_beg,g_fine,g_end )
len= length(g)

N_values = (2:6)
N_max = last(N_values)
energy_all = zeros(length(N_values),len, 3)
mag_all = zeros(length(N_values),len, 3)
for j in (1:length(N_values))
    N = N_values[j]
    println(" Working on N = ",N)
    energy = zeros(len, 3)
    mag = zeros(len,3)
    for i in (1:len)
        H, m_basis = tf_hamiltonian(N,g[i])
        e,v  = eigs(H, nev = 3, which=:SR)
        energy[i,1:3] = e/N
        mag[i,1:3] = magn(m_basis, v)
    end
    energy_all[j,1:len,1:3] = energy 
    mag_all[j,1:len,1:3] = mag
    save_object("Data/energy_all", energy_all)
    save_object("Data/mag_all", mag_all)
end
