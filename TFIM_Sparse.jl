using LinearAlgebra
using SparseArrays
using Arpack

function tf_hamiltonian_sparse(N,g)
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

