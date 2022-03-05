using LinearAlgebra


function ising_hamiltonian(N)
    dim = (2)^N
    H = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        for SpinIndex in (0:N-2)
            Si = 2*((ket>>SpinIndex)&1)-1
            Si_next = 2*((ket>>(SpinIndex+1))&1)-1
            Bond = Si*Si_next
            Diagonal += Bond
        H[ket+1,ket+1] = -1 * Diagonal
        end
        SpinIndex = N-1
        Si = 2*((ket>>SpinIndex)&1)-1
        Si_next = 2*((ket>>(0))&1)-1
        Bond = Si*Si_next
        H[ket+1,ket+1] += -1 * Bond
    end
    return H
end


function tf_hamiltonian(N,g)
    dim = (2)^N
    H = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        for SpinIndex in (0:N-1)
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
    return H
end


N=3
g = 1

H = tf_hamiltonian(N,g)
display(H)
e = eigvals(H)
v = eigvecs(H)

