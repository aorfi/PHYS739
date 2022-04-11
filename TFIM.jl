using LinearAlgebra
using PyPlot

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
    return H, m_basis/N
end


function magn(m_basis, eigvector)
    m = Float64(0)
    for i in (1:length(m_basis))
        m += conj(eigvector[i])*abs(m_basis[i])*eigvector[i]
    end
    return m
end



# N_values = range(2,6)
# maxN = last(N_values)
# energy = zeros(length(N_values),len, 2^maxN)
# mag = zeros(length(N_values),len,2^maxN)

# for N in N_values
#     dim = (2)^N
#     for i in (1:len)
#         H, m_basis = tf_hamiltonian(N,g[i])
#         e = eigvals(H)
#         energy[N-1,i,1:length(e)] = e/N
#         v = eigvecs(H)
#         m = magn(m_basis,e,v)
#         mag[N-1,i,1:length(e)] = m/N
#     end
# end


# N = 5
# plt.scatter(g, broadcast(abs,mag[N-1,1:len,1]))
# plt.scatter(g, broadcast(abs,mag[N-1,1:len,2]))
# plt.scatter(g, broadcast(abs,mag[N-1,1:len,3]))
# # plt.scatter(g, broadcast(abs,mag[N-1,1:len,4]))
# # plt.scatter(g, broadcast(abs,mag[N-1,1:len,5]))

# # plt.plot(g, energy[3,1:len,1])
# # plt.plot(g, energy[3,1:len,2])
# # plt.plot(g, energy[3,1:len,3])
# # plt.plot(g, energy[3,1:len,4])

# plt.title(string("Ground State N = ", N))
# plt.ylabel("mag")
# plt.xlabel("h/J")
# plt.show()