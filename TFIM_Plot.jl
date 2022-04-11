include("TFIM.jl")
using LaTeXStrings


# Plot energy/magnetization  of low lying states
len = 100
g = range(0; stop=2, length=len)

N = 3
energy = zeros(len, 2^N)
mag = zeros(len)


for i in (2:len)
    H, m_basis = tf_hamiltonian(N,g[i])
    e = eigvals(H)
    energy[i,1:2^N] = e/N
    v = eigvecs(H)
    # print(v[1:2^N])
    # print(square.(v[1:2^N]))
    m = magn(m_basis, v[1:2^N])
    mag[i] = abs(m)/N 
end

# plt.plot(g, energy[1:len,1])
# plt.plot(g, energy[1:len,2])
# plt.plot(g, energy[1:len,3])
# plt.title(string("Energy of Low-Lying States N = ", N))
# plt.ylabel(L"$\epsilon = E/N$")
# plt.xlabel(L"$h/J$")
# plt.show()
# plt.savefig("/Figures/1aEnergy.png")


plt.scatter(g, mag[1:len])
# plt.scatter(g, mag[1:len,2])
# plt.scatter(g, mag[1:len,3])
plt.title(string("Energy of Low-Lying States N = ", N))
plt.ylabel(L"$|m|/N$")
plt.xlabel(L"$h/J$")
# plt.savefig("1aEnergy.png")
plt.show()


