include("TFIM.jl")
using LaTeXStrings
using Arpack
using PyPlot
using JLD2

# len = 50
# g = range(0.01; stop=2, length=len)
# N_values = (2:12)
# N_max = last(N_values)
# energy_all = zeros(length(N_values),len, 3)
# mag_all = zeros(length(N_values),len, 3)
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     energy = zeros(len, 3)
#     mag = zeros(len,3)
#     for i in (1:len)
#         H, m_basis = tf_hamiltonian(N,g[i])
#         e,v  = eigs(H, nev = 3, which=:SR)
#         energy[i,1:3] = e/N
#         mag[i,1:3] = magn(m_basis, v)
#     end
#     energy_all[j,1:len,1:3] = energy 
#     mag_all[j,1:len,1:3] = mag
#     save_object("Data/energy_all", energy_all)
#     save_object("Data/mag_all", mag_all)
# end

# # Plot energy/magnetization  of low lying states 
# len = 50
# g = range(0.01; stop=2, length=len)
energy_all = load_object("Data/energy_all")
display(energy_all[:,1,1])
# N = size(energy_all,1) +1
# energy = energy_all[N-1,:,:]
# # Energy vs g plot
# plt.plot(g, energy[1:len,1], label = "Ground State")
# plt.plot(g, energy[1:len,2],label = "First Excited State")
# plt.plot(g, energy[1:len,3], label = "Second Excited State")
# plt.title(string("Energy of Low-Lying States N = ", N))
# plt.ylabel(L"$\epsilon = E/N$")
# plt.xlabel(L"$h/J$")
# plt.legend()
# plt.show()
# # plt.savefig("Figures/1aEnergy.png")

# # Mag vs g plot
# len = 50
# g = range(0.01; stop=2, length=len)
# mag_all = load_object("Data/mag_all")
# N = size(mag_all,1) +1
# mag = mag_all[N-1,:,:]
# plt.plot(g, mag[1:len,1], label = "Ground State")
# plt.plot(g, mag[1:len,2], label = "First Excited State")
# plt.plot(g, mag[1:len,3], label = "Second Excited State")
# plt.title(string("Energy of Low-Lying States N = ", N))
# plt.ylabel(L"$|m|/N$")
# plt.xlabel(L"$h/J$")
# plt.legend()
# # plt.savefig("Figures/1aMag.png")
# plt.show()

# # Ground-state energy for different N
# len = 20
# g = range(0.01; stop=2, length=len)
# N_values = (2:8)
# N_max = last(N_values)
# energy_all = zeros(length(N_values),len, 2^N_max)
# mag_all = zeros(length(N_values),len, 2^N_max)
# for j in (1:length(N_values))
#     N = N_values[j]
#     energy = zeros(len, 2^N)
#     mag = zeros(len,2^N)
#     for i in (1:len)
#         H, m_basis = tf_hamiltonian(N,g[i])
#         e = eigvals(H)
#         energy[i,1:2^N] = e/N
#         v = eigvecs(H)
#         mag[i,1:2^N] = magn(m_basis, v)
#     end
#     energy_all[j,1:len,1:2^N] = energy 
#     mag_all[j,1:len,1:2^N] = mag
# end
# for i in (1:length(N_values))
#     N = N_values[i]
#     plt.plot(g, mag_all[i,1:len,1], label = string("N= ", N))
# end
# plt.title("Energy of Ground State ")
# plt.ylabel(L"$|m|/N$")
# plt.xlabel(L"$h/J$")
# plt.legend()
# plt.savefig("Figures/1aMagScaling.png")
# plt.show()

# # Gap near critical point
# len = 3
# g = range(0.95; stop=1.05, length=len)
# N_values = (2:5)
# N_max = last(N_values)
# gap_all = zeros(length(N_values),len)
# for j in (1:length(N_values))
#     N = N_values[j]
#     print(" Current N Value: ", N)
#     gap = zeros(len)
#     mag = zeros(len,2^N)
#     for i in (1:len)
#         print(" Making H ")
#         H, m_basis = tf_hamiltonian(N,g[i])
#         print(" Eigenvalues ")
#         e = eigs(H, nev = 2, which=:SR)
#         gap[i] = (e[1][1]-e[1][2])/N
#     end
#     gap_all[j,1:len] = gap
# end


# for i in (1:len)
#     g_value = g[i]
#     plt.plot(1 ./ (N_values), gap_all[1:length(N_values),i,1], label = string("g = ", g_value))
# end
# plt.title("Energy Gap of Ground State ")
# plt.ylabel("gap")
# plt.xlabel(L"$1/N$")
# plt.legend()
# # plt.savefig("Figures/1cGapScaling.png")
# plt.show()
