include("XXZ.jl")
include("Entropy.jl")
using SparseArrays
using LinearAlgebra
using Arpack
# using PyPlot
using JLD2
using Random
using Statistics

rng = MersenneTwister(1234)

function rand_ham(N,Delta,J_perp,W)
    dim = 2^N
    H = zeros(dim,dim)
    rng = MersenneTwister(1234)
    random_hi_mag = rand!(rng, zeros(N))
    random_hi_sign = 2 .*bitrand(rng, N).-1
    random_his = W.*random_hi_mag.*random_hi_sign
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = 0
        random_field = 0
        for SpinIndex in (0:N-1)
            hi = random_his[SpinIndex+1]
            random_field += hi*(2*((ket>>SpinIndex)&1)-1)
            if SpinIndex == N-1
                bit_last = Int(2)^(SpinIndex)
                bit_first = Int(2)^(0)
                flipbit_last = ket ⊻ bit_last
                flipbit_first  = flipbit_last ⊻ bit_first
                Sy_flip_last = (2*((flipbit_first>>SpinIndex)&1)-1)im
                Sy_flip_first = (2*((flipbit_first>>(0))&1)-1)im
                Sy = Sy_flip_last*Sy_flip_first
                H[flipbit_first+1,ket+1] += -(J_perp/4)*(1 + Sy)
                Sz_last = 2*((ket>>SpinIndex)&1)-1
                Sz_first = 2*((ket>>(0))&1)-1
                Bondz = (Delta/4)*Sz_last*Sz_first
                Diagonal += Bondz
                break
            end
            bit = Int(2)^(SpinIndex)
            nn_bit = Int(2)^(SpinIndex+1)
            flipbit = ket ⊻ bit
            Sy_flip1 = (2*((flipbit>>SpinIndex)&1)-1)im
            flipnn  = flipbit ⊻ nn_bit
            Sy_flip2 = (2*((flipnn>>(SpinIndex+1))&1)-1)im
            Sy = Sy_flip1*Sy_flip2
            H[flipnn+1,ket+1] += -(J_perp/4)*(1 + Sy)
            Szi = 2*((ket>>SpinIndex)&1)-1
            Szi_next = 2*((ket>>(SpinIndex+1))&1)-1
            Bondz = (Delta/4)*Szi*Szi_next
            Diagonal += Bondz
        end
        H[ket+1,ket+1] += (-1 * (Diagonal)) - random_field
    end
    H_sparse = H |> sparse
    return H_sparse
end

function rand_ham_sz0(N,Delta,J_perp,W)
    dict = sz0_dict(N)
    dim = length(dict)
    H = zeros(dim,dim)
    random_his = W.-2*W.*rand!(rng, zeros(N))
    for ket in keys(dict)
        ket_binary = bitstring(ket)
        Diagonal = 0
        random_field = 0
        for SpinIndex in (0:N-1)
            hi = random_his[SpinIndex+1]
            random_field += hi*(2*((ket>>SpinIndex)&1)-1)
            if SpinIndex == N-1
                Sz_last = 2*((ket>>SpinIndex)&1)-1
                Sz_first = 2*((ket>>(0))&1)-1
                Bondz = -(Delta/4)*Sz_last*Sz_first
                Diagonal += Bondz
                bit_last = Int(2)^(SpinIndex)
                bit_first = Int(2)^(0)
                Sx_i = ket ⊻ bit_last
                Sx_j  = Sx_i ⊻ bit_first
                Sy_flip_last = (2*((Sx_j>>SpinIndex)&1)-1)im
                Sy_flip_first = (2*((Sx_j>>(0))&1)-1)im
                Sy = Sy_flip_last*Sy_flip_first
                value = -(J_perp/4)*(1 + Sy)
                if value != 0.0
                    index_original = dict[ket]
                    index_flipped = dict[Sx_j]
                    H[index_flipped,index_original] += value
                end
                break
            end
            bit = Int(2)^(SpinIndex)
            nn_bit = Int(2)^(SpinIndex+1)
            Sx_i = ket ⊻ bit
            Sy_flip1 = (2*((Sx_i>>SpinIndex)&1)-1)im
            Sx_j  = Sx_i ⊻ nn_bit
            Sy_flip2 = (2*((Sx_j>>(SpinIndex+1))&1)-1)im
            Sy = Sy_flip1*Sy_flip2
            value = -(J_perp/4)*(1 + Sy)
            if value != 0.0
                index_original = dict[ket]
                index_flipped = dict[Sx_j]
                H[index_flipped,index_original] += value
            end
            Szi = 2*((ket>>SpinIndex)&1)-1
            Szi_next = 2*((ket>>(SpinIndex+1))&1)-1
            Bondz = -(Delta/4)*Szi*Szi_next
            Diagonal += Bondz
        end
        index_ket = dict[ket]
        H[index_ket,index_ket]+= Diagonal - random_field
    end
    H_sparse = H |> sparse
    return H_sparse
end


Delta = 1
J_perp = 1
W = 0.5
runs = 100
N_values = range(8,16,step = 2)
# avg_S = zeros(length(N_values))
# std_S = zeros(length(N_values))
# for i in (1:length(N_values))
#     N = N_values[i]
#     println("Working on N: ", N)
#     S_all = zeros(runs)
#     for i in (1:runs)
#         x = Int(N/2)
#         println("Iteration number: ", i)
#         H = rand_ham_sz0(N,Delta,J_perp,W)
#         e,v  = eigs(H, nev = 1, which=:SR)
#         v_full = convert_full(N,v)
#         # e,v  = eigs(H, nev = x, which=:SR)
#         # v_ex = v[:,x]
#         # v_full = convert_full(N,v_ex)
#         S_all[i] = entropy(v_full,x)
#     end
#     avg_S[i] = mean(S_all)
#     std_S[i] = std(S_all)
# end
# save_object("Data/Q4/Sgs0.5", avg_S)
# save_object("Data/Q4/STDgs0.5", std_S)

avg_S = zeros(length(N_values))
std_S = zeros(length(N_values))
for i in (1:length(N_values))
    N = N_values[i]
    println("Working on N: ", N)
    S_all = zeros(runs)
    for i in (1:runs)
        x = Int(N/2)
        println("Iteration number: ", i)
        H = rand_ham_sz0(N,Delta,J_perp,W)
        # e,v  = eigs(H, nev = 1, which=:SR)
        # v_full = convert_full(N,v)
        ex_pos = trunc(Int, 0.7*length(H[1,:]))
        e,v  = eigs(H, nev = ex_pos, which=:SR)
        v_ex = v[:,ex_pos]
        v_full = convert_full(N,v_ex)
        S_all[i] = entropy(v_full,x)
    end
    avg_S[i] = mean(S_all)
    std_S[i] = std(S_all)
end
save_object("Data/Q4/Ses0.5", avg_S)
save_object("Data/Q4/STDes0.5", std_S)

# W = 9
# avg_S = zeros(length(N_values))
# std_S = zeros(length(N_values))
# for i in (1:length(N_values))
#     N = N_values[i]
#     println("Working on N: ", N)
#     S_all = zeros(runs)
#     for i in (1:runs)
#         x = Int(N/2)
#         println("Iteration number: ", i)
#         H = rand_ham_sz0(N,Delta,J_perp,W)
#         e,v  = eigs(H, nev = 1, which=:SR)
#         v_full = convert_full(N,v)
#         # e,v  = eigs(H, nev = x, which=:SR)
#         # v_ex = v[:,x]
#         # v_full = convert_full(N,v_ex)
#         S_all[i] = entropy(v_full,x)
#     end
#     avg_S[i] = mean(S_all)
#     std_S[i] = std(S_all)
# end
# save_object("Data/Q4/Sgs9", avg_S)
# save_object("Data/Q4/STDgs9", std_S)

W = 9
avg_S = zeros(length(N_values))
std_S = zeros(length(N_values))
for i in (1:length(N_values))
    N = N_values[i]
    println("Working on N: ", N)
    S_all = zeros(runs)
    for i in (1:runs)
        x = Int(N/2)
        println("Iteration number: ", i)
        H = rand_ham_sz0(N,Delta,J_perp,W)
        # e,v  = eigs(H, nev = 1, which=:SR)
        # v_full = convert_full(N,v)
        ex_pos = trunc(Int, 0.7*length(H[1,:]))
        e,v  = eigs(H, nev = ex_pos, which=:SR)
        v_ex = v[:,ex_pos]
        v_full = convert_full(N,v_ex)
        S_all[i] = entropy(v_full,x)
    end
    avg_S[i] = mean(S_all)
    std_S[i] = std(S_all)
end
save_object("Data/Q4/Ses9", avg_S)
save_object("Data/Q4/STDes9", std_S)
  



# N_values = range(8,16,step = 2)
# avg_S = load_object("Data/Simon/Long/Sgs0.5")
# std_S = load_object("Data/Simon/Long/STDgs0.5")
# avg_Se = load_object("Data/Q4/Ses0.5")
# std_Se = load_object("Data/Q4/STDes0.5")
# avg_S9 = load_object("Data/Simon/Long/Sgs9")
# std_S9 = load_object("Data/Simon/Long/STDgs9")
# avg_Se9 = load_object("Data/Simon/Long/Ses9")
# std_Se9 = load_object("Data/Simon/Long/STDes9")

# plt.scatter(N_values,avg_S,label= "Ground State W=0.5")
# plt.errorbar(N_values,avg_S, yerr = std_S)
# plt.scatter(N_values,avg_Se,label= "Excited State W=0.5")
# plt.errorbar(N_values,avg_Se, yerr = std_Se)
# plt.scatter(N_values,avg_S9,label= "Ground State W=9")
# plt.errorbar(N_values,avg_S9, yerr = std_S)
# plt.scatter(N_values,avg_Se9,label= "Excited State W=9")
# plt.errorbar(N_values,avg_Se9, yerr = std_Se)
# plt.title("Average Entropy Random Heisenberg")
# plt.ylabel(L"Average $S_A$")
# plt.xlabel(L"$N$")
# plt.legend()
# plt.grid()
# plt.show()



