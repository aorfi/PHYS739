using OrderedCollections
using SparseArrays
using Arpack
# using PyPlot
using JLD2
using LsqFit

function sz0_dict(N)
    dict = OrderedDict()
    dim = 2^N
    sz = zeros(dim)
    index = 1
    for ket in (0:dim-1)
        for SpinIndex in (0:N-1)
            sz[ket+1] += 2*((ket>>SpinIndex)&1)-1
            end
        if sz[ket+1]==0
            dict[ket] = index
            index +=1
        end
    end
    return dict
end

function XXZ(N,Delta,J_perp)
    dim = 2^N
    H = zeros(dim,dim)
    for ket in (0:dim-1)
        Diagonal = 0
        for SpinIndex in (0:N-1)
            # Boundary Condition
            if SpinIndex == N-1
                # SzSz term
                Sz_i = 2*((ket>>SpinIndex)&1)-1
                Sz_j = 2*((ket>>(0))&1)-1
                z_term = -(Delta/4)*Sz_i*Sz_j
                Diagonal += z_term
                # SxSx term
                bit_i = 2^(SpinIndex)
                bit_j = 2^(0)
                Sx_i = ket ⊻ bit_i
                Sx_j  = Sx_i ⊻ bit_j
                # SySy term
                Sy_i = (2*((Sx_j>>SpinIndex)&1)-1)im
                Sy_j = (2*((Sx_j>>(0))&1)-1)im
                Sy_all = Sy_i*Sy_j
                H[Sx_j+1,ket+1]  += -(J_perp/4)*(1 + Sy_all)
                break
            end
            # Get binary
            bit = Int(2)^(SpinIndex)
            nn_bit = Int(2)^(SpinIndex+1)
            # SzSz term
            Sz_i = 2*((ket>>SpinIndex)&1)-1
            Sz_j = 2*((ket>>(SpinIndex+1))&1)-1
            z_term = -(Delta/4)*Sz_i*Sz_j
            Diagonal += z_term
            # SxSx term
            Sx_i = ket ⊻ bit
            Sx_j  = Sx_i ⊻ nn_bit
            # SySy term
            Sy_i = (2*((Sx_i>>SpinIndex)&1)-1)im
            Sy_j = (2*((Sx_j>>(SpinIndex+1))&1)-1)im
            Sy = Sy_i*Sy_j
            H[Sx_j+1,ket+1]  += -(J_perp/4)*(1 + Sy)
        end
        # Update
        H[ket+1,ket+1] = Diagonal
    end
    return H |> sparse
end


function XXZ_sz0(N,Delta,J_perp)
    dict = sz0_dict(N)
    dim = length(dict)
    H = zeros(dim,dim)
    for ket in keys(dict)
        Diagonal = 0
        for SpinIndex in (0:N-1)
            # Boundary Condition
            if SpinIndex == N-1
                # SzSz term
                Sz_i = 2*((ket>>SpinIndex)&1)-1
                Sz_j = 2*((ket>>(0))&1)-1
                z_term = -(Delta/4)*Sz_i*Sz_j
                Diagonal += z_term
                # SxSx term
                bit_i = 2^(SpinIndex)
                bit_j = 2^(0)
                Sx_i = ket ⊻ bit_i
                Sx_j  = Sx_i ⊻ bit_j
                # SySy term
                Sy_i = (2*((Sx_j>>SpinIndex)&1)-1)im
                Sy_j = (2*((Sx_j>>(0))&1)-1)im
                Sy_all = Sy_i*Sy_j
                value = -(J_perp/4)*(1 + Sy_all)
                # Find index to change
                if value != 0.0
                    index_original = dict[ket]
                    index_flipped = dict[Sx_j]
                    H[index_flipped,index_original] += value
                end
                break
            end
            # Get binary
            bit = Int(2)^(SpinIndex)
            nn_bit = Int(2)^(SpinIndex+1)
            # SxSx term
            Sx_i = ket ⊻ bit
            Sx_j  = Sx_i ⊻ nn_bit
            # SySy term
            Sy_i = (2*((Sx_i>>SpinIndex)&1)-1)im
            Sy_j = (2*((Sx_j>>(SpinIndex+1))&1)-1)im
            Sy = Sy_i*Sy_j
            value = -(J_perp/4)*(1 + Sy)
            if value != 0.0
                index_original = dict[ket]
                index_flipped = dict[Sx_j]
                H[index_flipped,index_original] += value
            end
            # SzSz term
            Sz_i = 2*((ket>>SpinIndex)&1)-1
            Sz_j = 2*((ket>>(SpinIndex+1))&1)-1
            z_term = -(Delta/4)*Sz_i*Sz_j
            Diagonal += z_term
        end
        # Update
        index_ket = dict[ket]
        H[index_ket,index_ket] = Diagonal
    end
    return H |> sparse
end

function convert_full(N,v)
    dict = sz0_dict(N)
    positions = dict.keys
    v_orginal = zeros(2^N,1)
    for i in (1:length(positions))
        v_orginal[positions[i]+1] = v[i]
    end
    return v_orginal
end

function Hamiltonian_Prob2(N,Delta,J_perp)
    dim = (2)^N
    H = zeros(dim,dim)
    magz = zeros(dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        for SpinIndex in (0:N-1)
            magz[ket+1] += 2*((ket>>SpinIndex)&1)-1
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
        H[ket+1,ket+1] += -1 * (Diagonal)
    end
    H_sparse = H |> sparse
    return H_sparse,magz
end

# N = 4
# Delta = -1
# J_perp = -1
# H,magz = Hamiltonian_Prob2(N,Delta,J_perp)
# e,v  = eigs(H, nev = 1, which=:SR)
# display(v)
# H = XXZ_sz0(N,Delta,J_perp)
# e,v  = eigs(H, nev = 1, which=:SR)
# dict = sz0_dict(N)
# positions = dict.keys
# v_orginal = zeros(2^N,1)
# for i in (1:length(positions))
#     v_orginal[positions[i]+1] = v[i]
# end
# display(v_orginal)
# display(v-v_orginal)
# # v_conv = convert_full(N,v)
# # display(v_conv )




# Delta = -1
# J_perp = -1
# N_values = range(3,10)

# #Get data
# energy = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     H,magz = Hamiltonian_Prob2(N,Delta,J_perp)
#     e,v  = eigs(H, nev = 1, which=:SR)
#     energy[j] = e[1][1]/N 
# end
# plt.plot(N_values,energy[1:length(N_values)])
# plt.show()




# Delta = -1
# J_perp = -1
# N_values = range(2,16,step = 2)

# #Get data
# energy = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     if N%2 != 0
#         continue
#     end
#     H = XXZ_sz0(N,Delta,J_perp)
#     e,v  = eigs(H, nev = 1, which=:SR)
#     energy[j] = e[1][1]/N 
# end
# save_object("Data/Q2/energy1", energy)


# energy = load_object("Data/Q2/energy1")
# # @. model(x,p) = p[1]*exp(p[2]*N_values+p[3])+p[4]
# # p0 = [-0.1,-0.4,2.0,-0.3]
# # fit = curve_fit(model,N_values,energy[2:length(N_values)+1],p0)
# # param = fit.param
# # print(param)

# # Energy vs N plot
# plt.figure(figsize=(7,7))
# plt.plot(N_values, energy[1:length(N_values)])
# plt.scatter(N_values, energy[1:length(N_values)])
# plt.plot(N_values, (1/4-log(2))*ones(length(N_values)), linestyle = "dashed", color= "black")
# # plt.plot(N_values, model(N_values,param), linestyle = "dashed")
# plt.title(L"Energy of Ground-State $\Delta = -1$, $J = -1 $")
# plt.ylabel(L"$\epsilon = E/N$")
# plt.xlabel(L"$N$")

# plt.grid()
# plt.savefig("Figures/Q2/energy1.png")
# plt.show()