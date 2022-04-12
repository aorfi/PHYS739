using OrderedCollections
using SparseArrays
using Arpack
using PyPlot
using JLD2

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

function XXZ_sz0(N,Delta,J_perp)
    dict = sz0_dict(N)
    dim = length(dict)
    H = zeros(dim,dim)
    Diagonal = 0
    for ket in keys(dict)
        for SpinIndex in (0:N-1)
            # Boundary Condition
            if SpinIndex == N-1
                # SzSz term
                Sz_i = 2*((ket>>SpinIndex)&1)-1
                Sz_j = 2*((ket>>(0))&1)-1
                z_term = (Delta/4)*Sz_i*Sz_j
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
            bit = Int(2)^(SpinIndex)
            nn_bit = Int(2)^(SpinIndex+1)
            flipbit = ket ⊻ bit
            flipnn  = flipbit ⊻ nn_bit
            Sy_flip1 = (2*((flipbit>>SpinIndex)&1)-1)im
            Sy_flip2 = (2*((flipnn>>(SpinIndex+1))&1)-1)im
            Sy = Sy_flip1*Sy_flip2
            value = -(J_perp/4)*(1 + Sy)
            if value != 0.0
                index_original = dict[ket]
                index_flipped = dict[flipnn]
                H[index_flipped,index_original] += value
            end
            Szi = 2*((ket>>SpinIndex)&1)-1
            Szi_next = 2*((ket>>(SpinIndex+1))&1)-1
            z_term = (Delta/4)*Szi*Szi_next
            Diagonal += z_term
        end
        index_ket = dict[ket]
        H[index_ket,index_ket] = Diagonal
    end
    return H |> sparse
end

Delta = 0
J_perp = 1
N_values = range(2,16,step = 2)

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



energy = load_object("Data/Q2/energy1")
# Energy vs N plot
plt.plot(N_values, energy[1:length(N_values)])
plt.title(L"Energy of Ground-State $\Delta = 0$, $J = 1 $")
plt.ylabel(L"$\epsilon = E/N$")
plt.xlabel(L"$N$")
plt.show()
plt.savefig("Figures/3bGSEnergy.png")