# include("TFIM.jl")
using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

g_beg = range(0.001; stop=0.75, length=10)
g_fine = range(0.76,stop=1.25, length = 50)
g_end = range(1.325; stop=2, length=10)
g = vcat(g_beg,g_fine,g_end )
len= length(g)
# N = 16



# # Plot energy/magnetization  of low lying states 
# energy_all = load_object("Data/Final/energy_all")
# energy = energy_all[N-1,:,:]
# # Energy vs g plot
# plt.plot(g, energy[1:len,1], label = "Ground State")
# plt.plot(g, energy[1:len,2],label = "First Excited State")
# plt.plot(g, energy[1:len,3], label = "Second Excited State")
# plt.title(string("Energy of Low-Lying States N = ", N))
# plt.ylabel(L"$\epsilon = E/N$")
# plt.xlabel(L"$h/J$")
# plt.legend()
# plt.grid()
# plt.savefig("Figures/Q1/1aEnergy.png")
# plt.show()


# # # Mag vs g plot
# mag_all = load_object("Data/Final/mag_all")
# mag = mag_all[N-1,:,:]
# plt.plot(g, mag[1:len,1], label = "Ground State")
# plt.plot(g, mag[1:len,2], label = "First Excited State")
# plt.plot(g, mag[1:len,3], label = "Second Excited State")
# plt.title(string("Magnetization of Low-Lying States N = ", N))
# plt.ylabel(L"$|m|/N$")
# plt.xlabel(L"$h/J$")
# plt.legend()
# plt.grid()
# plt.savefig("Figures/Q1/1aMag.png")
# plt.show()

# # Ground-state energy for different N
# N_values = (2:16)
# N_max = last(N_values)
# mag_all = load_object("Data/Final/mag_all")
# for i in (1:length(N_values))
#     N = N_values[i]
#     plt.plot(g, mag_all[i,1:len,1], label = string("N= ", N))
# end
# plt.title("Magnetization of Ground State ")
# plt.ylabel(L"$|m|/N$")
# plt.xlabel(L"$h/J$")
# plt.legend()
# plt.grid()
# plt.savefig("Figures/1cMagScaling.png")
# plt.show()

# Get gap
N_values = (5:16)
N_max = last(N_values)
energy_all = load_object("Data/Final/energy_all")

indexCP = findall(x->x==1, g)[1]
g_values = [g[indexCP-2],g[indexCP-1],g[indexCP],g[indexCP+1],g[indexCP+2]]
gap_all = zeros(length(N_values),length(g_values))
for j in (1:length(N_values))
    N = N_values[j]
    for i in (1:length(g_values))
        gap_all[j,i] = abs(N*(energy_all[N-1,Int(indexCP-3+i),1]-energy_all[N-1,indexCP-3+i,2]))
    end
end

@. model(x,p) = p[1]*x+p[2]
p0 = [2.0,0]
xdata = 1 ./ (N_values)
params_all = zeros(length(g_values))
for i in (1:length(g_values))
    println("Working on ", i)
    fit = curve_fit(model,xdata,gap_all[1:length(N_values),i,1],p0)
    params_all[i]= fit.param[2]
end
fit = curve_fit(model,xdata,gap_all[1:length(N_values),3,1],p0)
param = fit.param
sigma = stderror(fit)

# for i in (1:length(g_values))
#     g_value = g_values[i]
#     plt.plot(1 ./ (N_values), gap_all[1:length(N_values),i,1], label = string("g = ", g_value))
# end

# plt.plot(1 ./ (N_values), gap_all[1:length(N_values),3,1])
# plt.scatter(1 ./ (N_values), gap_all[1:length(N_values),3,1])
# plt.plot(1 ./ (N_values),model(xdata,param), linestyle = "dashed", label = "Linear Fit")
# plt.title(L"Energy Gap of Ground State at $g=1$")
# plt.ylabel(L"$|E_1-E_0|$")
# plt.xlabel(L"$1/N$")
# plt.legend()
# plt.grid()
# plt.savefig("Figures/Q1/1dGapScaling.png")
# plt.show()


fit = curve_fit(model,g_values,params_all,p0)
param = fit.param
sigma = stderror(fit)
print(param)
print(sigma)
x = range(0.98,1.02,length=1000)
plt.figure(figsize=(7,7))
# plt.plot(g_values, params_all)
plt.scatter(g_values, params_all)
plt.plot(x, model(x,param), color = "green",linestyle = "dashed")

plt.title("Energy Gap of Ground State")
plt.ylabel("Fit y-intercept")
plt.xlabel(L"$h/J$")
# plt.legend()
plt.grid()
plt.savefig("Figures/Q1/1dyintercept.png")
plt.show()