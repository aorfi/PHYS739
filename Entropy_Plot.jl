using PyPlot
using JLD2
using LsqFit

# TFIM

# Plot all g values
# N = 15
# x_all = range(1,N-1)
# S_all = load_object("Data/Q3/TFIM_g0")
# S_all1 = load_object("Data/Q3/TFIM_g1")
# S_all2 = load_object("Data/Q3/TFIM_g2")
# plt.figure(figsize=(8,8))
# plt.plot(x_all, S_all[1:length(x_all)], label = L"g=0")
# plt.scatter(x_all, S_all[1:length(x_all)])
# plt.plot(x_all, S_all1[1:length(x_all)], label = L"g=1")
# plt.scatter(x_all, S_all1[1:length(x_all)])
# plt.plot(x_all, S_all2[1:length(x_all)], label = L"g=2")
# plt.scatter(x_all, S_all2[1:length(x_all)])
# plt.title(L"Entropy TFIM $N=15$", fontsize=18)
# plt.ylabel(L"$S_A$", fontsize=14)
# plt.xlabel(L"$x$", fontsize=14)
# plt.legend()
# plt.grid()
# plt.savefig("Figures/Q3/Entropy_TFIM.png")
# plt.show()

#Fitting 
# N = 15
# x_all = range(1,N-1)
# S_all = load_object("Data/Q3/TFIM_g1")
# # Fit Entropy Model
# @. model(x,p) = p[1]/3*log(N/pi*sin(pi*x/N)) + p[2]
# p0 = [0.5,0.0]
# fit = curve_fit(model,x_all,S_all[1:length(x_all)],p0)
# param = fit.param
# print(param)
# plt.figure(figsize=(8,8))
# plt.plot(x_all, S_all[1:length(x_all)])
# plt.scatter(x_all, S_all[1:length(x_all)])
# plt.plot(x_all, model(x_all,param),linestyle = "dashed", label = "Fit")
# plt.title(L"Entropy TFIM $N=15, g=1$", fontsize=18)
# plt.ylabel(L"$S_A$", fontsize=14)
# plt.xlabel(L"$x$", fontsize=14)
# plt.legend()
# plt.grid()
# # plt.savefig("Figures/Q3/Entropy_TFIM_fit.png")
# plt.show()

# # XXY 
# N = 14
# x_all = range(1,N-1)
# S_all = load_object("Data/Q3/XXY0")
# # Fit Entropy Model
# @. model(x,p) = p[1]/3*log(N/pi*sin(pi*(x/N)))+p[2]
# p0 = [0.5,0.0]
# fit = curve_fit(model,x_all,S_all[1:length(x_all)],p0)
# param = fit.param
# print(param)
# plt.figure(figsize=(8,8))
# plt.plot(x_all, S_all[1:length(x_all)])
# plt.scatter(x_all, S_all[1:length(x_all)])
# plt.plot(x_all, model(x_all,param),linestyle = "dashed", label = "Fit")
# plt.title(L"Entropy XY Model $N=16$", fontsize=18)
# plt.ylabel(L"$S_A$", fontsize=14)
# plt.xlabel(L"$x$", fontsize=14)
# plt.legend()
# plt.grid()
# plt.savefig("Figures/Q3/Entropy_XY_fit.png")
# plt.show()

# XXY 
N = 16
x_all = range(1,N-1)
S_all = load_object("Data/Q3/XXY1")
# Fit Entropy Model
@. model(x,p) = p[1]/3*log(N/pi*sin(pi*(x/N)))+p[2]
p0 = [0.5,0.0]
fit = curve_fit(model,x_all,S_all[1:length(x_all)],p0)
param = fit.param
print(param)
plt.figure(figsize=(8,8))
plt.plot(x_all, S_all[1:length(x_all)])
plt.scatter(x_all, S_all[1:length(x_all)])
plt.plot(x_all, model(x_all,param),linestyle = "dashed", label = "Fit")
plt.title(L"Entropy Heisenberg Model $N=16$", fontsize=18)
plt.ylabel(L"$S_A$", fontsize=14)
plt.xlabel(L"$x$", fontsize=14)
plt.legend()
plt.grid()
plt.savefig("Figures/Q3/Entropy_Hei_fit.png")
plt.show()