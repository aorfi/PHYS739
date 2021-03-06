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
N = 15
x_all = range(1,N-1)
S_all = load_object("Data/Q3/TFIM_g1")
# Fit Entropy Model
# @. model(x,p) = p[1]/3*log(N/pi*sin(pi*x/N)) + p[2]
# p0 = [0.5,0.0]
# fit = curve_fit(model,x_all,S_all[1:length(x_all)],p0)
# param = fit.param
# sigma = stderror(fit)
# println(param)
# println(sigma)
# x = range(1,N-1, length= 1000)
# plt.figure(figsize=(8,8))
# plt.scatter(x_all, S_all[1:length(x_all)], label = "ED Data")
# plt.plot(x, model(x,param),linestyle = "dashed", label = "Functional Form Fit", color = "green")
# plt.title(L"Entropy TFIM $N=15, g=1$", fontsize=18)
# plt.ylabel(L"$S_A$", fontsize=14)
# plt.xlabel(L"$x$", fontsize=14)
# plt.legend()
# plt.grid()
# plt.savefig("Figures/Q3/Entropy_TFIM_fit.png")
# plt.show()




# # XXY 
# Fitting
N = 16
x_all = range(1,N-1)
S_all = load_object("Data/Q3/XXY0")
# Fit Entropy Model
@. model(x,p) = p[1]/3*log(N/pi*sin(pi*(x/N)))+p[2]
p0 = [0.5,0.0]
fit = curve_fit(model,x_all,S_all[1:length(x_all)],p0)
param = fit.param
sigma = stderror(fit)
println(param)
println(sigma)
# plt.figure(figsize=(8,8))
x = range(1,N-1, length= 1000)
plt.scatter(x_all, S_all[1:length(x_all)], label = "ED Data")
plt.plot(x, model(x,param),linestyle = "dashed", label = "Functional Form Fit", color = "green")
plt.title(L"Entropy XY Model  $N=16$")
plt.ylabel(L"$S_A$")
plt.xlabel(L"$x$")
plt.legend()
plt.grid()
plt.savefig("Figures/Q3/Entropy_XY_fit.png")
plt.show()