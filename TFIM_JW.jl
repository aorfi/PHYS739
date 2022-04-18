using LinearAlgebra
using SparseArrays
using PyPlot

g = 1
J = 1
len = 200
k = range(-2*pi; stop=2*pi, length=len)
E = zeros(len)
E05 = zeros(len)
E15 = zeros(len)
for i in (1:len)
    E[i] = 2*J*sqrt((g-cos(k[i]))^2+(sin(k[i]))^2)
    E05[i] = 2*J*sqrt((0.5-cos(k[i]))^2+(sin(k[i]))^2)
    E15[i] = 2*J*sqrt((1.5-cos(k[i]))^2+(sin(k[i]))^2)
end
plt.plot(k, E, label = L"g=1")
plt.plot(k, E05, label = L"g=0.5")
plt.plot(k, E15, label = L"g=1.5")
plt.title("Jorgan-Wigner Energy")
plt.ylabel(L"$\epsilon_k$")
plt.xlabel(L"$k$")
plt.grid()
plt.legend()
plt.savefig("Figures/Q1/1JW.png")
plt.show()