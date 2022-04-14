using LinearAlgebra
using SparseArrays
using PyPlot

g = 1
J = 1
len = 200
k = range(-2*pi; stop=2*pi, length=len)
E = zeros(len)
for i in (1:len)
    E[i] = 2*J*sqrt((g-cos(k[i]))^2+(sin(k[i]))^2)
end
plt.plot(k, E)
plt.title("EJW Energy")
plt.ylabel(L"$\epsilon$")
plt.xlabel(L"$k$")
plt.legend()
# plt.savefig("Figures/1aMag.png")
plt.show()