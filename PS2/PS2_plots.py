from cProfile import label
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('seaborn')
matplotlib.rcParams.update({'font.size': 25})



U = 4

mu_values = np.arange(-5,5, 0.05)

T = 2
beta = 1/T
pT1 = (1/(1+2*np.exp(beta*mu_values)+np.exp(2*beta*mu_values-beta*U)))*(2*np.exp(beta*mu_values)+2*np.exp(2*beta*mu_values-beta*U))
T = 0.5
beta = 1/T
pT2 = (1/(1+2*np.exp(beta*mu_values)+np.exp(2*beta*mu_values-beta*U)))*(2*np.exp(beta*mu_values)+2*np.exp(2*beta*mu_values-beta*U))

plt.plot(mu_values, pT1, label = "T = 2")
plt.plot(mu_values, pT2, label = "T = 0.5")
plt.ylabel("Occupation "+ r"$\rho$")
plt.xlabel("Chemical Potential "+r"$\mu$")
plt.title("Single Site Limit Hubbard Model")
plt.legend()
plt.grid()
plt.show()

