import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os
from matplotlib2tikz import save as tikz_save


J = 1
k_B = 1

# Reading in data
if len(sys.argv) < 1:
	print("Please provide filename in commandline")
	sys.exit()
else:
	if sys.argv[1].endswith(".dat"):
		filename = "../Data/" + sys.argv[1]
	else:
		filename = "../Data/"  + sys.argv[1] + ".dat"

with open(filename, "r") as file:
    first_line = file.readline()
    T_min, T_max, dT = first_line.split()

T_list = np.arange(float(T_min), float(T_max), float(dT))

N, T, E, M, C_V, chi, absM, accepted_states = np.loadtxt(filename, unpack = True, skiprows=2)
N_samples = int(N[-1])

# Finding average expectation values for each temperature 
E_list = []
C_V_list = []
M_list = []
chi_list = []
for i in range(len(T_list)):
	E_list.append(np.sum(E[i*N_samples:(i+1)*N_samples])/N_samples)
	C_V_list.append(np.sum(C_V[i*N_samples:(i+1)*N_samples])/N_samples)
	M_list.append(np.sum(absM[i*N_samples:(i+1)*N_samples])/N_samples)
	chi_list.append(np.sum(chi[i*N_samples:(i+1)*N_samples])/N_samples)


def Z_analytical(T):
    beta = 1/(k_B*T)
    return (4*np.cosh(8*J*beta) + 12)

def E_analytical(T):
    beta = 1/(k_B*T)
    return -32/Z_analytical(T) *np.sinh(8*J*beta)

def C_V_analytical(T):
	beta = 1/(k_B*T)
	Z = Z_analytical(T)
	E = E_analytical(T)
	return beta/T*((256*J**2)/Z * np.cosh(8*beta*J) - E**2)

def M_analytical(T):
	beta = 1/(k_B*T)
	return (16+8*np.exp(8*beta*J))/Z_analytical(T)

def chi_analytical(T):
	beta = 1/(k_B*T)
	M = M_analytical(T)
	Z = Z_analytical(T)
	return beta*(32/Z * (1+np.exp(8*J*beta)) - M**2)



plt.figure()
plt.plot(T_list, E_list, 'ro', markersize=5, alpha=0.8, label='Monte Carlo')
plt.plot(T_list, E_analytical(T_list)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'$\langle E \rangle $')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_E.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.figure()
plt.plot(T_list, C_V_list, 'ro', markersize=5, alpha=0.8, label='Monte Carlo')
plt.plot(T_list, C_V_analytical(T_list)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'Specific heat $C_V$')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_Cv.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.figure()
plt.plot(T_list, M_list, 'ro', markersize=5, alpha=0.8, label='Monte Carlo')
plt.plot(T_list, M_analytical(T_list)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'$\langle |M| \rangle $')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_M.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.figure()
plt.plot(T_list, chi_list, 'ro', markersize=5, alpha=0.8, label='Monte Carlo')
plt.plot(T_list, chi_analytical(T_list)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'Susceptibility $\chi$')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_chi.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.show()
