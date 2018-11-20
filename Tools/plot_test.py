import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os
from matplotlib2tikz import save as tikz_save


J = 1
k_B = 1

if len(sys.argv) < 1:
	print("Please provide filename in commandline")
	sys.exit()
else:
	if sys.argv[1].endswith(".dat"):
		filename = "../Data/" + sys.argv[1]
	else:
		filename = "../Data/"  + sys.argv[1] + ".dat"


T, E, C_V, M,  chi, absM= np.loadtxt(filename, unpack = True, skiprows=2)

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
plt.plot(T, E, 'ro', label='Monte Carlo')
plt.plot(T, E_analytical(T)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'$\langle E \rangle $')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_E.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.figure()
plt.plot(T, C_V, 'ro', label='Monte Carlo')
plt.plot(T, C_V_analytical(T)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'Specific heat $C_V$')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_Cv.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.figure()
plt.plot(T, absM, 'ro', label='Monte Carlo')
plt.plot(T, M_analytical(T)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'$\langle |M| \rangle $')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_M.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.figure()
plt.plot(T, chi, 'ro', label='Monte Carlo')
plt.plot(T, chi_analytical(T)/4, label='Analytical')
plt.xlabel('Temperature')
plt.ylabel(r'Susceptibility $\chi$')
plt.legend()
plt.grid()
#tikz_save("../Figures/2Dlattice_chi.tex", figureheight="\\figureheight", figurewidth="\\figureheight")

plt.show()

# plt.figure()
# plt.subplot(221)
# plt.plot(T_list, E_list, 'ro', label='Monte Carlo')
# plt.plot(T_list, E_analytical(T_list)/4, label='Analytical')
# plt.xlabel('Temperature')
# plt.ylabel('Energy')
# plt.grid()
#
# plt.subplot(222)
# plt.plot(T_list, C_V_list, 'ro', label='Monte Carlo')
# plt.plot(T_list, C_V_analytical(T_list)/4, label='Analytical')
# plt.xlabel('Temperature')
# plt.ylabel('Specific heat')
# plt.grid()
#
# plt.subplot(223)
# plt.plot(T_list, M_list, 'ro', label='Monte Carlo')
# plt.plot(T_list, M_analytical(T_list)/4, label='Analytical')
# plt.xlabel('Temperature')
# plt.ylabel('Magnetization')
# plt.grid()
#
# plt.subplot(224)
# plt.plot(T_list, chi_list, 'ro', label='Monte Carlo')
# plt.plot(T_list, chi_analytical(T_list)/4, label='Analytical')
# plt.xlabel('Temperature')
# plt.ylabel('Susceptibility')
# plt.grid()
#
# tikz_save("../Figures/2Dlattice_100000.tex", figureheight="\\figureheight", figurewidth="\\figureheight")
# plt.tight_layout()
# plt.subplots_adjust(top=0.9)
# plt.suptitle('N=%d' %(N_samples))
# plt.show()
