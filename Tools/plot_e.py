import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os
from matplotlib2tikz import save as tikz_save


# Load in data
T, E_40, C_V_40, M_40, chi_40, absM_40 		= np.loadtxt('../Data/40_lattice.dat', unpack = True)
T, E_60, C_V_60, M_60, chi_60, absM_60		= np.loadtxt('../Data/60_lattice.dat', unpack = True)
T, E_80, C_V_80, M_80, chi_80, absM_80 		= np.loadtxt('../Data/80_lattice.dat', unpack = True)
T, E_100, C_V_100, M_100, chi_100, absM_100 = np.loadtxt('../Data/100_lattice.dat', unpack = True)


# Finding possible T_c values for each lattice size
Tc_possible = []

for C_V in [C_V_40, C_V_60, C_V_80, C_V_100]:
	index = C_V.tolist().index(np.max(C_V))
	Tc_possible.append(T[index])

print(Tc_possible)

Tc_40_60 = (60*Tc_possible[1]-40*Tc_possible[0])/(60-40)
Tc_60_80 = (80*Tc_possible[2]-60*Tc_possible[1])/(80-60)
Tc_80_100 = (100*Tc_possible[3]-80*Tc_possible[2])/(100-80)
Tc_40_80 = (80*Tc_possible[2]-40*Tc_possible[0])/(80-40)
Tc_40_100 = (100*Tc_possible[3]-40*Tc_possible[0])/(100-40)
Tc_60_100 = (100*Tc_possible[3]-60*Tc_possible[1])/(100-60)

Tc_list = [Tc_40_60, Tc_60_80, Tc_80_100, Tc_40_80, Tc_40_100, Tc_60_100]

# Finding the average off all the possible critical temperatures 
print(sum(Tc_list)/len(Tc_list))



plt.figure()
plt.plot(T, E_40, 'o', markersize=4, label=r"$40 \times 40$ lattice")
plt.plot(T, E_60, 'o', markersize=4, label=r"$60 \times 60$ lattice")
plt.plot(T, E_80, 'o', markersize=4, label=r"$80 \times 80$ lattice")
plt.plot(T, E_100,'o', markersize=4, label=r"$100 \times 100$ lattice")
plt.legend()
plt.xlabel('Temperature')
plt.ylabel(r'$\langle E \rangle $')
plt.grid()
#tikz_save("../Figures/4e/phase_E.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.plot(T, C_V_40, 'o', markersize=4, label=r"$40 \times 40$ lattice")
plt.plot(T, C_V_60, 'o', markersize=4, label=r"$60 \times 60$ lattice")
plt.plot(T, C_V_80, 'o', markersize=4, label=r"$80 \times 80$ lattice")
plt.plot(T, C_V_100,'o', markersize=4, label=r"$100 \times 100$ lattice")
plt.legend()
plt.xlabel('Temperature')
plt.ylabel(r'Specific heat $C_V$')
plt.grid()
#tikz_save("../Figures/4e/phase_C_V.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.plot(T, absM_40, 'o', markersize=4, label=r"$40 \times 40$ lattice")
plt.plot(T, absM_60, 'o', markersize=4, label=r"$60 \times 60$ lattice")
plt.plot(T, absM_80, 'o', markersize=4, label=r"$80 \times 80$ lattice")
plt.plot(T, absM_100,'o', markersize=4, label=r"$100 \times 100$ lattice")
plt.legend()
plt.xlabel('Temperature')
plt.ylabel(r'$\langle |M| \rangle $')
plt.grid()
#tikz_save("../Figures/4e/phase_absM.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.plot(T, chi_40, 'o', markersize=4, label=r"$40 \times 40$ lattice")
plt.plot(T, chi_60, 'o', markersize=4, label=r"$60 \times 60$ lattice")
plt.plot(T, chi_80, 'o', markersize=4, label=r"$80 \times 80$ lattice")
plt.plot(T, chi_100,'o', markersize=4, label=r"$100 \times 100$ lattice")
plt.legend()
plt.xlabel('Temperature')
plt.ylabel(r'Susceptibility $\chi$')
plt.grid()
#tikz_save("../Figures/4e/phase_chi.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")


plt.figure()
plt.plot(T, M_40, 'o', markersize=4, label=r"$40 \times 40$ lattice")
plt.plot(T, M_60, 'o', markersize=4, label=r"$60 \times 60$ lattice")
plt.plot(T, M_80, 'o', markersize=4, label=r"$80 \times 80$ lattice")
plt.plot(T, M_100,'o', markersize=4, label=r"$100 \times 100$ lattice")
plt.legend()
plt.xlabel('Temperature')
plt.ylabel(r'$\langle M \rangle $')
plt.grid()
#tikz_save("../Figures/4e/phase_M.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

#plt.show()
