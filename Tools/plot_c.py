import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os
from matplotlib2tikz import save as tikz_save


# Reading in data, note that the disordered filename has to be longer than the ordered filename
if len(sys.argv) < 2:
	print("Please provide filename for ordered and disordered initial lattice data in commandline")
	sys.exit()
else:
	if sys.argv[1].endswith(".dat"):
		filename1 = "../Data/" + sys.argv[1]
	if sys.argv[2].endswith(".dat"):
		filename2 = "../Data/" + sys.argv[2]
	else:
		filename1 = "../Data/"  + sys.argv[1] + ".dat"
		filename2 = "../Data/"  + sys.argv[2] + ".dat"


if len(filename2) > len(filename1):
	# Ordered:
	N_ordered, T_ordered, E_ordered, M_ordered, C_V_ordered, \
	chi_ordered, absM_ordered, accepted_states_ordered = np.loadtxt(filename1, unpack = True)
	T = T_ordered[0]

	#Disordered
	N_samples = int(N_ordered[-1])
	N_disordered, T_disordered, E_disordered, M_disordered, C_V_disordered, \
	chi_disordered, absM_disordered, accepted_states_disordered = np.loadtxt(filename2, unpack = True)

if len(filename1) > len(filename2):
	# Ordered:
	N_ordered, T_ordered, E_ordered, M_ordered, C_V_ordered, \
	chi_ordered, absM_ordered, accepted_states_ordered = np.loadtxt(filename2, unpack = True)
	T = T_ordered[0]

	#Disordered
	N_samples = int(N_ordered[-1])
	N_disordered, T_disordered, E_disordered, M_disordered, C_V_disordered, \
	chi_disordered, absM_disordered, accepted_states_disordered = np.loadtxt(filename1, unpack = True)


# Plotting every j'th point 
j = 100
plt.plot(N_ordered[::j], E_ordered[::j], 'b', label='Ordered initial lattice')
plt.plot(N_disordered[::j], E_disordered[::j], 'r', label='Disordered initial lattice')
plt.title('Temperature T=%.1f' %T)
plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'$\langle E \rangle $')
plt.legend()
plt.grid()
#tikz_save("../Figures/20D_lattice_E2.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")
plt.show()

plt.plot(N_ordered[::j], absM_ordered[::j], 'b', label='Ordered initial lattice')
plt.plot(N_disordered[::j], absM_disordered[::j], 'r', label='Disordered initial lattice')
plt.title('Temperature T=%.1f' %T)
plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'$\langle |M| \rangle $')
plt.legend()
plt.grid()
#tikz_save("../Figures/20D_lattice_M2.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")
plt.show()

plt.plot(N_ordered[::j], accepted_states_ordered[::j], 'b', label='Ordered initial lattice')
plt.plot(N_disordered[::j], accepted_states_disordered[::j], 'r', label='Disordered initial lattice')
plt.title('Temperature T=%.1f' %T)
plt.xlabel('Monte Carlo cycles')
plt.ylabel('Accepted states')
plt.legend()
plt.grid()
#tikz_save("../Figures/20D_lattice_accept2.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")
plt.show()
