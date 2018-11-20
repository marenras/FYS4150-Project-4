import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os
from matplotlib2tikz import save as tikz_save


# Reading in data
folder = "../Data/"
if len(sys.argv) < 1:
	print("Please provide filename in commandline")
	sys.exit()
else:
	if sys.argv[1].endswith(".dat"):
		filename = sys.argv[1]
	else:
		filename = sys.argv[1] + ".dat"


N, T, E, M, C_V, chi, absM, accepted_states = np.loadtxt(folder + filename, unpack = True, skiprows=2)
N_samples = int(N[-1])
T = T[0]
k_B = 1
sigma_E = C_V*k_B*T**2


# Plotting every j'th point 
j = 100
plt.figure()
plt.hist(E[int(N_samples*0.2):], bins=20)
plt.grid()
plt.xlabel("Energy")
plt.ylabel("Probability")
plt.title('Temperature T=%.1f' %T)
#tikz_save("../Figures/Project_4d/" + filename[:-4] + "_P.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.plot(N[int(N_samples*0.2):][::j], sigma_E[int(N_samples*0.2):][::j])
plt.grid()
plt.ylabel(r'$\sigma_E^2$')
plt.xlabel('Monte Carlo cycles')
plt.title('Temperature T=%.1f' %T)
#tikz_save("../Figures/Project_4d/" + filename[:-4] + "_sigma.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")
plt.show()
