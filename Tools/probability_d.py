import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os
from matplotlib2tikz import save as tikz_save


if len(sys.argv) < 1:
	print("Please provide filename in commandline")
	sys.exit()
else:
	if sys.argv[1].endswith(".dat"):
		filename = "../Data/" + sys.argv[1]
	else:
		filename = "../Data/"  + sys.argv[1] + ".dat"

N, T, E, M, C_V, chi, absM, accepted_states = np.loadtxt(filename, unpack = True, skiprows=2)
N_samples = int(N[-1])
T = T[0]
k_B = 1
sigma_E = C_V*k_B*T**2


plt.figure()
plt.hist(E[int(N_samples*0.2):], bins=40)
plt.grid()
plt.xlabel("Energy")
plt.ylabel("Probability")
plt.title('Probability of energy at T=%.1f' %T)
#tikz_save("../Figures/" + filename[:-4] + "_P.tex", figureheight="\\figureheight", figurewidth="\\figurewidth")

plt.figure()
plt.plot(N[int(N_samples*0.2):], sigma_E[int(N_samples*0.2):])
plt.grid()
plt.ylabel(r'$\sigma_E$')
plt.xlabel('Monte Carlo cycles')
plt.show()
