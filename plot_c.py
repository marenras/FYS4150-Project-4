import numpy as np
import sys as sys
import matplotlib.pyplot as plt
import os


if len(sys.argv) < 1:
	print("Please provide filename in commandline")
	sys.exit()
else:
	if sys.argv[1].endswith(".dat"):
		filename = "Data/" + sys.argv[1]
	else:
		filename = "Data/"  + sys.argv[1] + ".dat"

#filename = '../test20'

N, T, E, M, C_V, chi, absM, accepted_states = np.loadtxt(filename, unpack = True)
T = T[0]
N_samples = int(N[-1])
#
# N, T, E, E2, M, M2, Mabs = np.loadtxt(filename, unpack = True)
# T = T[0]

plt.plot(N, E)
plt.title('Temperature T=%.1f' %T)
plt.grid()
plt.show()

# plt.plot(N, M)
# plt.title('Temperature T=%.1f' %T)
# plt.grid()
# plt.show()
