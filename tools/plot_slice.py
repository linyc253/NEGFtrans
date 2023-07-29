import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

filename = "grid2_slice.dat"

DATA = np.loadtxt(filename, skiprows=1)

X = DATA[:, 0]
V_real = DATA[:, 1]
V_new_real = DATA[:, 2]

plt.plot(X, V_real, label="V_real")
plt.plot(X, V_new_real, "--", label="V_new")
plt.legend()
plt.title("New version ENCUT = 150.0 eV")
plt.xlabel("X")
plt.ylabel("Potential")
plt.savefig("grid2_slice2")
