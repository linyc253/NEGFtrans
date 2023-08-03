import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

filename = "grid3_slice.dat"

DATA = np.loadtxt(filename, skiprows=1)

X = DATA[:, 0]
V_real = DATA[:, 1]
V_new_real = DATA[:, 2]

plt.plot(X, V_real, marker='o', label="Non-Local")
plt.plot(X, V_new_real, "--", marker='x', label=r"Local / $N_x N_y$")
plt.legend()
plt.title("Comparison of Local and Non-Local code")
#plt.title("New version ENCUT = 150.0 eV")
plt.xlabel(r"$G_y$ ($G_x = 0$)")
plt.ylabel("Potential")
plt.savefig("grid3_slice1")
