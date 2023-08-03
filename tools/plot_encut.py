import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

filename = "encut_test.dat"

DATA = np.loadtxt(filename, skiprows=1)

X = DATA[:, 0]
Y = DATA[:, 1]

plt.plot(X, Y)
plt.title(r"ENCUT vs Error (Cell: 10.0 bohr X 8.0 bohr)")
plt.xlabel("ENCUT (eV)")
plt.ylabel("Error (hartree)")
plt.ylim(0, 39)
plt.savefig("encut")
