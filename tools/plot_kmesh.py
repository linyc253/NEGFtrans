import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

a_0 = 0.529177210903 # Bohr radius (in Angstrom)

filename = "OUTPUT"

DATA = np.loadtxt(filename, skiprows=43, max_rows=170)

KX = DATA[:, 0]
KY = DATA[:, 1]


plt.scatter(KY, KX)
plt.title(r"K-point Mesh: $10 \times 17$")

plt.ylabel(r"$k_x (bohr^{-1})$")
plt.ylim(- np.pi / 10 * a_0, np.pi / 10 * a_0)

plt.xlabel(r"$k_y (bohr^{-1})$")
plt.xlim(- np.pi / 6 * a_0, np.pi / 6 * a_0)

plt.axhline(0, c='k', linewidth=1)
plt.axvline(0, c='k', linewidth=1)
plt.gca().set_aspect('equal')
plt.savefig("kmesh")
plt.show()