import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

hartree = 27.211386245988 # in eV
k_B = 8.617333262E-5 # in eV
k_B = k_B / hartree # convert to atomic unit
a_0 = 0.529177210903 # Bohr radius (in Angstrom)
pi = 3.141592653589793


# Data to be filled
# &INPUT
# TEMPERATURE = 20
# LX = 10.0
# LY = 6.0
# LZ = 60.0
# MU = 5.56800678
# NKX = 10
# NKY = 17
# /
K_LIST = np.array([1, 8, 40, 170, 680])
CHARGE = np.array([172.38153240995572, 216.65332014805426, 215.26492587575447, 215.36234083747746, 214.98307978429531])
rs = 3.0
volume = 10.0 * 6.0 * 60.0# (Angstrom^3)

rho = 1 / (4 / 3 * pi * rs**3) / a_0 ** 3
RHO_CALCULATED = CHARGE / volume

plt.plot(K_LIST, RHO_CALCULATED / rho, marker='o', color='tab:blue')
plt.axhline(1.0, linestyle='--', color='tab:orange')
plt.xlabel("Number of K-Point")
plt.ylabel(r"$\rho / \rho_{theory}$")
plt.title("K-Point Convergence Test")
plt.show()
