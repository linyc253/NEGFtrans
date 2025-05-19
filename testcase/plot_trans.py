import numpy as np
import matplotlib.pyplot as plt
import re
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams['font.family'] = "STIXGeneral"
plt.rcParams['font.size'] = 12

filename = "TRANSMISSION"

DATA = np.loadtxt(filename, skiprows=1)

E = DATA[:, 0]
TAU = DATA[:, 1]

for i in range(7):
    plt.axhline(i, color="gray", linewidth=0.5)
plt.plot(E, TAU)


plt.title("Transmission Coefficient")
plt.xlabel("E (eV)")
plt.ylabel(r"$\tau$")
plt.savefig("transmission")
