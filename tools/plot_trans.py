import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

filename = "TRANSMISSION"

DATA = np.loadtxt(filename, skiprows=1)

E = DATA[:, 0]
TAU = DATA[:, 1]

plt.plot(E[:100], TAU[:100])
plt.title("Transmission Coefficient")
plt.xlabel("E (eV)")
plt.ylabel(r"$\tau$")
plt.savefig("plot")
