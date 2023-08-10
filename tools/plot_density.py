import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

filename = "DENSITY"
N_x, N_y, N_z = np.loadtxt(filename, skiprows=0, max_rows=1, dtype=int)
Lx, Ly, Lz = 10.0, 8.0, 15.0
Y = np.linspace(0.0, Ly, N_y)
DATA = np.loadtxt(filename, skiprows=1)
DATA = np.reshape(DATA, np.size(DATA))
DATA = np.reshape(DATA, [N_x, N_y, N_z], order='F')
print("Total Charge: ", np.sum(DATA))
plt.plot(Y, DATA[N_x // 2, :, N_z // 2])
plt.savefig("density_slice")
