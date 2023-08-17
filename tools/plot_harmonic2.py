import numpy as np
import matplotlib.pyplot as plt
import re
plt.rc("font", size=14)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("figure", titlesize=16)

Lx, Ly, Lz = 10.0, 8.0, 15.0


filename = "POTENTIAL"
N_x, N_y, N_z = np.loadtxt(filename, skiprows=0, max_rows=1, dtype=int)
POTENTIAL = np.loadtxt(filename, skiprows=1)
POTENTIAL = np.reshape(POTENTIAL, np.size(POTENTIAL))
POTENTIAL = np.reshape(POTENTIAL, [N_x, N_y, N_z], order='F')

filename = "../test_4/DENSITY"
DENSITY_1 = np.loadtxt(filename, skiprows=1)
DENSITY_1 = np.reshape(DENSITY_1, np.size(DENSITY_1))
DENSITY_1 = np.reshape(DENSITY_1, [N_x, N_y, N_z], order='F')

filename = "DENSITY"
DENSITY_2 = np.loadtxt(filename, skiprows=1)
DENSITY_2 = np.reshape(DENSITY_2, np.size(DENSITY_2))
DENSITY_2 = np.reshape(DENSITY_2, [N_x, N_y, N_z], order='F')

DENSITY = DENSITY_2 - DENSITY_1

# Theoretical calculation
Y_dense = np.linspace(0, Ly, 1000)
omega = 0.47316
y_0 = omega ** (-1/2)
THEORY = 2 * ((Y_dense - Ly / 2) / 0.529177210903 / y_0) * \
    np.exp(- ((Y_dense - Ly / 2) / 0.529177210903 / y_0) ** 2 / 2) / (np.pi ** (3/4)) / (y_0 ** (3/2)) / 2 ** (1/2)
THEORY = 2 * (THEORY ** 2)
THEORY /= 0.529177210903 ** 3

Y = np.arange(1, N_y + 1) * Ly / N_y
print("Total Charge: ", np.sum(DENSITY) * (Lx * Ly * Lz) / (N_x * N_y * N_z))

fig, ax1 = plt.subplots()

ax1.plot(Y, DENSITY[N_x // 2 - 1, :, N_z // 2 - 1], marker='x', c='tab:orange', label='Calculated Density')
ax1.plot(Y_dense, THEORY, '--', c='k', label='Theoretical Density', linewidth=1)
ax1.axhline(y=0.0, linestyle='--', color='tab:gray', linewidth=0.5)
ax2 = ax1.twinx()
ax2.plot(Y, POTENTIAL[N_x // 2 - 1, :, N_z // 2 - 1], c='tab:blue', label='Potential')

ax1.set_xlabel(r'Y ($\AA$)')
ax1.set_ylabel(r'Density ($1/\AA^3$)', c='tab:orange')
ax1.tick_params(axis='y', labelcolor='tab:orange')
ax1.set_ylim(-0.05, 1.4)
ax2.set_ylabel('Potential (eV)', c='tab:blue')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.legend()
plt.savefig("multi_slice")
plt.show()
