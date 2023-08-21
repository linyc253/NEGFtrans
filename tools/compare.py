import numpy as np
filename_1 = "DENSITY"
filename_2 = "DENSITY_standard"
DATA_1 = np.loadtxt(filename_1, skiprows=1)
DATA_2 = np.loadtxt(filename_2, skiprows=1)
print("ERROR = ", np.sum(np.abs(DATA_1 - DATA_2)))