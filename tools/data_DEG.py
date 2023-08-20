
hartree = 27.211386245988 # in eV
k_B = 8.617333262E-5 # in eV
k_B = k_B / hartree # convert to atomic unit
a_0 = 0.529177210903 # Bohr radius (in Angstrom)
pi = 3.141592653589793

RS_LIST = [3.0, 4.0, 5.0, 6.0]

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
CHARGE_LIST = [215.36234083747746, 90.739865169592548, 46.251666620514953, 26.981247639128554]
volume = 10.0 * 6.0 * 60.0# (Angstrom^3)


print(" RS      RHO (1/Angstrom^3)    MU (eV)      Calculated RHO       Error(%)")
for i in range(len(RS_LIST)):
    rs = RS_LIST[i]
    charge = CHARGE_LIST[i]
    rho = 1 / (4 / 3 * pi * rs**3) / a_0 ** 3
    rho_calculated = charge / volume
    error = (rho_calculated - rho) / rho * 100
    mu = (3 * rho * pi**2) ** (2 / 3) / 2
    print("{:<0.2f}         {:<0.6f}         {:<0.3f}          {:<0.6f}             {:<0.2f} ".format(rs, rho, mu * hartree, rho_calculated, error))
