import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 🔹 Full BCC Data from LAMMPS Output
lattice_constants = np.array([
    2.80, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.90, 2.91,
    2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.00, 3.01, 3.02, 3.03,
    3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.10, 3.11, 3.12, 3.13, 3.14, 3.15,
    3.16, 3.17, 3.18, 3.19, 3.20, 3.21
])
potential_energies = np.array([
    -871.121430, -873.030900, -874.603960, -875.851500, -876.784120, -877.412840, -877.766310,
    -877.910780, -877.777380, -877.375590, -876.714640, -875.803520, -874.650990, -873.265600,
    -871.655660, -869.829280, -867.794350, -865.558550, -863.125240, -860.665800, -858.166720,
    -855.487890, -852.636330, -849.618850, -846.442090, -843.112490, -839.636340, -836.019740,
    -832.268610, -828.388740, -824.385730, -820.265030, -816.031940, -811.691620, -807.249080,
    -802.709170, -798.076620, -793.356040, -788.551450, -783.669460, -778.710190, -773.635020
])

# 🔹 Define a Quadratic Function
def quadratic(l, p0, p1, p2):
    return p0 + p1 * l + p2 * l**2

# 🔹 Fit the Data
popt, _ = curve_fit(quadratic, lattice_constants, potential_energies)

# 🔹 Generate Smooth Fit Curve
l_fit = np.linspace(min(lattice_constants), max(lattice_constants), 300)
w_fit = quadratic(l_fit, *popt)

# 🔹 Plot the Data and Fit
plt.figure(figsize=(7, 5))
plt.scatter(lattice_constants, potential_energies, color='dodgerblue', label="Actual Data")
plt.plot(l_fit, w_fit, 'r--', label="Quadratic Fit")

# 🔹 Labels & Legends
plt.xlabel("Variation of Lattice Constant (Å)")
plt.ylabel("Potential Energy (eV)")
plt.legend()
plt.grid(True)
plt.title("BCC Lattice Energy vs. Lattice Constant")

# 🔹 Save or Show Plot
plt.savefig("bcc_quadratic_fit.png", dpi=300)
plt.show()
