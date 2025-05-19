import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 🔹 Enter Your Data from LAMMPS Output
lattice_constants = np.array([
    3.40, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49,
    3.50, 3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59,
    3.60, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69,
    3.70, 3.71
])
potential_energies = np.array([
    -1670.798200, -1680.653700, -1689.911700, -1698.586800, -1706.693200, -1714.244700,
    -1721.255000, -1727.737500, -1733.705200, -1739.171100, -1744.144700, -1748.645900,
    -1752.682000, -1756.263800, -1759.404100, -1762.114500, -1764.406000, -1766.289600,
    -1767.776000, -1768.875700, -1769.599000, -1769.955800, -1769.956200, -1769.609800,
    -1768.926100, -1767.914200, -1766.583400, -1764.942400, -1763.000100, -1760.764900,
    -1758.245300, -1755.449300
])

# Convert potential energy to relative scale for visualization
potential_energies = potential_energies - min(potential_energies)

# 🔹 Define a Quadratic Function for Curve Fitting
def quadratic(l, p0, p1, p2):
    return p0 + p1 * l + p2 * l**2

# 🔹 Fit the Data with a Quadratic Curve
popt, _ = curve_fit(quadratic, lattice_constants, potential_energies)

# 🔹 Generate Smooth Line for Fitted Curve
l_fit = np.linspace(min(lattice_constants), max(lattice_constants), 300)

w_fit = quadratic(l_fit, *popt)

p0, p1, p2 = popt
print(f"p0: {p0}, p1: {p1}, p2: {p2}")

# 🔹 Plot the Data
plt.figure(figsize=(7, 5))
plt.scatter(lattice_constants, potential_energies, color='dodgerblue', label="Actual Data")
plt.plot(l_fit, w_fit, 'r--', label="Fitted Data")

# 🔹 Labels & Legends
plt.xlabel("Variation of Lattice Constant (Å)")
plt.ylabel("Potential Energy (eV)")
plt.legend()
plt.grid(True)

# 🔹 Save or Show the Plot
plt.savefig("fcc_curve_fit.png", dpi=300)
plt.show()
