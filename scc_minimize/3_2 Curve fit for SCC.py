import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 🔹 High-resolution SC Data from LAMMPS Output (2.36–2.415 Å)
lattice_constants = np.array([
    2.36000, 2.36250, 2.36500, 2.36750, 2.37000, 2.37250, 2.37500,
    2.37750, 2.38000, 2.38250, 2.38500, 2.38750, 2.39000, 2.39250,
    2.39500, 2.39750, 2.40000, 2.40250, 2.40500, 2.40750, 2.41000,
    2.41250, 2.41500
])

potential_energies = np.array([
    -387.613370, -387.726240, -387.827130, -387.916150, -387.993440,
    -388.059120, -388.113310, -388.156140, -388.187730, -388.208190,
    -388.217660, -388.216250, -388.204070, -388.181250, -388.147890,
    -388.104120, -388.050060, -387.985800, -387.911470, -387.827170,
    -387.733020, -387.629130, -387.515610
])

# 🔹 Normalize energy for visualization
potential_energies = potential_energies - min(potential_energies)

# 🔹 Define Quadratic Function for Curve Fitting
def quadratic(l, p0, p1, p2):
    return p0 + p1 * l + p2 * l**2

# 🔹 Fit the Data with a Quadratic Curve
popt, _ = curve_fit(quadratic, lattice_constants, potential_energies)

# 🔹 Generate Smooth Line for Fitted Curve
l_fit = np.linspace(min(lattice_constants), max(lattice_constants), 300)
w_fit = quadratic(l_fit, *popt)

# 🔹 Print coefficients and equilibrium lattice constant
p0, p1, p2 = popt
l_min = -p1 / (2 * p2)
print(f"p0: {p0:.6f}, p1: {p1:.6f}, p2: {p2:.6f}")
print(f"Estimated equilibrium lattice constant (a₀): {l_min:.5f} Å")

# 🔹 Plot the Data
plt.figure(figsize=(7, 5))
plt.scatter(lattice_constants, potential_energies, color='dodgerblue', label="Actual Data")
plt.plot(l_fit, w_fit, 'r--', label="Fitted Quadratic")

# 🔹 Labels & Legends
plt.xlabel("Lattice Constant (Å)")
plt.ylabel("Relative Potential Energy (eV)")
plt.title("Quadratic Fit: SC Cu Potential Energy vs. Lattice Constant")
plt.legend()
plt.grid(True)

# 🔹 Save or Show the Plot
plt.savefig("sc_quadratic_fit_highres.png", dpi=300)
plt.show()
