import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# === Step 1: Input Data ===
strain_values = np.array([
    -0.080, -0.079, -0.078, -0.077, -0.076, -0.075, -0.074, -0.073, -0.072, -0.071,
    -0.070, -0.069, -0.068, -0.067, -0.066, -0.065, -0.064, -0.063, -0.062, -0.061,
    -0.060, -0.059, -0.058, -0.057, -0.056, -0.055, -0.054, -0.053, -0.052, -0.051,
    -0.050, -0.049, -0.048, -0.047, -0.046, -0.045, -0.044, -0.043, -0.042, -0.041,
    -0.040, -0.039, -0.038, -0.037, -0.036, -0.035, -0.034, -0.033, -0.032, -0.031,
    -0.030, -0.029, -0.028, -0.027, -0.026, -0.025, -0.024, -0.023, -0.022, -0.021,
    -0.020, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.013, -0.012, -0.011,
    -0.010, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001,
     0.000,  0.001,  0.002,  0.003,  0.004,  0.005,  0.006,  0.007,  0.008,  0.009,
     0.010,  0.011,  0.012,  0.013,  0.014,  0.015,  0.016,  0.017,  0.018,  0.019,
     0.020,  0.021,  0.022,  0.023,  0.024,  0.025,  0.026,  0.027,  0.028,  0.029,
     0.030,  0.031,  0.032,  0.033,  0.034,  0.035,  0.036,  0.037,  0.038,  0.039,
     0.040,  0.041,  0.042,  0.043,  0.044,  0.045,  0.046,  0.047,  0.048,  0.049,
     0.050,  0.051,  0.052,  0.053,  0.054,  0.055,  0.056,  0.057,  0.058,  0.059,
     0.060,  0.061,  0.062,  0.063,  0.064,  0.065,  0.066,  0.067,  0.068,  0.069,
     0.070,  0.071,  0.072,  0.073,  0.074,  0.075,  0.076,  0.077,  0.078,  0.079,
     0.080
])

energy_values = np.array([
    -1678.6897, -1681.1776, -1683.6237, -1686.0285, -1688.3921, -1690.7151, -1692.9976, -1695.2400, -1697.4427, -1699.6060,
    -1701.7301, -1703.8159, -1705.8635, -1707.8725, -1709.8424, -1711.7725, -1713.6646, -1715.5229, -1717.3517, -1719.1617,
    -1720.9381, -1722.6787, -1724.3837, -1726.0535, -1727.6884, -1729.2886, -1730.8544, -1732.3862, -1733.8841, -1735.3485,
    -1736.7786, -1738.1759, -1739.5441, -1740.8856, -1742.1881, -1743.4359, -1744.6249, -1745.7696, -1746.8851, -1747.9742,
    -1749.0341, -1750.0631, -1751.0623, -1752.0320, -1752.9724, -1753.8837, -1754.7664, -1755.6209, -1756.4445, -1757.2392,
    -1758.0249, -1758.7947, -1759.5368, -1760.2515, -1760.9389, -1761.5993, -1762.2329, -1762.8399, -1763.4207, -1763.9752,
    -1764.5039, -1765.0069, -1765.4844, -1765.9366, -1766.3638, -1766.7661, -1767.1438, -1767.4971, -1767.8261, -1768.1311,
    -1768.4122, -1768.6697, -1768.9038, -1769.1146, -1769.3024, -1769.4674, -1769.6096, -1769.7294, -1769.8269, -1769.9023,
    -1769.9558, -1769.9876, -1769.9978, -1769.9867, -1769.9544, -1769.9010, -1769.8268, -1769.7320, -1769.6167, -1769.4811,
    -1769.3253, -1769.1496, -1768.9541, -1768.7389, -1768.5043, -1768.2503, -1767.9773, -1767.6852, -1767.3744, -1767.0449,
    -1766.6969, -1766.3305, -1765.9460, -1765.5435, -1765.1230, -1764.6849, -1764.2292, -1763.7561, -1763.2657, -1762.7581,
    -1762.2336, -1761.6923, -1761.1343, -1760.5598, -1759.9688, -1759.3616, -1758.7383, -1758.0990, -1757.4439, -1756.7731,
    -1756.0867, -1755.3849, -1754.6678, -1753.9355, -1753.1882, -1752.4260, -1751.6490, -1750.8574, -1750.0513, -1749.2308,
    -1748.3960, -1747.5471, -1746.6842, -1745.8074, -1744.9169, -1744.0127, -1743.0950, -1742.1639, -1741.2195, -1740.2619,
    -1739.2913, -1738.3078, -1737.3115, -1736.3024, -1735.2808, -1734.2467, -1733.2003, -1732.1416, -1731.0707, -1729.9878,
    -1728.8930, -1727.7864, -1726.6681, -1725.5382, -1724.3968, -1723.2440, -1722.0799, -1720.9047, -1719.7183, -1718.5210,
    -1717.3128
])

# === Step 2: Normalize Energy ===
energy_values -= min(energy_values)  # Set minimum energy to 0

# === Step 3: Define Quadratic Model ===
def quadratic(eps, p0, p1, p2):
    return p0 + p1 * eps + p2 * eps**2

# === Step 4: Fit Only in a Small Range Around 0 ===
fit_range = 0.015  # ±0.015 strain range
mask = np.abs(strain_values) <= fit_range
fit_strain = strain_values[mask]
fit_energy = energy_values[mask]

# Fit the quadratic model
popt, _ = curve_fit(quadratic, fit_strain, fit_energy)
p0, p1, p2 = popt

print("Quadratic Coefficients from ±0.015 range:")
print(f"p0 = {p0:.6f},  p1 = {p1:.6f},  p2 = {p2:.6f}")

# === Step 5: Create Fit Curve Over Wider Range for Visualization ===
strain_fit = np.linspace(min(strain_values), max(strain_values), 300)
energy_fit = quadratic(strain_fit, *popt)

# === Step 6: Plot Full Curve ===
plt.figure(figsize=(10, 6))
plt.plot(strain_values, energy_values, 'ko', markersize=3, label="LAMMPS Data")
plt.plot(strain_fit, energy_fit, 'r--', linewidth=2, label="Quadratic Fit")
plt.xlabel("Average Strain (ε)")
plt.ylabel("Strain Energy Density")
plt.title("Strain Energy Density vs Strain with Quadratic Fit")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("task36_quadratic_fullrange.png", dpi=300)
plt.show()

# === Step 7: Plot Zoomed-In View to Identify Valid Quadratic Region ===
zoom_range = 0.04
zoom_mask = np.abs(strain_values) <= zoom_range
zoom_strain = strain_values[zoom_mask]
zoom_energy = energy_values[zoom_mask]

plt.figure(figsize=(10, 6))
plt.scatter(zoom_strain, zoom_energy, color='black', s=25, label="LAMMPS Data")
plt.plot(strain_fit, energy_fit, 'r--', linewidth=2, label="Quadratic Fit")
plt.xlabel("Average Strain (ε)")
plt.ylabel("Strain Energy Density")
plt.title("Zoomed-In View: Validity of Quadratic Assumption")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("task36_quadratic_zoomed.png", dpi=300)
plt.show()

# === Step 8: Suggest Valid Quadratic Region ===
print("\n➡️ From the zoomed-in plot, the quadratic fit closely matches the LAMMPS data roughly within the range ±0.02 to ±0.025.")
print("➡️ Beyond that, deviations become visually significant, so the quadratic assumption is less accurate.")
