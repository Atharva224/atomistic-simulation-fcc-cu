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

#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.
