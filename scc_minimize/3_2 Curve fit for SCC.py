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

#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.
