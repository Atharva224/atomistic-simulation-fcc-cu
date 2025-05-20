import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Path to LAMMPS executable
LAMMPS_EXE = "lmp"  # Change if needed (e.g., "lmp_serial")

# Lattice constant range (in Angstroms)
lattice_constants = np.arange(3.40, 3.71, 0.01)
energies = []

# Loop over lattice constants
for a in lattice_constants:
    input_filename = f"input_{a:.2f}.in"
    log_filename = f"log_{a:.2f}.lammps"

    # Write input file
    with open(input_filename, "w") as f:
        f.write(f"""
units metal
dimension 3
boundary p p p
atom_style atomic

#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.


