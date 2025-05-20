import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Path to LAMMPS executable
LAMMPS_EXE = "lmp"  # Change if needed (e.g., "lmp_serial")

# Lattice constant range (in Angstroms)
lattice_constants = np.arange(2.36, 2.415, 0.0025)  # SC range for Cu

energies = []

# Loop over lattice constants
for a in lattice_constants:
    input_filename = f"sc_input_{a:.5f}.in"
    log_filename = f"log_sc_{a:.5f}.lammps"

    # Write input file
    with open(input_filename, "w") as f:
        f.write(f""" 
units metal
dimension 3
boundary p p p
atom_style atomic

#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.
