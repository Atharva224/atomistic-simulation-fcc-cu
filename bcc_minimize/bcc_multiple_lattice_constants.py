import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Path to LAMMPS executable
LAMMPS_EXE = "lmp"  # Change if needed (e.g., "lmp_serial")

# Lattice constant range (in Angstroms)
lattice_constants = np.arange(2.80, 3.21, 0.01)  # Typical BCC Cu range
energies = []

# Loop over lattice constants
for a in lattice_constants:
    input_filename = f"bcc_input_{a:.2f}.in"
    log_filename = f"log_bcc_{a:.2f}.lammps"

    # Write input file
    with open(input_filename, "w") as f:
        f.write(f"""
        
#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.
