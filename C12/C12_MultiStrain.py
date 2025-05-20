import os
import subprocess
import numpy as np

# Path to LAMMPS executable
LAMMPS_EXE = "lmp"  

# Fixed lattice constant (Angstrom)
lattice_constant = 3.61

# Strain values from -0.08 to +0.08 in steps of 0.001
strain_values = np.arange(-0.080, 0.0801, 0.001)

# Loop over strain values
for strain in strain_values:
    # Label for filenames
    strain_label = f"{strain:+.3f}".replace("+", "p").replace("-", "m")
    input_filename = f"strain_C12_{strain_label}.in"
    log_filename = f"log_C12_strain_{strain_label}.lammps"

    # Biaxial deformation in x and y directions
    deform_xx = 1.0 + strain
    deform_yy = 1.0 + strain

    # Write the LAMMPS input file
    with open(input_filename, "w") as f:
        f.write(f"""
units metal
dimension 3
boundary p p p
atom_style atomic

#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.
