import os
import subprocess
import numpy as np

# --- Configuration ---
LAMMPS_EXE = "lmp"  
lattice_constant = 3.61
strain_values = np.arange(-0.08, 0.081, 0.001)

# --- Loop to create and run LAMMPS input for each shear strain ---
for strain in strain_values:
    strain_label = f"{strain:+.3f}".replace("+", "p").replace("-", "m")
    input_filename = f"strain_C44_{strain_label}.in"
    log_filename = f"log_C44_strain_{strain_label}.lammps"

    with open(input_filename, "w") as f:
        f.write(f"""
units metal
dimension 3
boundary p p p
atom_style atomic

variable a equal {lattice_constant}
lattice fcc ${{a}}
region box block 0 5 0 5 0 5 units lattice
create_box 1 box
create_atoms 1 box

#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the file and mention the proper usecase.

