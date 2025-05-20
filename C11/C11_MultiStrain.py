import os
import subprocess
import numpy as np

# Path to LAMMPS executable
LAMMPS_EXE = "lmp"  

# Lattice constant (optimized)
lattice_constant = 3.61

# Strain range for uniaxial deformation
strain_values = np.arange(-0.080, 0.0801, 0.001)

for strain in strain_values:
    strain_label = f"{strain:+.3f}".replace("+", "p").replace("-", "m")
    input_filename = f"uniaxial_strain_C11_{strain_label}.in"
    log_filename = f"log_C11_uniaxial_{strain_label}.lammps"

    deform_xx = 1.0 + strain
    deform_yy = 1.0  # No strain in y
    deform_zz = 1.0  # No strain in z

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
