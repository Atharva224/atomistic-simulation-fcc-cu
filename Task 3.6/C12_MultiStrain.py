import os
import subprocess
import numpy as np

# Path to LAMMPS executable
LAMMPS_EXE = "lmp"  # Update this if your LAMMPS command is different

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

variable a equal {lattice_constant}
lattice fcc ${{a}}
region box block 0 5 0 5 0 5 units lattice
create_box 1 box
create_atoms 1 box

pair_style eam
pair_coeff * * Cu_u3.eam

variable strain equal {strain}
variable deform_xx equal {deform_xx}
variable deform_yy equal {deform_yy}

# Apply biaxial strain (z-direction remains unchanged)
change_box all x scale ${{deform_xx}} y scale ${{deform_yy}} remap

fix 1 all box/relax iso 0.1 vmax 0.001
min_style cg
minimize 1.0e-8 1.0e-10 5000 10000

thermo_style custom step pe vol

thermo 1

write_data C12_strained_{strain_label}.lmp
""")

    # Run the simulation
    subprocess.run([LAMMPS_EXE, "-in", input_filename], stdout=open(log_filename, "w"))
