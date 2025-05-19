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

variable a equal {a}
lattice fcc ${{a}}
region box block 0 5 0 5 0 5 units lattice
create_box 1 box
create_atoms 1 box

pair_style eam
pair_coeff * * Cu_u3.eam

fix 1 all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1.0e-8 1.0e-10 5000 10000

thermo_style custom step pe
thermo 1
print "Final Energy for a = ${{a}} : ${{pe}}"

write_data data_{a:.2f}.lmp
""")

    # Run LAMMPS
    subprocess.run([LAMMPS_EXE, "-in", input_filename], stdout=open(log_filename, "w"))

    # Extract final energy from log file
    with open(log_filename) as f:
        for line in f:
            if f"Final Energy for a = {a:.2f}" in line:
                energy = float(line.strip().split()[-1])
                energies.append(energy)
                print(f"a = {a:.2f} Å --> Energy = {energy:.6f} eV")
                break

# Save to file
with open("energy_vs_lattice.txt", "w") as f:
    f.write("LatticeConstant(EnergyUnits) PotentialEnergy(eV)\n")
    for a, e in zip(lattice_constants, energies):
        f.write(f"{a:.3f} {e:.6f}\n")


