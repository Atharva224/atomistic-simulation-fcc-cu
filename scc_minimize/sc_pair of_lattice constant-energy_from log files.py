import os
import re

# ✅ Use raw string to avoid unicode errors in Windows paths
folder_path = r"C:\Users\Atharva\Downloads\Practicals\Atomistic Practical\Task_3_2"

lattice_constants = []
e_pair_values = []

# Loop through all log files in the directory
for filename in os.listdir(folder_path):
    if filename.startswith("log_sc_") and filename.endswith(".lammps"):
        # Extract the lattice constant from the filename
        match = re.search(r"log_sc_(\d+\.\d{1,5})\.lammps", filename)

        if match:
            lattice_const = float(match.group(1))
            full_path = os.path.join(folder_path, filename)

            with open(full_path, 'r') as file:
                lines = file.readlines()
                for i, line in enumerate(lines):
                    if "Step" in line and "E_pair" in line:
                        # Read the next line after header
                        if i + 1 < len(lines):
                            data_line = lines[i + 1].strip()
                            values = data_line.split()
                            try:
                                e_pair = float(values[2])  # E_pair is the 3rd column (index 2)
                                lattice_constants.append(lattice_const)
                                e_pair_values.append(e_pair)
                            except (IndexError, ValueError):
                                print(f"Could not extract E_pair from {filename}")
                        break

# ✅ Display results
for a, e in zip(lattice_constants, e_pair_values):
    print(f"Lattice constant: {a:.5f} Å -> E_pair: {e:.6f} eV")

# ✅ Save results to file
output_file = os.path.join(folder_path, "sc_lattice_vs_epair.txt")
with open(output_file, "w") as f:
    f.write("LatticeConstant(EnergyUnits) E_pair(eV)\n")
    for a, e in zip(lattice_constants, e_pair_values):
        f.write(f"{a:.5f} {e:.6f}\n")

print(f"\n✅ Results saved to: {output_file}")