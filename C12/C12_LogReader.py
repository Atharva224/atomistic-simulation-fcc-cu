import os
import re

# 📁 Set folder path to your log files
folder_path = r"C:\Users\Atharva\Downloads\Practicals\Atomistic Practical\Task_3_2"

strain_values = []
e_pair_values = []

# Pattern for filenames like log_C12_strain_m0.010.lammps
pattern = re.compile(r"log_C12_strain_([mp]\d+\.\d+)\.lammps")

def convert_label_to_strain(label):
    sign = -1.0 if label.startswith("m") else 1.0
    return sign * float(label[1:])

for filename in os.listdir(folder_path):
    match = pattern.match(filename)
    if match:
        label = match.group(1)
        strain = convert_label_to_strain(label)
        full_path = os.path.join(folder_path, filename)

        with open(full_path, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if "Step" in line and "E_pair" in line:
                    if i + 1 < len(lines):
                        values = lines[i + 1].strip().split()
                        try:
                            e_pair = float(values[2])
                            strain_values.append(strain)
                            e_pair_values.append(e_pair)
                        except (IndexError, ValueError):
                            print(f"❌ Could not extract E_pair from {filename}")
                    break

# Sort values before saving
strain_values, e_pair_values = zip(*sorted(zip(strain_values, e_pair_values)))

# Print values
for s, e in zip(strain_values, e_pair_values):
    print(f"Strain: {s:+.3f} -> E_pair: {e:.6f} eV")

# Save to file
output_file = os.path.join(folder_path, "C12_strain_vs_epair.txt")
with open(output_file, "w") as f:
    f.write("Strain E_pair(eV)\n")
    for s, e in zip(strain_values, e_pair_values):
        f.write(f"{s:+.3f} {e:.6f}\n")

print(f"\n✅ Results saved to: {output_file}")
