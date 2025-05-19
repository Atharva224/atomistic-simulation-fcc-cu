
# Atomistic Simulation of FCC Copper Elastic Constants

This repository presents a complete workflow for simulating and analyzing the **elastic constants (C11, C12, C44)** of **FCC copper** using **atomistic simulations**. The simulations are run using **LAMMPS**, and the data is postprocessed with **Python** to extract mechanical properties via energy-strain fitting.

---

## 📁 Repository Structure

```text
atomistic-simulation-fcc-copper-elastic-constants/
│
├── Atomistic_Practical_report.pdf      # Final report with methodology and results
├── Practical 4 - Atomistic.pdf         # Original exercise instructions
│
├── C11/                                # Folder for computing elastic constant C11
├── C12/                                # Folder for computing elastic constant C12
├── C44/                                # Folder for computing elastic constant C44
│
├── fcc_minimize/                       # FCC lattice energy minimization scripts
├── bcc_minimize/                       # BCC lattice energy minimization scripts
├── scc_minimize/                       # SCC lattice energy minimization scripts
│
├── Task 3.6/                           # Optional: Quadratic strain energy validation
│
└── README.md


📘 Project Overview
Simulation Tool: LAMMPS

Postprocessing: Python (numpy, matplotlib, scipy)

Material: FCC Copper (Cu)

Objective: Compute C11, C12, and C44 from energy-strain relationships

✅ Getting Started
Prerequisites
Python 3.7+
LAMMPS installed and in your PATH

Python packages:
pip install numpy matplotlib scipy
