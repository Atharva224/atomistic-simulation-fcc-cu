# atomistic-simulation-fcc-cu
Atomistic simulation of FCC copper using LAMMPS to compute elastic constants (C11, C12, C44) via energy-strain analysis.



# Atomistic Simulation of FCC Copper Elastic Constants

This repository presents a complete workflow for simulating and analyzing the **elastic constants (C11, C12, C44)** of **FCC copper** using **atomistic simulations**. The simulations are run using **LAMMPS**, and the data is processed and analyzed with **Python** to extract mechanical properties via energy-strain fitting.

---

## 📁 Repository Structure

atomistic-simulation-fcc-copper-elastic-constants/
│
├── C11/ # Uniaxial strain simulation for C11
│ ├── C11_MultiStrain.py # Generates LAMMPS inputs for multiple strains
│ ├── C11_LogReader.py # Extracts potential energy from LAMMPS logs
│ ├── C11.py # Fits energy-strain curve and computes C11
│ ├── fcc_elastic_C11.in # Single-strain LAMMPS input (example)
│ └── Cu_u3.eam # EAM potential file for Copper (add manually)
│
├── report/
│ ├── Atomistic_Practical_report.pdf # Final report with methods and results
│ └── Practical 4 - Atomistic.pdf # Original exercise instructions
│
├── README.md
└── .gitignore




---

## 🧪 Project Overview

- **Simulation Tool**: [LAMMPS](https://www.lammps.org)
- **Postprocessing**: Python (`numpy`, `matplotlib`, `scipy`)
- **Material**: FCC Copper (Cu)
- **Objective**: Compute C11, C12, and C44 from energy-strain relationships

---

## ⚙️ Getting Started

### ✅ Prerequisites

- [Python 3.7+](https://www.python.org/)
- [LAMMPS installed and in your PATH](https://lammps.org/)
- Python packages:
  ```bash
  pip install numpy matplotlib scipy

  
🚀 How to Run
Step 1: Generate strained LAMMPS inputs and run simulations
cd C11
python C11_MultiStrain.py
This creates .in files for different strain values and automatically runs LAMMPS for each.

Step 2: Parse simulation logs to extract energy values
python C11_LogReader.py
This script will create a text file with strain vs. potential energy (C12_strain_vs_epair.txt).

Step 3: Fit energy-strain data to get C11
python C11.py
Outputs the cubic fit and calculates the elastic modulus from second derivative.

📊 Results
The computed value of C11 (and similarly C12, C44 in other folders) is compared against literature values at 0 K. Energy-strain data follows a quadratic trend in the small-strain regime (±0.02), validating the linear elasticity assumption.

📚 References
Zaiser, M. Atomistic Simulation Practical, FAU WW8.

Jacobsen, K. W. et al., Effective Medium Theory, Surf. Sci. 366, 394 (1996).

Overton & Gaffney, Elastic Constants of Copper, Phys. Rev. 98, 969 (1955).

LAMMPS Documentation

🧾 License
This project is provided for educational and academic use only.

👨‍💻 Author
Atharva Sinnarkar
M.Sc. Computational Engineering, FAU Erlangen-Nürnberg
📧 atharvasinnarkar@gmail.com
