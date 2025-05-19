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
```

---

## 📘 Project Overview

- **Simulation Tool**: [LAMMPS](https://www.lammps.org)
- **Postprocessing**: Python (`numpy`, `matplotlib`, `scipy`)
- **Material**: FCC Copper (Cu)
- **Objective**: Compute C11, C12, and C44 from energy-strain relationships

---

## ✅ Getting Started

### Prerequisites

- [Python 3.7+](https://www.python.org/)
- [LAMMPS installed and in your PATH](https://lammps.org/)
- Python packages:
  ```bash
  pip install numpy matplotlib scipy
  ```

---

## 🧪 How to Run

Each of the `C11`, `C12`,`C44` and `fcc/bcc/sc minimization`, folders contains:

1. Python script to generate multiple LAMMPS `.in` files with different strains  
2. A LAMMPS loop or batch runner  
3. A log parser script to extract `E_pair` values  
4. Python code to fit energy-strain curves and compute the elastic constants

## 🧪 SCC Copper Lattice Minimization & Fitting

The `scc_minimize` folder includes:

- `sc_minimize.in`: A basic LAMMPS input script to minimize the lattice energy of simple cubic (SC) copper.
- `sc_multiple_lattice_constants.py`: Automates LAMMPS runs across a range of SC lattice constants.
- `sc_pair of_lattice constant-energy_from log files.py`: Extracts `E_pair` values from multiple LAMMPS log files.
- `3_2 Curve fit for SCC.py`: Performs quadratic fitting of potential energy vs. lattice constant to estimate equilibrium values.

### 🔄 Example Workflow:

**1. Run single SC relaxation:**
```bash
lmp -in sc_minimize.in
```

**2. Perform batch runs across lattice constants:**
```bash
python sc_multiple_lattice_constants.py
```

**3. Extract `E_pair` values from log files:**
```bash
python "sc_pair of_lattice constant-energy_from log files.py"
```

**4. Fit quadratic curve and visualize:**
```bash
python "3_2 Curve fit for SCC.py"
```

Outputs include:
- `sc_energy_vs_lattice.txt`: Raw lattice constant vs energy data.
- `sc_quadratic_fit_highres.png`: Plot of fitted curve.


### 🔄 Example for `C11`

**Step 1**: Generate strained LAMMPS inputs and run simulations  
```bash
cd C11
python C11_MultiStrain.py
```

**Step 2**: Parse LAMMPS log files  
```bash
python C11_LogReader.py
```

**Step 3**: Fit energy-strain curve and compute C11  
```bash
python C11.py
```

Repeat similar steps for `C12` and `C44`.

---

## 📊 Results

- Elastic constants are extracted via second derivatives of polynomial fits
- Quadratic energy-strain assumption is validated in `Task 3.6` folder
- Results are compared with literature values at 0 K for verification

---

## 📚 References

- Zaiser, M. *Atomistic Simulation Practical*, FAU WW8  
- Jacobsen, K. W. et al., *Effective Medium Theory*, Surf. Sci. 366, 394 (1996)  
- Overton & Gaffney, *Elastic Constants of Copper*, Phys. Rev. 98, 969 (1955)  
- [LAMMPS Documentation](https://docs.lammps.org)

---

## 🧾 License

This project is provided for educational and academic use only.

---

## 👨‍💻 Author

**Atharva Sinnarkar**  
M.Sc. Computational Engineering, FAU Erlangen-Nürnberg  
📧 [atharvasinnarkar@gmail.com](mailto:atharvasinnarkar@gmail.com)
---
