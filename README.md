# Molecular Coarse-Graining and ANM Analysis

This repository contains scripts for performing **High&Low Resolution (HLR) processing**, **Mixed Coarse-Graining (MCG)**, **Anisotropic Network Model (ANM) calculations**, and **B-factor analysis** for protein structures.

## 📌 Overview
These scripts facilitate the analysis of protein structures by:
- Converting **full-atom** representation into **mixed coarse-grained** model.
- Computing **normal modes** of **mixed coase-grained** model using ANM.
- Comparing **theoretical and experimental B-factors** to evaluate model accuracy.

## 🛠 Dependencies
- **Python 3.x** (for `HLR.py`)
- **MATLAB** (for `mcg.mlx`, `ANM.mlx`, `bfac.mlx`)
- **`sign_flip.m`** (required for MCG)
- **PDB file** as input

---

## 🚀 Usage

### **1️. Hierarchical Low-Resolution (HLR) Processing**
Assigns **coarse-grained (CG)** or **full-atom (FA)** labels to atoms within a specified **cutoff distance**.

**Command-Line Usage:**
```bash
python HLR.py input.pdb HLR.out A 165 CG 30
```
**Parameters:**
- `input.pdb` → Input PDB file.
- `HLR.out` → Processed PDB data with CG/FA assignments.
- `A` → Spherical high-resolution region center chain ID.
- `165` → Spherical high-resolution region center residue ID.
- `FA` → Label for the selected atoms (full-atom FA as in the PDB or coarse-grained CG).
- `30` → Spherical high-resolution region radius (Å).

**Output:**
- **PDB data** with updated CG/FA labels.
- The CG/FA labels in the output file can be manually modified if needed.
    Adjustments can be made to short segments or specific residues based on analysis requirements.
    Modifications should be done carefully to ensure the calculations remain consistent.

---

### **2️. Mixed Coarse-Graining (MCG)**
Processes the **HLR output**, generating new **CG and FA nodes** for further analysis.

**Usage:**
- Run `mcg.mlx` in **MATLAB**.
- Requires **HLR output**.
- Uses `sign_flip.m` for mass matrix correction.

**Output:**
- `new_coords.pdb` → Updated mixed coarse-grained PDB.
- `HES.dat` → Hessian matrix for ANM analysis.
- `FA_new.out`, `CG_new.out` → FA and CG node assignments.
- `mwt.out` → Weight of each node (residue) taken as the number of atoms."

---

### **3️. Anisotropic Network Model (ANM) Analysis**
Computes **ANM normal modes** from the Hessian matrix.

**Usage:**
- Run `ANM.mlx` in **MATLAB**.
- Requires: `HES.dat`, `new_coords.pdb`, `mwt.out`, `FA_new.out`, `CG_new.out`.
- The **`C` parameter** (**default = 50**) controls the **deformation scale**.

**Output:**
- `eigenvalues.txt` → ANM **eigenvalues**.
- `eigenvectors.txt` → First **10 eigenvectors**.
- `mXa/b.pdb` → Mode shapes at slow mode X.
- `mXa/b_hr.pdb` → High-resolution region at slow mode X.
- `mXa/b_lr.pdb` → Low-resolution region at slow mode X.

**Chain-Specific Mode Deformations**
- The ANM script generates two deformed structures per mode:
    "a" structure represents the positive displacement along the normal mode.
    "b" structure represents the negative displacement along the normal mode.
- For each chain present in the structure, a separate deformed structure is written.
    The script automatically detects the chains (e.g., A, B, C, etc.) and assigns deformations accordingly.
    If a structure contains multiple chains, the chain identifiers must be properly adjusted within the script to ensure correct labeling in the output files.
    Users may need to manually modify chain names in the script if working with structures that deviate from the standard A/B labeling.

Ensuring proper chain assignment is critical to maintaining structural integrity when interpreting mode-based deformation

---

### **4️. B-factor Calculation & Plotting**
Compares **theoretical vs experimental B-factors** using normal modes from ANM.

**Usage:**
- Run `bfac.mlx` in **MATLAB**.
- Requires: `HES.dat`, `new_coords.pdb`.

**Output:**
- `Bfacs.png` → **Experimental vs theoretical B-factor plot**.
- Pearson & Spearman **correlation coefficients**.

---

## 🔍 **Example Input & Output Files**
The `example_data/` directory includes:
- **TIM Interface**:
  - `HLR/` → HLR input & output.
  - `mcg/` → MCG input & output.
  - `ANM_Bfc/` → ANM & Theoretical B-factor calculation input & output.

---

## 📜 **License**
This project is licensed under the **MIT License**. See the `LICENSE` file for details.

---

## 📌 **Citation**
If you use this work in your research, please cite it appropriately.

---


