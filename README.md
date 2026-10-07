# Automated High-Throughput Molecular Docking & Interaction Pipeline

A fully automated structural biology workflow that executes batch molecular docking, identifies binding cavities, and automatically compiles protein-ligand interaction data into professional, client-ready PDF reports. 

This pipeline integrates open-source chemoinformatics tools with automated bash scripting to eliminate the manual bottleneck of virtual screening, converting raw structural files (`.pdb`, `.sdf`, `.mol2`) into visualized docking profiles in a single command.

## Pipeline Architecture
The workflow handles everything from ligand sanitization to final visual reporting without requiring intermediate manual input.

* **Receptor & Ligand Preparation:** `OpenBabel` (Automatically sanitizes `.pdb` structures, extracts ATOM records, and batch-converts `.sdf`/`.mol2` libraries into 3D `.pdbqt` formats).
* **Automated Cavity Detection:** `fpocket` (Calculates binding pocket volumes and extracts precise geometric center coordinates).
* **Molecular Docking:** `AutoDock Vina` (Executes high-throughput virtual screening using an automated grid box).
* **Interaction Profiling:** `PLIP` (Protein-Ligand Interaction Profiler) & `PyMOL` (Calculates hydrogen bonds, hydrophobic interactions, and pi-stacking, generating 3D `.pse` session files).
* **Report Generation:** `Python (FPDF)` (Parses Vina logs and PLIP text files to generate a unified 2-page PDF summary for every ligand).

## Key Engineering Features
* **Intelligent Pose Extraction:** Prevents catastrophic downstream visualization failures by mathematically extracting only the top-scoring thermodynamic binding pose (`-f 1 -l 1`) before passing the complex to PLIP.
* **Dynamic Ligand Handling:** Automatically detects and converts diverse ligand inputs (`.sdf`, `.mol2`, `.pdb`) instead of hard-failing on non-standard formats.
* **Phantom Pocket Failsafe:** Python-based coordinate validation prevents AutoDock Vina from docking into empty space if `fpocket` fails to find a valid cavity in tightly folded proteins.

## How to Run the Pipeline

### 1. Setup the Environment
Clone this repository and build the dedicated Conda environment to resolve all biophysics dependencies automatically:
```bash
conda env create -f environment_docking.yml
conda activate docking_env
