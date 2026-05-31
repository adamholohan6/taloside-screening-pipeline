# Taloside Screening Pipeline – Phase 2

**Virtual library generation, drug‑likeness filtering, and PAINS screening for taloside derivatives.**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://python.org)
[![RDKit](https://img.shields.io/badge/RDKit-2022.09.1-green.svg)](https://rdkit.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) <!-- Replace with your Zenodo DOI after release -->

---

## 📌 Overview

This pipeline automates the discovery of new taloside‑based drug candidates for galectin‑3 (Gal‑3), a protein implicated in cancer, fibrosis, and inflammation. Starting from a validated taloside scaffold (from previous MSc work), the pipeline:

1. **Generates a virtual library** via click chemistry (8 building blocks × 2 regioisomers = 16 compounds)
2. **Calculates physicochemical descriptors** (MW, LogP, TPSA, H‑bond donors/acceptors, rotatable bonds)
3. **Applies Lipinski’s Rule of 5** with carbohydrate‑adjusted thresholds (MW ≤600, LogP ≤4, HBD ≤6, HBA ≤12)
4. **Screens for PAINS** (Pan‑Assay Interference Compounds) using RDKit’s validated substructure filters

All results are exported as CSV files for downstream analysis, docking, or experimental validation.

---

## ✨ Key Results

| Metric | Value |
|--------|-------|
| **Compounds generated** | 16 (8 building blocks × 2 regioisomers) |
| **Lipinski pass rate** | 100% (carbohydrate‑adjusted) |
| **PAINS flagged** | 0 (all compounds clean) |
| **Lead list** | `phase2_output/04_lipinski_clean_no_pains.csv` |

All compounds have molecular weights between 490–579 Da, LogP 0.02–1.36, and 2 hydrogen‑bond donors – favourable for an extracellular target like galectin‑3.

---

## 🛠️ Installation

### Requirements
- Python 3.8+
- RDKit ≥2022.09.1
- pandas, numpy (for descriptor calculation)

### Setup

```bash
# Clone the repository
git clone https://github.com/adamholohan6/taloside-screening-pipeline.git
cd taloside-screening-pipeline

# (Optional) Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
🚀 Quick Start
Run the full pipeline from the src/taloside_pipeline directory:

bash
cd src/taloside_pipeline
python phase2_integration.py
Output CSV files are saved in phase2_output/:

File	Content
01_all_generated_compounds.csv	All 16 compounds with descriptors
02_lipinski_passed.csv	Drug‑like compounds (all 16)
03_lipinski_failed.csv	(empty – no violations)
04_lipinski_clean_no_pains.csv	High‑priority leads (no PAINS)
05_lipinski_with_pains.csv	(empty – no PAINS detected)
For interactive visualisation (descriptor distributions, building block analysis, etc.), open the Jupyter notebook:

bash
jupyter notebook notebooks/Phase2_VirtualLibraryExpansion.ipynb
📁 Project Structure
text
taloside-screening-pipeline/
├── src/
│   └── taloside_pipeline/
│       ├── phase2_integration.py      # Main pipeline script
│       ├── glycolibrary_generator.py  # Combinatorial library generator
│       ├── descriptor_calculator.py   # Legacy Phase 1 (descriptors only)
│       └── __init__.py
├── notebooks/
│   └── Phase2_VirtualLibraryExpansion.ipynb   # Visualisation notebook
├── phase2_output/                     # Generated CSV files (created at runtime)
├── tests/                             # Unit tests (pytest)
├── requirements.txt
├── LICENSE
└── README.md
🧪 Testing
Run the unit tests from the project root:

bash
pytest tests/
📄 License
This project is licensed under the MIT License – see the LICENSE file for details.

🔬 References
Lipinski, C. A. et al. (1997). Adv. Drug Deliv. Rev. – Rule of 5

Baell, J. B. & Holloway, G. A. (2010). J. Med. Chem. – PAINS substructures

RDKit: https://www.rdkit.org/

🙏 Acknowledgements
This work builds on MSc research conducted at the University of Galway under the supervision of Prof. Helen Blanchard and Dr. Chandan Kishor.

📧 Contact
Adam Holohan – GitHub – ORCID (add your ORCID if you have one)

For questions or suggestions, please open an issue on GitHub.

📌 How to Cite
If you use this pipeline in your research, please cite the Zenodo archive:

Holohan, A. (2026). Taloside Screening Pipeline – Phase 2 (Version 1.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

(Replace XXXXXXX with your actual Zenodo DOI after release.)
