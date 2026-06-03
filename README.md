# Taloside Screening Pipeline – Phase 2 & 3

**Virtual library generation, drug-likeness filtering, PAINS screening, and AutoDock Vina docking for taloside-triazole derivatives against Galectin-3 (PDB: 3ZSJ).**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://python.org)
[![RDKit](https://img.shields.io/badge/RDKit-2022.09.1+-green.svg)](https://rdkit.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://doi.org/10.5281/zenodo.20476421)

---

## Overview

This pipeline automates the generation and prioritisation of taloside-triazole
derivatives as Galectin-3 (Gal-3) ligand candidates. Gal-3 is a β-galactoside-
binding lectin overexpressed in cancer, fibrosis, and inflammation, and is
currently in clinical evaluation as a drug target. Starting from a validated
taloside azide scaffold, the pipeline:

1. **Generates a virtual library** via CuAAC and RuAAC click chemistry
   (8 building blocks × 2 regioisomers = 16 geometrically verified compounds)
2. **Filters for drug-likeness** using carbohydrate-adjusted Lipinski thresholds
3. **Screens for PAINS** using RDKit's validated FilterCatalog (PAINS_A/B/C)
4. **Ranks compounds** by a composite lead score incorporating physicochemical
   descriptors and a synthetic-accessibility term
5. **Docks all leads** against the 3ZSJ Gal-3 CRD crystal structure using
   AutoDock Vina with dynamic grid centring on the crystal lactose pose

---

## Key Results

| Metric | Value |
|---|---|
| **Compounds generated** | 16 (8 BBs × 2 regioisomers, all geometry-verified) |
| **Lipinski pass rate** | 87.5% (14 / 16) |
| **PAINS flagged** | 0 — all 14 clean |
| **Lead score range** | 0.399 – 0.905 |
| **Docking success** | 14 / 14 (100%) |
| **Vina score range** | −4.98 to −6.14 kcal/mol |
| **Top hit (combined score)** | SCAF-001_BB-004-4F_CuAAC_1 (0.922) |

MW range: 490–579 Da · LogP range: 0.54–1.91 · TPSA: 168–211 Å²

The two 2-nitrophenyl analogues fail Lipinski (HBA = 13 > threshold 12) and
are retained in `03_lipinski_failed.csv` for reference.

---

## Installation

### Python dependencies

```bash
git clone https://github.com/adamholohan6/taloside-screening-pipeline.git
cd taloside-screening-pipeline

python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

pip install -r requirements.txt
pip install -e .
```

For development and testing:

```bash
pip install -r requirements-dev.txt
pip install -e .
```

### External dependencies (Phase 3)

Phase 3 docking requires two executables on `PATH`:

| Tool | Purpose | Download |
|---|---|---|
| **AutoDock Vina** | Molecular docking | https://vina.scripps.edu/downloads/ |
| **Open Babel 3.x** | Ligand PDBQT preparation (Gasteiger charges) | https://openbabel.org/wiki/Get_Open_Babel |

Verify both are available before running Phase 3:

```bash
obabel -V        # should print: Open Babel 3.x.x
vina --version   # should print version info
```

---

## Running the pipeline

### Phase 2 — Library generation and lead scoring

```bash
python -m taloside_pipeline.phase2_integration
```

Outputs are written to `phase2_output/`:

| File | Content |
|---|---|
| `01_all_generated_compounds.csv` | All 16 geometry-verified compounds with descriptors |
| `02_lipinski_passed.csv` | 14 compounds passing carbohydrate-adjusted Lipinski |
| `03_lipinski_failed.csv` | 2 compounds failing (2-nitrophenyl series, HBA > 12) |
| `04_lipinski_clean_no_pains.csv` | 14 Lipinski-passed, PAINS-clean compounds |
| `05_lipinski_with_pains.csv` | PAINS-flagged compounds (empty for this library) |
| `06_pains_undetermined.csv` | Compounds with undetermined PAINS status (empty) |
| `07_lead_scored.csv` | 14 compounds ranked by composite lead score |

### Phase 3 — AutoDock Vina docking

```python
from pathlib import Path
from taloside_pipeline.phase3_docking import DockingConfig, run_phase3_pipeline

cfg = DockingConfig(
    receptor_pdb=Path("data/docking/3ZSJ.pdb"),
    receptor_pdbqt=Path("data/docking/3ZSJ_clean.pdbqt"),
    vina_executable="vina",          # Windows: r"C:\path\to\vina.exe"
    lead_csv=Path("phase2_output/07_lead_scored.csv"),
)
run_phase3_pipeline(cfg, validate=True)
```

Output is written to `phase3_output/08_docking_results.csv`. Running with
`validate=True` auto-centres the docking grid on the crystal lactose pose
(BGC+GAL residues from 3ZSJ.pdb) and performs a validation redock before
screening the lead compounds.

### Interactive notebook

```bash
jupyter notebook notebooks/Phase2_VirtualLibraryExpansion.ipynb
```

---

## Project structure

```text
taloside-screening-pipeline/
├── src/taloside_pipeline/
│   ├── __init__.py
│   ├── glycolibrary_generator.py   # SMARTS-based combinatorial library engine
│   ├── phase2_integration.py       # Lipinski → PAINS → lead scoring workflow
│   ├── phase3_docking.py           # 3D prep, Open Babel, AutoDock Vina docking
│   ├── descriptor_calculator.py    # Standalone descriptor utilities
│   └── library_generator.py       # Compatibility wrapper
├── notebooks/
│   └── Phase2_VirtualLibraryExpansion.ipynb
├── tests/                          # 53 unit tests (pytest)
├── data/docking/                   # 3ZSJ.pdb, 3ZSJ_clean.pdbqt (not committed)
├── CHANGELOG.md
├── requirements.txt
├── requirements-dev.txt
├── setup.py
└── README.md
```

---

## Testing

```bash
pytest                  # run all 53 tests
pytest -m unit          # unit tests only (no external tools required)
```

All tests are self-contained. Phase 3 unit tests mock the Vina subprocess;
only `test_embed_and_pdbqt_benzene` and `test_ligand_pdbqt_has_charges_and_torsions`
require Open Babel to be on PATH.

---

## Limitations

- **Docking scores are not binding affinities.** AutoDock Vina scores are
  used for within-library ranking only. No experimental K_d or IC₅₀ data are
  available for these compounds.
- **TPSA > 140 Å² for all compounds.** The Veber oral-bioavailability criterion
  is not met, which is expected and acceptable for an extracellular target.
- **Library size.** 16 compounds is a proof-of-concept enumeration. No
  structure–activity relationships can be established at this scale.

---

## References

- Lipinski, C. A. et al. (1997). *Adv. Drug Deliv. Rev.* — Rule of 5
- Baell, J. B. & Holloway, G. A. (2010). *J. Med. Chem.* — PAINS substructures
- Eberhardt, J. et al. (2021). *J. Chem. Inf. Model.* — AutoDock Vina 1.2
- Macedo, J. N. A. et al. (2025). *Mol. Med.* — Galectin-3 inhibitor review
- Öberg, C. T. et al. (2011). *ChemMedChem* — Taloside inhibitors of Gal-1/3
- RDKit: https://www.rdkit.org/

---

## Acknowledgements

This work builds on MSc research conducted at the University of Galway under
the supervision of Prof. Helen Blanchard and Dr. Chandan Kishor.

---

## Contact

Adam Holohan — [GitHub](https://github.com/adamholohan6)

For questions or suggestions, please open an issue on GitHub.

---

## How to cite

Holohan, A. (2026). Taloside Screening Pipeline (Version 1.0). Zenodo.
https://doi.org/10.5281/zenodo.20476422
