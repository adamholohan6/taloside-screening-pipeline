# Taloside Screening Pipeline

A Python-based computational chemistry pipeline for calculating physicochemical descriptors of taloside compounds. This tool automates the calculation of ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) properties using RDKit, supporting structure-activity relationship (SAR) analysis.

## 📌 Overview

Taloidases are complex carbohydrates with significant pharmaceutical potential. This pipeline:
- Converts SMILES strings to 3D molecular structures
- Calculates key drug-likeness descriptors
- Exports results to CSV for downstream analysis
- Provides data for Lipinski's Rule of Five screening

## ✨ Features

- **Batch Processing:** Calculate descriptors for 10+ taloside variants
- **ADMET Calculations:**
  - Molecular Weight (MW)
  - Partition Coefficient (LogP)
  - Hydrogen Bond Donors/Acceptors
  - Topological Polar Surface Area (TPSA)
  - Rotatable Bonds
- **Error Handling:** Gracefully skips invalid SMILES strings
- **CSV Export:** Automated output for statistical analysis
- **CLI Support:** Command-line interface with flexible output options

## 🛠️ Installation

### Requirements
- **Python 3.6+**
- **RDKit** (cheminformatics toolkit)
- **Pandas** (data manipulation)

### Setup

1. **Clone the repository:**
```bash
git clone https://github.com/adamholohan6/taloside-screening-pipeline.git
cd taloside-screening-pipeline
```

2. **Install dependencies:**
```bash
pip install -r requirements.txt
```

Or manually install:
```bash
pip install rdkit pandas
```

## 🚀 Quick Start

### Run with Defaults
```bash
python generate_library.py
```

This generates `taloside_descriptors.csv` with all 10 taloside compounds.

### Custom Output File
```bash
python generate_library.py -o my_descriptors.csv
```

### View Help
```bash
python generate_library.py -h
```

## 📊 Example Output

The script generates a CSV with the following structure:

| Compound | SMILES | Molecular_Weight | LogP | H_Donors | H_Acceptors | TPSA | Rotatable_Bonds |
|----------|--------|------------------|------|----------|-------------|------|------------------|
| Talo-1 | O=C(O[C@H]1... | 530.49 | 1.15 | 2 | 12 | 177.53 | 10 |
| Talo-2 | O=C(O[C@H]1... | 530.49 | 1.15 | 2 | 12 | 177.53 | 10 |

## 📋 Descriptor Reference

| Descriptor | Interpretation | Drug-like Range |
|---|---|---|
| **MW** | Molecular weight | < 500 g/mol |
| **LogP** | Lipophilicity | 0 - 5 |
| **H-Donors** | H-bond donors | ≤ 5 |
| **H-Acceptors** | H-bond acceptors | ≤ 10 |
| **TPSA** | Polar surface area | 20 - 130 Ų |
| **Rotatable Bonds** | Molecular flexibility | ≤ 10 |

## 🔧 Customization

### Add New Compounds

Edit `generate_library.py` and add to the `load_compounds_from_dict()` function:

```python
"Talo-11": "YOUR_SMILES_STRING_HERE"
```

### Load from External File

Modify the script to read SMILES from a CSV:

```python
df_input = pd.read_csv("compounds.csv")
compounds = dict(zip(df_input['Name'], df_input['SMILES']))
```

## 📁 Project Structure

```
taloside-screening-pipeline/
├── generate_library.py           # Main calculation script
├── taloside_descriptors.csv      # Output data
├── requirements.txt              # Python dependencies
├── README.md                     # This file
└── CONTRIBUTING.md              # Contribution guidelines
```

## 🧪 Testing

Run the script on the included compound library:

```bash
python generate_library.py -o test_output.csv
```

Check the output:
```bash
cat taloside_descriptors.csv
```

Expected: 10 compounds with calculated descriptors.

## 📈 Future Development

- [ ] **Phase 1:** ADMET descriptor calculation (Current)
- [ ] **Phase 2:** Docking integration with AutoDock Vina
- [ ] **Phase 3:** QSAR modeling for binding affinity prediction
- [ ] **Phase 4:** Virtual library expansion (50-100 variants)
- [ ] **Phase 5:** R-based statistical analysis and visualization

## 🤝 Contributing

Contributions welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📄 License

This project is provided as-is for research and educational purposes.

## 📧 Questions?

Open a GitHub issue for questions, bug reports, or feature requests.

## 🔬 References

- **RDKit:** https://www.rdkit.org/
- **Lipinski's Rule of Five:** Lipinski et al., Adv. Drug Deliv. Rev. (1997)
- **SMILES Notation:** https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
