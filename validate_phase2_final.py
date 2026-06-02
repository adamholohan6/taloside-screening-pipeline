"""
Final validation of Phase 2 products after targeted charge normalization.
"""

import pandas as pd
from rdkit import Chem

# Load regenerated products
df = pd.read_csv("phase2_output/01_all_generated_compounds.csv")

print("=" * 80)
print("PHASE 2 VALIDATION AFTER TARGETED CHARGE NORMALIZATION")
print("=" * 80)
print("Total products:", len(df))
print("=" * 80)

# Check for [N+] or [N-] in SMILES
charged_smiles = []
for i, row in df.iterrows():
    smiles = row["product_smiles"]
    compound_id = row["compound_id"]
    
    if "[N+]" in smiles or "[N-]" in smiles:
        charged_smiles.append((compound_id, smiles))

print("\nSMILES with charged nitrogens:", len(charged_smiles))
if charged_smiles:
    print("\nCharged SMILES found:")
    for compound_id, smiles in charged_smiles:
        print("  {}: {}".format(compound_id, smiles))
else:
    print("No charged nitrogens found in SMILES")

# MolFromSmiles validation
print("\n" + "=" * 80)
print("MOLFROMSMILES VALIDATION")
print("=" * 80)

failures = []

for i, row in df.iterrows():
    smiles = row["product_smiles"]
    compound_id = row["compound_id"]
    
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        failures.append((compound_id, smiles))
        print("FAIL:", compound_id, "-", smiles)
    else:
        print("PASS:", compound_id)

print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)
print("Total products:", len(df))
print("Successfully parsed:", len(df) - len(failures))
print("Failed:", len(failures))

if failures:
    print("\nFAILING SMILES:")
    for compound_id, smiles in failures:
        print("  {}: {}".format(compound_id, smiles))
else:
    print("\nAll products successfully parsed!")

print("=" * 80)
