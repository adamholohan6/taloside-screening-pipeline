"""
Validate regenerated Phase 2 products with Chem.MolFromSmiles
"""

import pandas as pd
from rdkit import Chem

# Load regenerated products
df = pd.read_csv("phase2_output/01_all_generated_compounds.csv")

print("=" * 80)
print("VALIDATING REGENERATED PHASE 2 PRODUCTS")
print("=" * 80)
print("Total products:", len(df))
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
        print(f"  {compound_id}: {smiles}")
else:
    print("\nAll products successfully parsed!")

print("=" * 80)
