"""
Analyze Phase 2 SMARTS patterns and generated products to determine root cause
"""

# Phase 2 SMARTS patterns from glycolibrary_generator.py
TRIAZOLE_1_4_CuAAC = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"
TRIAZOLE_1_5_RuAAC = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n:3]-[n+:2]=[n-:1]-1"

print("=" * 80)
print("PHASE 2 SMARTS PATTERN ANALYSIS")
print("=" * 80)

print("\n1. CuAAC SMARTS Pattern:")
print(TRIAZOLE_1_4_CuAAC)
print("\nProduct template: [C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1")
print("Key observation: Uses AROMATIC nitrogen atoms [n:1], [n+:2], [n-:3]")

print("\n2. RuAAC SMARTS Pattern:")
print(TRIAZOLE_1_5_RuAAC)
print("\nProduct template: [C:5]1=[C:4]-[n:3]-[n+:2]=[n-:1]-1")
print("Key observation: Uses AROMATIC nitrogen atoms [n:3], [n+:2], [n-:1]")

print("\n" + "=" * 80)
print("GENERATED PRODUCT SMILES ANALYSIS")
print("=" * 80)

# Sample products from 01_all_generated_compounds.csv
products = {
    "SCAF-001_BB-001-Ph_1 (1,4-CuAAC)": "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OCN2[N+]=[N-]C=C2C2=CC=CC=C2)[C@@H]1OC(=O)C1=CC=CC=C1[N+](=O)[O-]",
    "SCAF-001_BB-001-Ph_2 (1,5-RuAAC)": "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OCN2C=C(C3=CC=CC=C3)[N-]=[N+]2)[C@@H]1OC(=O)C1=CC=CC=C1[N+](=O)[O-]",
}

for name, smiles in products.items():
    print("\n" + name)
    print("SMILES:", smiles)
    
    # Extract triazole portion
    if "1,4-CuAAC" in name:
        triazole_part = smiles.split("OCN2")[1].split("C2=CC")[0]
        print("Triazole portion:", "OCN2" + triazole_part + "C2")
        print("Pattern: [N+]=[N-]C=C - NON-AROMATIC charged triazole")
    else:
        triazole_part = smiles.split("OCN2")[1].split(")[C@@H]")[0]
        print("Triazole portion:", "OCN2" + triazole_part + ")")
        print("Pattern: C=C...[N-]=[N+] - NON-AROMATIC charged triazole")

print("\n" + "=" * 80)
print("ROOT CAUSE IDENTIFICATION")
print("=" * 80)

print("\nDISCREPANCY:")
print("- SMARTS patterns specify AROMATIC nitrogens: [n:1], [n+:2], [n-:3]")
print("- Generated products contain NON-AROMATIC charged nitrogens: [N+]=[N-]")
print("\nThis indicates that RDKit's reaction engine is:")
print("1. Correctly applying the SMARTS to create aromatic triazole products")
print("2. BUT the subsequent SMILES export (MolToSmiles) is KEKULIZING the aromatic ring")
print("3. Kekulization converts aromatic [n] to alternating single/double bonds with charges")
print("4. Result: [n+]=[n-] becomes [N+]=[N-] in the exported SMILES")

print("\n" + "=" * 80)
print("WHY THIS CAUSES PHASE 3 FAILURE")
print("=" * 80)

print("\nWhen Phase 3 tries to parse the kekulized SMILES:")
print("1. MolFromSmiles sees [N+]=[N-] pattern (non-aromatic)")
print("2. RDKit attempts to re-aromatize during sanitization")
print("3. The valence calculation conflicts with the explicit charges")
print("4. Result: 'Explicit valence for atom #14 N, 3, is greater than permitted'")

print("\n" + "=" * 80)
print("PROPOSED PHASE 2 FIX")
print("=" * 80)

print("\nOption 1: Disable kekulization in MolToSmiles")
print("- Change: Chem.MolToSmiles(mol, isomericSmiles=True)")
print("- To: Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=False)")
print("- This preserves aromatic notation in exported SMILES")

print("\nOption 2: Use canonical aromatic SMILES export")
print("- Add: doKekule=False parameter to MolToSmiles")
print("- This ensures aromatic rings remain aromatic in SMILES")

print("\nOption 3: Post-process to re-aromatize before export")
print("- After reaction, call Chem.Kekulize(mol, clearAromaticFlags=False)")
print("- Then export with aromatic flags preserved")

print("\nRECOMMENDED: Option 1 (disable kekulization)")
print("- Minimal change")
print("- Preserves aromatic chemistry")
print("- Generates RDKit-sanitizable SMILES")
print("- Maintains regiochemistry information")

print("\n" + "=" * 80)
