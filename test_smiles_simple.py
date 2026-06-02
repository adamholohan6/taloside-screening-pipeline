"""
Simple test to parse a single SMILES string with RDKit
"""

from rdkit import Chem

# Test SMILES from Phase 2 output
smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OCN2[N+]=[N-]C=C2C2=CC=CC=C2)[C@@H]1OC(=O)C1=CC=CC=C1[N+](=O)[O-]"

print("=" * 70)
print("TESTING SMILES PARSING")
print("=" * 70)
print("SMILES:", smiles)
print("=" * 70)

# Step 1: Try standard MolFromSmiles
print("\nStep 1: MolFromSmiles with default sanitization")
try:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("  FAILED - returned None")
    else:
        print("  SUCCESS -", mol.GetNumAtoms(), "atoms")
except Exception as e:
    print("  EXCEPTION:", type(e).__name__, ":", e)

# Step 2: Try MolFromSmiles with sanitize=False
print("\nStep 2: MolFromSmiles with sanitize=False")
try:
    mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol_raw is None:
        print("  FAILED - returned None")
    else:
        print("  SUCCESS -", mol_raw.GetNumAtoms(), "atoms")
        print("\n  Atom table (before sanitization):")
        for atom in mol_raw.GetAtoms():
            print("    idx=", atom.GetIdx(), " sym=", atom.GetSymbol(), " charge=", atom.GetFormalCharge(), " degree=", atom.GetDegree(), " explicitHs=", atom.GetNumExplicitHs(), " aromatic=", atom.GetIsAromatic())
except Exception as e:
    print("  EXCEPTION:", type(e).__name__, ":", e)

# Step 3: Try sanitization with different flags
print("\nStep 3: Sanitization with SANITIZE_ALL ^ SANITIZE_PROPERTIES")
try:
    mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol_raw is not None:
        mask = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        Chem.SanitizeMol(mol_raw, sanitizeOps=mask)
        print("  SUCCESS")
        print("\n  Atom table (after partial sanitization):")
        for atom in mol_raw.GetAtoms():
            print("    idx=", atom.GetIdx(), " sym=", atom.GetSymbol(), " charge=", atom.GetFormalCharge(), " degree=", atom.GetDegree(), " explicitHs=", atom.GetNumExplicitHs(), " aromatic=", atom.GetIsAromatic())
except Exception as e:
    print("  EXCEPTION:", type(e).__name__, ":", e)
    import traceback
    traceback.print_exc()

# Step 4: Try UpdatePropertyCache
print("\nStep 4: UpdatePropertyCache(strict=True)")
try:
    mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol_raw is not None:
        mask = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        Chem.SanitizeMol(mol_raw, sanitizeOps=mask)
        mol_raw.UpdatePropertyCache(strict=True)
        print("  SUCCESS")
except Exception as e:
    print("  EXCEPTION:", type(e).__name__, ":", e)
    import traceback
    traceback.print_exc()

# Step 5: Try full sanitization
print("\nStep 5: Full SanitizeMol")
try:
    mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol_raw is not None:
        mask = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
        Chem.SanitizeMol(mol_raw, sanitizeOps=mask)
        mol_raw.UpdatePropertyCache(strict=True)
        Chem.SanitizeMol(mol_raw)
        print("  SUCCESS")
        print("\n  Atom table (after full sanitization):")
        for atom in mol_raw.GetAtoms():
            print("    idx=", atom.GetIdx(), " sym=", atom.GetSymbol(), " charge=", atom.GetFormalCharge(), " degree=", atom.GetDegree(), " explicitHs=", atom.GetNumExplicitHs(), " aromatic=", atom.GetIsAromatic())
except Exception as e:
    print("  EXCEPTION:", type(e).__name__, ":", e)
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
