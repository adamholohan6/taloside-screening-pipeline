from rdkit import Chem

smiles = r"PASTE_THE_FAILING_SMILES_HERE"

print("SMILES:")
print(smiles)
print()

m = Chem.MolFromSmiles(smiles, sanitize=False)

print("Loaded:", m is not None)

if m is None:
    raise SystemExit("Could not load molecule")

print("\nTrying sanitization...")

try:
    Chem.SanitizeMol(m)
    print("Sanitization succeeded")
except Exception as e:
    print("Sanitization FAILED")
    print(type(e).__name__)
    print(e)

print("\nNitrogen atoms:")

for atom in m.GetAtoms():
    if atom.GetAtomicNum() == 7:
        print(
            f"idx={atom.GetIdx()} "
            f"charge={atom.GetFormalCharge()} "
            f"degree={atom.GetDegree()} "
            f"aromatic={atom.GetIsAromatic()}"
        )