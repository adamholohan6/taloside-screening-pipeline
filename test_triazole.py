from rdkit import Chem

smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OCN2[N+]=[N-]C=C2C2=CC=CC=C2)[C@@H]1OC(=O)C1=CC=CC=C1[N+](=O)[O-]"

print("Testing:")
print(smiles)
print()

try:
    mol = Chem.MolFromSmiles(smiles)
    print("MolFromSmiles result:", mol)

except Exception as e:
    print("Exception:")
    print(type(e))
    print(e)

print()
print("Trying sanitize=False")

mol = Chem.MolFromSmiles(smiles, sanitize=False)

print("mol =", mol)

for atom in mol.GetAtoms():
    print(
        atom.GetIdx(),
        atom.GetSymbol(),
        "charge=", atom.GetFormalCharge(),
        "degree=", atom.GetDegree(),
        "aromatic=", atom.GetIsAromatic()
    )

try:
    Chem.SanitizeMol(mol)
    print("Sanitize success")
except Exception as e:
    print("Sanitize failed:")
    print(type(e))
    print(e)