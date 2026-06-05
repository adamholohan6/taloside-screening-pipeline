from rdkit import Chem

# Pyranose ring atoms only (C1, C2, C3, C4, C5, O5)
ring_smarts = "C1(O)C(O)C(O)C(CO)OC1"
pat = Chem.MolFromSmarts(ring_smarts)

taloside_smiles = "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)C1=CC=CC=C1[N+](=O)[O-]"
mol = Chem.MolFromSmiles(taloside_smiles)
match = mol.GetSubstructMatch(pat)
if match:
    print(f"Matches for pyranose ring: {match}")
    for idx in match:
        print(f"Atom {idx}: {mol.GetAtomWithIdx(idx).GetSymbol()}")
