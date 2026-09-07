from rdkit import Chem

# Pyranose core of galactoside/taloside
pyranose_smarts = "[C@H]1([O,N,C])[C@@H]([O,N,C])[C@@H]([O,N,C])[C@@H](CO)O[C@H]1(OC)"
pat = Chem.MolFromSmarts(pyranose_smarts)

taloside_smiles = "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)C1=CC=CC=C1[N+](=O)[O-]"
mol = Chem.MolFromSmiles(taloside_smiles)
match = mol.GetSubstructMatch(pat)
print(f"Matches for pyranose core: {match}")
