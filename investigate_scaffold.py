"""
Investigate scaffold molecule atom mapping.
"""

from rdkit import Chem

# Scaffold from phase2_integration.py
SCAFFOLD_SMILES = (
    "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)"
    "C4=C([N+]([O-])=O)C=CC=C4"
)

print("Scaffold SMILES:", SCAFFOLD_SMILES)
print()

# Parse scaffold
scaffold = Chem.MolFromSmiles(SCAFFOLD_SMILES, sanitize=False)
if scaffold is None:
    print("Failed to parse scaffold")
    exit(1)

print("Total atoms:", scaffold.GetNumAtoms())
print()

# Print atom table with map numbers
print("ATOM TABLE:")
print("Index | Symbol | MapNum | Charge | Degree")
print("-" * 50)
for i in range(scaffold.GetNumAtoms()):
    atom = scaffold.GetAtomWithIdx(i)
    map_num = atom.GetAtomMapNum()
    print("{:5} | {:6} | {:6} | {:6} | {:6}".format(
        i,
        atom.GetSymbol(),
        map_num,
        atom.GetFormalCharge(),
        atom.GetDegree()
    ))

# Count mapped atoms
mapped_atoms = [i for i in range(scaffold.GetNumAtoms()) if scaffold.GetAtomWithIdx(i).GetAtomMapNum() > 0]
print("\nMapped atoms in scaffold:", len(mapped_atoms))
print("Mapped atom indices:", mapped_atoms)

# Find azide nitrogens
print("\nAzide nitrogens:")
for i in range(scaffold.GetNumAtoms()):
    atom = scaffold.GetAtomWithIdx(i)
    if atom.GetSymbol() == "N":
        neighbor_atoms = atom.GetNeighbors()
        neighbor_symbols = [n.GetSymbol() for n in neighbor_atoms]
        print("  N{}: charge={}, degree={}, neighbors={}".format(
            i, atom.GetFormalCharge(), atom.GetDegree(), neighbor_symbols
        ))
