"""
Diagnostic script to gather evidence for triazole valence error
"""

from rdkit import Chem

# The exact failing product_smiles for the first compound from Phase 2
smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OCN2[N+]=[N-]C=C2C2=CC=CC=C2)[C@@H]1OC(=O)C1=CC=CC=C1[N+](=O)[O-]"

print("=" * 80)
print("DIAGNOSTIC EVIDENCE FOR TRIAZOLE VALENCE ERROR")
print("=" * 80)
print("\n1. EXACT FAILING PRODUCT_SMILES:")
print(smiles)
print("\nCompound ID: SCAF-001_BB-001-Ph_1")
print("Regioisomer: 1,4-CuAAC")
print("=" * 80)

# Load with sanitize=False to see the raw atom table
print("\n2. ATOM TABLE (Chem.MolFromSmiles(smiles, sanitize=False)):")
print("=" * 80)

mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)

if mol_raw is None:
    print("ERROR: MolFromSmiles(sanitize=False) returned None")
else:
    print("Total atoms:", mol_raw.GetNumAtoms())
    print("\n{:<6} {:<6} {:<12} {:<8} {:<10} {:<14} {:<14}".format(
        "Index", "Elem", "FormalChg", "Degree", "Aromatic", "ExplicitValence", "ExplicitHs"
    ))
    print("-" * 80)
    
    for atom in mol_raw.GetAtoms():
        idx = atom.GetIdx()
        elem = atom.GetSymbol()
        fc = atom.GetFormalCharge()
        degree = atom.GetDegree()
        aromatic = atom.GetIsAromatic()
        explicit_valence = atom.GetExplicitValence()
        explicit_hs = atom.GetNumExplicitHs()
        
        print("{:<6} {:<6} {:+<12} {:<8} {:<10} {:<14} {:<14}".format(
            idx, elem, fc, degree, aromatic, explicit_valence, explicit_hs
        ))

print("\n" + "=" * 80)
print("3. ATOM #14 DETAILS:")
print("=" * 80)

if mol_raw and mol_raw.GetNumAtoms() > 14:
    atom_14 = mol_raw.GetAtomWithIdx(14)
    print("Index:", atom_14.GetIdx())
    print("Element:", atom_14.GetSymbol())
    print("Formal Charge:", atom_14.GetFormalCharge())
    print("Degree (number of bonded neighbors):", atom_14.GetDegree())
    print("Aromatic:", atom_14.GetIsAromatic())
    print("Explicit Valence:", atom_14.GetExplicitValence())
    print("Explicit Hs:", atom_14.GetNumExplicitHs())
    print("Total Valence:", atom_14.GetTotalValence())
    print("Implicit Valence:", atom_14.GetImplicitValence())
    
    # Get bonded atoms
    print("\nBonded neighbors:")
    for neighbor in atom_14.GetNeighbors():
        print("  - Atom", neighbor.GetIdx(), ":", neighbor.GetSymbol(), "(charge=", neighbor.GetFormalCharge(), ")")
    
    # Get bonds
    print("\nBonds:")
    for bond in atom_14.GetBonds():
        other_atom = bond.GetOtherAtom(atom_14)
        bond_type = bond.GetBondType()
        print("  - To atom", other_atom.GetIdx(), "(", other_atom.GetSymbol(), "):", bond_type)

print("\n" + "=" * 80)
print("4. TRIAZOLE RING ATOMS:")
print("=" * 80)

# Find the triazole ring atoms (pattern: N2[N+]=[N-]C=C2)
# The SMILES contains: OCN2[N+]=[N-]C=C2
# So the triazole ring atoms are the ones in the N2[N+]=[N-]C=C2 pattern

print("Searching for triazole ring pattern in SMILES...")
print("Pattern: N2[N+]=[N-]C=C2")
print("\nTriazole ring atoms (indices based on SMILES parsing):")

# The triazole ring in this SMILES is: N2[N+]=[N-]C=C2
# This corresponds to atoms around the N2...C=C2 closure
# Let's identify them by their properties

triazole_atoms = []
for atom in mol_raw.GetAtoms():
    if atom.GetSymbol() == "N" and atom.GetIsAromatic():
        triazole_atoms.append(atom.GetIdx())

print("Aromatic nitrogen atoms (potential triazole):", triazole_atoms)

# Also look for the aromatic carbons in the triazole
triazole_carbons = []
for atom in mol_raw.GetAtoms():
    if atom.GetSymbol() == "C" and atom.GetIsAromatic():
        # Check if it's in a small ring (triazole is 5-membered)
        for ring in Chem.GetSymmSSSR(mol_raw):
            if len(ring) == 5 and atom.GetIdx() in ring:
                triazole_carbons.append(atom.GetIdx())
                break

print("Aromatic carbon atoms in 5-membered rings:", triazole_carbons)

# Get the 5-membered ring containing the triazole
print("\n5-membered rings:")
for ring in Chem.GetSymmSSSR(mol_raw):
    if len(ring) == 5:
        print("  Ring atoms:", ring)
        # Check if this is the triazole (should have 3 N and 2 C)
        n_count = sum(1 for idx in ring if mol_raw.GetAtomWithIdx(idx).GetSymbol() == "N")
        c_count = sum(1 for idx in ring if mol_raw.GetAtomWithIdx(idx).GetSymbol() == "C")
        print("    Composition:", n_count, "N,", c_count, "C")
        if n_count == 3 and c_count == 2:
            print("    -> THIS IS THE TRIAZOLE RING")
            print("\n    Detailed atom table for triazole ring:")
            for idx in ring:
                atom = mol_raw.GetAtomWithIdx(idx)
                print("      Atom", idx, ":", atom.GetSymbol(), ", charge=", atom.GetFormalCharge(), ", ",
                      "degree=", atom.GetDegree(), ", aromatic=", atom.GetIsAromatic(), ", ",
                      "explicit_valence=", atom.GetExplicitValence(), ", explicit_Hs=", atom.GetNumExplicitHs())

print("\n" + "=" * 80)
print("5. WHY RDKIT THINKS ATOM #14 IS OVER-VALENT:")
print("=" * 80)

if mol_raw and mol_raw.GetNumAtoms() > 14:
    atom_14 = mol_raw.GetAtomWithIdx(14)
    
    print("Atom #14 is", atom_14.GetSymbol(), "with:")
    print("  - Formal charge:", atom_14.GetFormalCharge())
    print("  - Degree (bonds):", atom_14.GetDegree())
    print("  - Explicit valence:", atom_14.GetExplicitValence())
    print("  - Explicit Hs:", atom_14.GetNumExplicitHs())
    
    print("\nRDKit's valence rules for nitrogen:")
    print("  - Neutral N (charge 0): max valence = 3 (3 bonds + Hs)")
    print("  - Positively charged N (charge +1): max valence = 4")
    print("  - Negatively charged N (charge -1): max valence = 2")
    
    print("\nAtom #14 analysis:")
    if atom_14.GetFormalCharge() == 0:
        max_allowed = 3
        print("  - Neutral nitrogen: max allowed valence =", max_allowed)
        print("  - Current explicit valence =", atom_14.GetExplicitValence())
        if atom_14.GetExplicitValence() > max_allowed:
            print("  - ERROR: Valence", atom_14.GetExplicitValence(), "exceeds max", max_allowed)
            print("  - This is the 'Explicit valence for atom #14 N, 3, is greater than permitted' error")
    elif atom_14.GetFormalCharge() == 1:
        max_allowed = 4
        print("  - Positively charged nitrogen: max allowed valence =", max_allowed)
        print("  - Current explicit valence =", atom_14.GetExplicitValence())
    elif atom_14.GetFormalCharge() == -1:
        max_allowed = 2
        print("  - Negatively charged nitrogen: max allowed valence =", max_allowed)
        print("  - Current explicit valence =", atom_14.GetExplicitValence())
        if atom_14.GetExplicitValence() > max_allowed:
            print("  - ERROR: Valence", atom_14.GetExplicitValence(), "exceeds max", max_allowed)

print("\n" + "=" * 80)
