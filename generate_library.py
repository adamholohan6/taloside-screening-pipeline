import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# 1. Paste your exact extracted SMILES strings from ChemDraw here
real_thesis_compounds = {
    "Talo-1": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC(OC)=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-2": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=C(OC)C=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-3": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-4": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC=C(OC)C=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-5": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CN=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-6": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC(F)=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-7": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC(Cl)=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-8": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC(I)=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-9": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC(Br)=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4",
    "Talo-10": "O=C(O[C@H]1[C@@H](OCN2C=C(C3=CC(N)=CC=C3)N=N2)[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4"
}

# 2. Create an empty list to store our calculated data
results = []

# 3. Loop through each molecule to calculate its descriptors
for name, smiles in real_thesis_compounds.items():
    
    # Convert the SMILES text string into a 3D RDKit Molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    # If the SMILES string is valid, calculate the properties
    if mol is not None:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Save the results for this molecule
        results.append({
            "Compound": name,
            "SMILES": smiles,
            "Molecular_Weight": mw,
            "LogP": logp,
            "H_Donors": h_donors,
            "H_Acceptors": h_acceptors,
            "TPSA": tpsa,
            "Rotatable_Bonds": rot_bonds
        })
    else:
        print(f"Error: Could not process {name}. Check the SMILES string.")

# 4. Convert the data into a Pandas DataFrame (a virtual spreadsheet)
df = pd.DataFrame(results)

# 5. Save the spreadsheet to your computer
output_filename = "taloside_descriptors.csv"
df.to_csv(output_filename, index=False)

print(f"Success! Calculated descriptors for {len(df)} compounds.")
print(f"Data saved to {output_filename}")