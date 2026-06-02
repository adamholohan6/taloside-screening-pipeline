"""
Chemistry Investigation: Azide-Alkyne Cycloaddition SMARTS Analysis

Investigates:
1. Canonical RDKit representation of neutral 1,2,3-triazole
2. Whether azide formal charges are preserved in triazole products
3. Atom mappings in current SMARTS that cause charge preservation
4. Comparison of valid triazole vs SMARTS-generated products
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def analyze_neutral_triazole():
    """
    Manually construct a neutral 1,2,3-triazole and analyze its RDKit representation.
    """
    print("=" * 80)
    print("1. CANONICAL RDKIT REPRESENTATION OF NEUTRAL 1,2,3-TRIAZOLE")
    print("=" * 80)
    
    # Construct neutral 1,2,3-triazole manually
    # Using SMILES: c1nnc[nH]1 (1H-1,2,3-triazole)
    neutral_triazole_smiles = "c1nnc[nH]1"
    
    print("\nConstructing neutral 1,2,3-triazole from SMILES:", neutral_triazole_smiles)
    
    mol = Chem.MolFromSmiles(neutral_triazole_smiles)
    if mol is None:
        print("FAILED to parse neutral triazole SMILES")
        return None
    
    print("Parsing: SUCCESS")
    print("Total atoms:", mol.GetNumAtoms())
    
    # Sanitize
    try:
        Chem.SanitizeMol(mol)
        print("Sanitization: SUCCESS")
    except Exception as e:
        print("Sanitization: FAILED -", str(e))
        return None
    
    # Get canonical SMILES
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    print("\nCanonical SMILES:", canonical_smiles)
    
    # Print atom table
    print("\nAtom table:")
    print("Index | Symbol | Charge | Degree | Aromatic | Hybridization")
    print("-" * 70)
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        print("{:5} | {:6} | {:6} | {:6} | {:8} | {}".format(
            i,
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetDegree(),
            atom.GetIsAromatic(),
            atom.GetHybridization()
        ))
    
    # Print bond table
    print("\nBond table:")
    print("Atom1 | Atom2 | Bond Type | Aromatic")
    print("-" * 50)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        aromatic = bond.GetIsAromatic()
        print("{:5} | {:5} | {:9} | {}".format(a1, a2, str(bond_type), aromatic))
    
    return mol, canonical_smiles


def analyze_azide():
    """
    Analyze the azide reactant structure.
    """
    print("\n" + "=" * 80)
    print("2. AZIDE REACTANT ANALYSIS")
    print("=" * 80)
    
    # Azide SMILES with charges
    azide_smiles = "N=[N+]=[N-]"
    
    print("\nAzide SMILES:", azide_smiles)
    
    mol = Chem.MolFromSmiles(azide_smiles, sanitize=False)
    if mol is None:
        print("FAILED to parse azide SMILES")
        return None
    
    print("Parsing (sanitize=False): SUCCESS")
    print("Total atoms:", mol.GetNumAtoms())
    
    # Print atom table
    print("\nAtom table:")
    print("Index | Symbol | Charge | Degree | Aromatic")
    print("-" * 50)
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        print("{:5} | {:6} | {:6} | {:6} | {}".format(
            i,
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetDegree(),
            atom.GetIsAromatic()
        ))
    
    return mol


def analyze_current_smarts():
    """
    Analyze the current SMARTS patterns and their atom mappings.
    """
    print("\n" + "=" * 80)
    print("3. CURRENT SMARTS PATTERN ANALYSIS")
    print("=" * 80)
    
    current_cuaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"
    current_ruaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n:3]-[n+:2]=[n-:1]-1"
    
    print("\nCurrent CuAAC SMARTS:")
    print(current_cuaac)
    print("\nAtom mapping analysis:")
    print("Reactant side:")
    print("  [N:1] - neutral nitrogen (azide terminal)")
    print("  [N+:2] - positively charged nitrogen (azide central)")
    print("  [N-:3] - negatively charged nitrogen (azide terminal)")
    print("  [C:4] - alkyne carbon 1")
    print("  [C:5] - alkyne carbon 2")
    print("\nProduct side (CuAAC):")
    print("  [C:4]1 - maps to alkyne C4, becomes ring carbon")
    print("  =[C:5] - maps to alkyne C5, becomes ring carbon with double bond")
    print("  -[n:1] - maps to azide N1, becomes aromatic nitrogen")
    print("  -[n+:2] - maps to azide N2, becomes aromatic nitrogen WITH CHARGE PRESERVED")
    print("  =[n-:3] - maps to azide N3, becomes aromatic nitrogen WITH CHARGE PRESERVED")
    print("  -1 - ring closure")
    
    print("\n" + "-" * 80)
    print("\nCurrent RuAAC SMARTS:")
    print(current_ruaac)
    print("\nAtom mapping analysis:")
    print("Reactant side: (same as CuAAC)")
    print("\nProduct side (RuAAC):")
    print("  [C:5]1 - maps to alkyne C5, becomes ring carbon")
    print("  =[C:4] - maps to alkyne C4, becomes ring carbon with double bond")
    print("  -[n:3] - maps to azide N3, becomes aromatic nitrogen")
    print("  -[n+:2] - maps to azide N2, becomes aromatic nitrogen WITH CHARGE PRESERVED")
    print("  =[n-:1] - maps to azide N1, becomes aromatic nitrogen WITH CHARGE PRESERVED")
    print("  -1 - ring closure")
    
    print("\n" + "=" * 80)
    print("KEY FINDING:")
    print("=" * 80)
    print("The SMARTS patterns explicitly map charged azide nitrogens [N+:2] and [N-:3]")
    print("to charged aromatic triazole nitrogens [n+:2] and [n-:3] in the product.")
    print("This is why the formal charges from the azide are preserved in the triazole.")
    print("=" * 80)


def test_neutral_smarts():
    """
    Test what happens if we use neutral nitrogens in the product template.
    """
    print("\n" + "=" * 80)
    print("4. TESTING NEUTRAL NITROGEN PRODUCT TEMPLATE")
    print("=" * 80)
    
    # Proposed SMARTS with neutral nitrogens in product
    neutral_cuaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n:2]=[n:3]-1"
    
    print("\nProposed CuAAC SMARTS (neutral product nitrogens):")
    print(neutral_cuaac)
    print("\nAtom mapping analysis:")
    print("Reactant side:")
    print("  [N:1] - neutral nitrogen (azide terminal)")
    print("  [N+:2] - positively charged nitrogen (azide central)")
    print("  [N-:3] - negatively charged nitrogen (azide terminal)")
    print("\nProduct side:")
    print("  [n:1] - maps to azide N1, becomes aromatic nitrogen (charge NOT preserved)")
    print("  [n:2] - maps to azide N2, becomes aromatic nitrogen (charge NOT preserved)")
    print("  [n:3] - maps to azide N3, becomes aromatic nitrogen (charge NOT preserved)")
    
    print("\nTesting this SMARTS with reactants...")
    
    # Build reaction
    rxn = AllChem.ReactionFromSmarts(neutral_cuaac)
    print("Reaction built: SUCCESS")
    
    # Parse reactants
    azide = Chem.MolFromSmiles("N=[N+]=[N-]CO", sanitize=False)
    phenyl_acetylene = Chem.MolFromSmiles("C#Cc1ccccc1", sanitize=False)
    
    if azide is None or phenyl_acetylene is None:
        print("Reactant parsing: FAILED")
        return
    
    print("Reactants parsed: SUCCESS")
    
    # Run reaction
    products = rxn.RunReactants((azide, phenyl_acetylene))
    print("RunReactants(): SUCCESS")
    print("Number of product tuples:", len(products))
    
    if not products:
        print("No products generated")
        return
    
    # Get first product
    product = products[0][0]
    print("Product obtained: SUCCESS")
    
    # Print atom table
    print("\nProduct atom table:")
    print("Index | Symbol | Charge | Degree | Aromatic")
    print("-" * 50)
    for i in range(product.GetNumAtoms()):
        atom = product.GetAtomWithIdx(i)
        print("{:5} | {:6} | {:6} | {:6} | {}".format(
            i,
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetDegree(),
            atom.GetIsAromatic()
        ))
    
    # Attempt sanitization
    print("\nAttempting sanitization...")
    try:
        Chem.SanitizeMol(product, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        print("Sanitization: SUCCESS")
    except Exception as e:
        print("Sanitization: FAILED -", str(e))


def main():
    print("\n")
    print("*" * 80)
    print("CHEMISTRY INVESTIGATION: AZIDE-ALKYNE CYCLOADDITION SMARTS")
    print("*" * 80)
    print("\n")
    
    # 1. Analyze neutral triazole
    neutral_triazole, canonical_smiles = analyze_neutral_triazole()
    
    # 2. Analyze azide
    azide = analyze_azide()
    
    # 3. Analyze current SMARTS
    analyze_current_smarts()
    
    # 4. Test neutral SMARTS
    test_neutral_smarts()
    
    # 5. Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("\n1. Canonical neutral 1,2,3-triazole SMILES:", canonical_smiles if canonical_smiles else "N/A")
    print("\n2. Azide structure: N=[N+]=[N-] (terminal N neutral, central N+, terminal N-)")
    print("\n3. Current SMARTS issue:")
    print("   - Reactant: [N+:2] and [N-:3] have formal charges")
    print("   - Product: [n+:2] and [n-:3] preserve these charges")
    print("   - Result: Charged aromatic triazole nitrogens (chemically invalid)")
    print("\n4. Root cause:")
    print("   The SMARTS atom mappings explicitly preserve azide charges in the product.")
    print("   In a real 1,3-dipolar cycloaddition, charges are redistributed during")
    print("   ring formation, resulting in a neutral aromatic triazole.")
    print("\n5. Solution:")
    print("   Product template should use neutral nitrogens [n:1][n:2][n:3]")
    print("   to allow RDKit to assign correct charges during sanitization.")
    print("=" * 80)


if __name__ == "__main__":
    main()
