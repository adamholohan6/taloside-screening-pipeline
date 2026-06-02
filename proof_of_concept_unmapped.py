"""
Proof-of-Concept: Unmapped Product Nitrogens for Azide-Alkyne Cycloaddition

Tests whether removing atom mapping numbers from product nitrogens in the
actual CuAAC and RuAAC SMARTS patterns produces neutral, RDKit-valid triazoles.

Based on finding: Product templates without atom maps produce neutral triazoles.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


def test_unmapped_smarts(smarts_pattern, pattern_name):
    """
    Test a SMARTS pattern and report results.
    """
    print("=" * 80)
    print("TESTING:", pattern_name)
    print("SMARTS:", smarts_pattern)
    print("=" * 80)
    
    # Build reaction
    try:
        rxn = AllChem.ReactionFromSmarts(smarts_pattern)
        print("Reaction built: SUCCESS")
    except Exception as e:
        print("Reaction built: FAILED -", str(e))
        return False
    
    # Parse reactants (use actual molecules from pipeline)
    # Azide: N=[N+]=[N-]CO (from pipeline)
    # Phenyl acetylene: C#Cc1ccccc1 (from pipeline)
    azide = Chem.MolFromSmiles("N=[N+]=[N-]CO", sanitize=False)
    phenyl_acetylene = Chem.MolFromSmiles("C#Cc1ccccc1", sanitize=False)
    
    if azide is None or phenyl_acetylene is None:
        print("Reactant parsing: FAILED")
        return False
    
    print("Reactants parsed: SUCCESS")
    
    # Run reaction
    try:
        products = rxn.RunReactants((azide, phenyl_acetylene))
        print("RunReactants(): SUCCESS")
        print("Number of product tuples:", len(products))
    except Exception as e:
        print("RunReactants(): FAILED -", str(e))
        return False
    
    if not products:
        print("No products generated")
        return False
    
    # Get first product
    product = products[0][0]
    print("Product obtained: SUCCESS")
    print("Total atoms:", product.GetNumAtoms())
    
    # 1. Product nitrogen charges after RunReactants
    print("\n1. PRODUCT NITROGEN CHARGES AFTER RunReactants():")
    print("Index | Symbol | Charge | Degree | Aromatic")
    print("-" * 50)
    for i in range(product.GetNumAtoms()):
        atom = product.GetAtomWithIdx(i)
        if atom.GetSymbol() == "N":
            print("{:5} | {:6} | {:6} | {:6} | {}".format(
                i,
                atom.GetSymbol(),
                atom.GetFormalCharge(),
                atom.GetDegree(),
                atom.GetIsAromatic()
            ))
    
    # 2. Sanitization success
    print("\n2. SANITIZATION TEST:")
    try:
        Chem.SanitizeMol(product, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        print("Sanitization: SUCCESS")
        sanitization_success = True
    except Exception as e:
        print("Sanitization: FAILED -", str(e))
        sanitization_success = False
    
    if not sanitization_success:
        return False
    
    # 3. MolToSmiles round-trip
    print("\n3. MOLTOSMILES ROUND-TRIP TEST:")
    try:
        smiles = Chem.MolToSmiles(product, isomericSmiles=True)
        print("MolToSmiles(product):", smiles)
    except Exception as e:
        print("MolToSmiles(product): FAILED -", str(e))
        return False
    
    try:
        reimported = Chem.MolFromSmiles(smiles)
        if reimported is None:
            print("MolFromSmiles(smiles): FAILED - returned None")
            return False
        else:
            print("MolFromSmiles(smiles): SUCCESS")
            print("Re-imported atoms:", reimported.GetNumAtoms())
            roundtrip_success = True
    except Exception as e:
        print("MolFromSmiles(smiles): FAILED -", str(e))
        roundtrip_success = False
    
    # 4. Chemical correctness of triazole
    print("\n4. CHEMICAL CORRECTNESS CHECK:")
    
    # Find triazole ring
    triazole_ring = None
    for ring in Chem.GetSymmSSSR(product):
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if product.GetAtomWithIdx(idx).GetSymbol() == "N")
            c_count = sum(1 for idx in ring if product.GetAtomWithIdx(idx).GetSymbol() == "C")
            if n_count == 3 and c_count == 2:
                triazole_ring = ring
                break
    
    if triazole_ring:
        print("Triazole ring found: atoms", list(triazole_ring))
        
        # Check if all nitrogens are neutral
        all_neutral = True
        for idx in triazole_ring:
            atom = product.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "N" and atom.GetFormalCharge() != 0:
                all_neutral = False
                print("  WARNING: Nitrogen {} has charge {}".format(idx, atom.GetFormalCharge()))
        
        if all_neutral:
            print("  All triazole nitrogens are neutral: CORRECT")
        else:
            print("  Triazole has charged nitrogens: INCORRECT")
        
        # Check if aromatic
        all_aromatic = all(product.GetAtomWithIdx(idx).GetIsAromatic() for idx in triazole_ring)
        if all_aromatic:
            print("  All triazole atoms are aromatic: CORRECT")
        else:
            print("  Triazole is not fully aromatic: INCORRECT")
        
        chemical_correctness = all_neutral and all_aromatic
    else:
        print("Triazole ring not found")
        chemical_correctness = False
    
    print("\n" + "=" * 80)
    print("RESULT:")
    print("  Nitrogen charges: Neutral" if all_neutral else "  Nitrogen charges: Charged")
    print("  Sanitization: SUCCESS" if sanitization_success else "  Sanitization: FAILED")
    print("  Round-trip: SUCCESS" if roundtrip_success else "  Round-trip: FAILED")
    print("  Chemical correctness: CORRECT" if chemical_correctness else "  Chemical correctness: INCORRECT")
    print("=" * 80)
    
    return sanitization_success and roundtrip_success and chemical_correctness


def main():
    print("\n")
    print("*" * 80)
    print("PROOF-OF-CONCEPT: UNMAPPED PRODUCT NITROGENS")
    print("*" * 80)
    print("\n")
    
    # Original SMARTS from pipeline (with mapped product nitrogens)
    original_cuaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"
    original_ruaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n:3]-[n+:2]=[n-:1]-1"
    
    # Modified SMARTS (unmapped product nitrogens)
    # Remove atom map numbers from product nitrogens only
    modified_cuaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n]-[n]=[n]-1"
    modified_ruaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n]-[n]=[n]-1"
    
    print("\n" + "#" * 80)
    print("# ORIGINAL SMARTS (with mapped product nitrogens)")
    print("#" * 80)
    
    original_cuaac_result = test_unmapped_smarts(original_cuaac, "Original CuAAC")
    original_ruaac_result = test_unmapped_smarts(original_ruaac, "Original RuAAC")
    
    print("\n" + "#" * 80)
    print("# MODIFIED SMARTS (with unmapped product nitrogens)")
    print("#" * 80)
    
    modified_cuaac_result = test_unmapped_smarts(modified_cuaac, "Modified CuAAC (unmapped)")
    modified_ruaac_result = test_unmapped_smarts(modified_ruaac, "Modified RuAAC (unmapped)")
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("\nOriginal CuAAC:", "PASS" if original_cuaac_result else "FAIL")
    print("Original RuAAC:", "PASS" if original_ruaac_result else "FAIL")
    print("Modified CuAAC (unmapped):", "PASS" if modified_cuaac_result else "FAIL")
    print("Modified RuAAC (unmapped):", "PASS" if modified_ruaac_result else "FAIL")
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
