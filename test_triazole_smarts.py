"""
Standalone diagnostic test for triazole SMARTS patterns.

Tests current vs proposed neutral-nitrogen SMARTS patterns to determine
whether proposed patterns produce RDKit-valid triazoles.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Test molecules
AZIDE_SMILES = "N=[N+]=[N-]CO"
PHENYL_ACETYLENE_SMILES = "C#Cc1ccccc1"

# Current SMARTS patterns (with charged aromatic nitrogens)
CURRENT_CUAAC_SMARTS = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"
CURRENT_RUAAC_SMARTS = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n:3]-[n+:2]=[n-:1]-1"

# Proposed SMARTS patterns (with neutral nitrogens)
PROPOSED_CUAAC_SMARTS = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5][N:1][N:2]=[N:3]1"
PROPOSED_RUAAC_SMARTS = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4][N:1][N:2]=[N:3]1"


def test_smarts_pattern(smarts_pattern, pattern_name):
    """
    Test a single SMARTS pattern with azide + phenyl acetylene.
    """
    print("=" * 80)
    print("TESTING:", pattern_name)
    print("SMARTS:", smarts_pattern)
    print("=" * 80)
    
    # Build reaction
    try:
        rxn = Chem.AllChem.ReactionFromSmarts(smarts_pattern)
        print("Reaction built: SUCCESS")
    except Exception as e:
        print("Reaction built: FAILED -", str(e))
        return False
    
    # Parse reactants (use sanitize=False for charged molecules)
    azide = Chem.MolFromSmiles(AZIDE_SMILES, sanitize=False)
    phenyl_acetylene = Chem.MolFromSmiles(PHENYL_ACETYLENE_SMILES, sanitize=False)
    
    if azide is None:
        print("Azide parsing: FAILED")
        return False
    if phenyl_acetylene is None:
        print("Phenyl acetylene parsing: FAILED")
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
    first_product_tuple = products[0]
    if not first_product_tuple:
        print("First product tuple is empty")
        return False
    
    product = first_product_tuple[0]
    print("Product obtained: SUCCESS")
    print("Total atoms:", product.GetNumAtoms())
    
    # Print raw SMILES
    print("\n=== RAW PRODUCT ===")
    try:
        raw_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
        print("MolToSmiles(product):", raw_smiles)
    except Exception as e:
        print("MolToSmiles(product): FAILED -", str(e))
        raw_smiles = None
    
    # Print atom table
    print("\n=== ATOM TABLE ===")
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
    print("\n=== SANITIZATION TEST ===")
    try:
        Chem.SanitizeMol(product, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        print("SanitizeMol (no kekulize): SUCCESS")
    except Exception as e:
        print("SanitizeMol (no kekulize): FAILED -", str(e))
        return False
    
    # Test SMILES export with different settings
    print("\n=== SMILES EXPORT TESTS ===")
    
    smiles_default = None
    smiles_kekule_true = None
    smiles_kekule_false = None
    
    try:
        smiles_default = Chem.MolToSmiles(product, isomericSmiles=True)
        print("MolToSmiles(default):", smiles_default)
    except Exception as e:
        print("MolToSmiles(default): FAILED -", str(e))
    
    try:
        smiles_kekule_true = Chem.MolToSmiles(product, isomericSmiles=True, kekuleSmiles=True)
        print("MolToSmiles(kekuleSmiles=True):", smiles_kekule_true)
    except Exception as e:
        print("MolToSmiles(kekuleSmiles=True): FAILED -", str(e))
    
    try:
        smiles_kekule_false = Chem.MolToSmiles(product, isomericSmiles=True, kekuleSmiles=False)
        print("MolToSmiles(kekuleSmiles=False):", smiles_kekule_false)
    except Exception as e:
        print("MolToSmiles(kekuleSmiles=False): FAILED -", str(e))
    
    # Export canonical smiles
    print("\n=== CANONICAL SMILES EXPORT ===")
    if smiles_default:
        # Re-parse the SMILES to get a molecule, then canonicalize
        canonical_mol = Chem.MolFromSmiles(smiles_default, sanitize=False)
        if canonical_mol is not None:
            canonical_smiles = Chem.MolToSmiles(canonical_mol, isomericSmiles=True, canonical=True)
            print("Canonical SMILES:", canonical_smiles)
        else:
            print("Cannot re-parse SMILES for canonicalization")
            return False
    else:
        print("Cannot export canonical SMILES (no default SMILES available)")
        return False
    
    # Attempt re-import
    print("\n=== RE-IMPORT TEST ===")
    try:
        reimported_mol = Chem.MolFromSmiles(canonical_smiles)
        if reimported_mol is None:
            print("MolFromSmiles(canonical_smiles): FAILED - returned None")
            return False
        else:
            print("MolFromSmiles(canonical_smiles): SUCCESS")
            print("Re-imported atoms:", reimported_mol.GetNumAtoms())
    except Exception as e:
        print("MolFromSmiles(canonical_smiles): FAILED -", str(e))
        return False
    
    print("\n=== RESULT ===")
    print("PASS")
    print("=" * 80)
    return True


def main():
    print("\n")
    print("*" * 80)
    print("TRIAZOLE SMARTS DIAGNOSTIC TEST")
    print("*" * 80)
    print("\n")
    
    results = {}
    
    # Test current SMARTS patterns
    print("\n\n")
    print("#" * 80)
    print("# CURRENT SMARTS PATTERNS")
    print("#" * 80)
    
    results["current_cuaac"] = test_smarts_pattern(CURRENT_CUAAC_SMARTS, "Current CuAAC")
    results["current_ruaac"] = test_smarts_pattern(CURRENT_RUAAC_SMARTS, "Current RuAAC")
    
    # Test proposed SMARTS patterns
    print("\n\n")
    print("#" * 80)
    print("# PROPOSED SMARTS PATTERNS")
    print("#" * 80)
    
    results["proposed_cuaac"] = test_smarts_pattern(PROPOSED_CUAAC_SMARTS, "Proposed CuAAC")
    results["proposed_ruaac"] = test_smarts_pattern(PROPOSED_RUAAC_SMARTS, "Proposed RuAAC")
    
    # Summary
    print("\n\n")
    print("*" * 80)
    print("SUMMARY")
    print("*" * 80)
    print("\n")
    print("Pattern                    | Result")
    print("-" * 80)
    print("Current CuAAC              |", "PASS" if results["current_cuaac"] else "FAIL")
    print("Current RuAAC              |", "PASS" if results["current_ruaac"] else "FAIL")
    print("Proposed CuAAC             |", "PASS" if results["proposed_cuaac"] else "FAIL")
    print("Proposed RuAAC             |", "PASS" if results["proposed_ruaac"] else "FAIL")
    print("\n")


if __name__ == "__main__":
    main()
