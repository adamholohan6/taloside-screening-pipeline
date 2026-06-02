"""
Atom Mapping Investigation: Azide-Alkyne Cycloaddition

Investigates minimum atom mapping requirements for RDKit reactions to fire
while preventing charge carry-over from azide to triazole.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


def test_smarts_mapping(smarts_pattern, pattern_name):
    """
    Test a SMARTS pattern and report atom mapping details.
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
        return
    
    # Parse reactants
    azide = Chem.MolFromSmiles("N=[N+]=[N-]", sanitize=False)
    alkyne = Chem.MolFromSmiles("c1ccccc1C#C", sanitize=False)
    
    if azide is None or alkyne is None:
        print("Reactant parsing: FAILED")
        return
    
    print("Reactants parsed: SUCCESS")
    
    # Run reaction
    try:
        products = rxn.RunReactants((azide, alkyne))
        print("RunReactants(): SUCCESS")
        print("Number of product tuples:", len(products))
    except Exception as e:
        print("RunReactants(): FAILED -", str(e))
        return
    
    if not products:
        print("No products generated")
        return
    
    # Get first product
    product = products[0][0]
    print("Product obtained: SUCCESS")
    print("Total atoms:", product.GetNumAtoms())
    
    # Print atom table with map numbers and charges
    print("\nPRODUCT ATOM TABLE:")
    print("Index | Symbol | MapNum | Charge | Degree | Aromatic")
    print("-" * 60)
    for i in range(product.GetNumAtoms()):
        atom = product.GetAtomWithIdx(i)
        map_num = atom.GetAtomMapNum()
        print("{:5} | {:6} | {:6} | {:6} | {:6} | {}".format(
            i,
            atom.GetSymbol(),
            map_num,
            atom.GetFormalCharge(),
            atom.GetDegree(),
            atom.GetIsAromatic()
        ))
    
    # Count mapped atoms
    mapped_atoms = [i for i in range(product.GetNumAtoms()) if product.GetAtomWithIdx(i).GetAtomMapNum() > 0]
    print("\nMapped atoms in product:", len(mapped_atoms))
    print("Mapped atom indices:", mapped_atoms)
    
    # Count charged atoms
    charged_atoms = [(i, product.GetAtomWithIdx(i).GetFormalCharge()) for i in range(product.GetNumAtoms()) if product.GetAtomWithIdx(i).GetFormalCharge() != 0]
    print("\nCharged atoms in product:", len(charged_atoms))
    for idx, charge in charged_atoms:
        atom = product.GetAtomWithIdx(idx)
        map_num = atom.GetAtomMapNum()
        print("  Atom {}: {}, charge={}, map_num={}".format(idx, atom.GetSymbol(), charge, map_num))


def main():
    print("\n")
    print("*" * 80)
    print("ATOM MAPPING INVESTIGATION")
    print("*" * 80)
    print("\n")
    
    # Original SMARTS from pipeline
    original_cuaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"
    original_ruaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n:3]-[n+:2]=[n-:1]-1"
    
    # Unmapped triazole nitrogens SMARTS (what I tried)
    unmapped_cuaac = "N=[N+]=[N-].[C:4]#[C:5]>>[C:4]1=[C:5]-[n]-[n]=[n]-1"
    unmapped_ruaac = "N=[N+]=[N-].[C:4]#[C:5]>>[C:5]1=[C:4]-[n]-[n]=[n]-1"
    
    # Partially unmapped (keep carbon maps, remove nitrogen maps)
    partial_cuaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n]-[n]=[n]-1"
    partial_ruaac = "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4]-[n]-[n]=[n]-1"
    
    print("\n" + "#" * 80)
    print("# ORIGINAL SMARTS (with mapped product nitrogens)")
    print("#" * 80)
    
    test_smarts_mapping(original_cuaac, "Original CuAAC")
    test_smarts_mapping(original_ruaac, "Original RuAAC")
    
    print("\n" + "#" * 80)
    print("# UNMAPPED TRIAZOLE NITROGENS SMARTS (what I tried)")
    print("#" * 80)
    
    test_smarts_mapping(unmapped_cuaac, "Unmapped CuAAC")
    test_smarts_mapping(unmapped_ruaac, "Unmapped RuAAC")
    
    print("\n" + "#" * 80)
    print("# PARTIALLY UNMAPPED SMARTS (keep carbon maps, remove nitrogen maps)")
    print("#" * 80)
    
    test_smarts_mapping(partial_cuaac, "Partial CuAAC")
    test_smarts_mapping(partial_ruaac, "Partial RuAAC")
    
    print("\n" + "=" * 80)
    print("INVESTIGATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
