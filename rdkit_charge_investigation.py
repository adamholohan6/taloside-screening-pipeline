"""
Focused RDKit Investigation: Azide-Alkyne Cycloaddition Charge Preservation

Investigates:
A) Standard RDKit reaction SMARTS for azide-alkyne cycloaddition
B) Can formal charges be reset in product template despite atom mapping?
C) Do RDKit reactions require post-reaction normalization for cycloadditions?
D) Minimal test: azide + alkyne with neutral product template
"""

from rdkit import Chem
from rdkit.Chem import AllChem


def test_minimal_azide_alkyne():
    """
    Minimal test: [N:1]=[N+:2]=[N-:3] + C#C with neutral triazole product template.
    """
    print("=" * 80)
    print("MINIMAL TEST: Azide + Alkyne with Neutral Product Template")
    print("=" * 80)
    
    # Reactants
    azide_smiles = "[N:1]=[N+:2]=[N-:3]"
    alkyne_smiles = "C#C"
    
    print("\nReactants:")
    print("  Azide:", azide_smiles)
    print("  Alkyne:", alkyne_smiles)
    
    # Parse reactants with atom maps
    azide = Chem.MolFromSmiles(azide_smiles, sanitize=False)
    alkyne = Chem.MolFromSmiles(alkyne_smiles, sanitize=False)
    
    if azide is None or alkyne is None:
        print("FAILED to parse reactants")
        return
    
    print("\nAzide atom properties (BEFORE reaction):")
    print("Index | Symbol | Charge | Degree")
    print("-" * 40)
    for i in range(azide.GetNumAtoms()):
        atom = azide.GetAtomWithIdx(i)
        print("{:5} | {:6} | {:6} | {:6}".format(
            i,
            atom.GetSymbol(),
            atom.GetFormalCharge(),
            atom.GetDegree()
        ))
    
    # Test different product templates
    product_templates = [
        ("Charged aromatic (current)", "[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"),
        ("Neutral aromatic", "[C:4]1=[C:5]-[n:1]-[n:2]=[n:3]-1"),
        ("Charged non-aromatic", "[C:4]1=[C:5]-[N:1]-[N+:2]=[N-:3]-1"),
        ("Neutral non-aromatic", "[C:4]1=[C:5]-[N:1]-[N:2]=[N:3]-1"),
        ("Explicit charge reset", "[C:4]1=[C:5]-[n:1]-[n:2:0]=[n:3:0]-1"),
    ]
    
    for template_name, product_template in product_templates:
        print("\n" + "=" * 80)
        print("Testing:", template_name)
        print("Product template:", product_template)
        print("=" * 80)
        
        # Build SMARTS
        smarts = azide_smiles + "." + alkyne_smiles + ">>" + product_template
        print("\nFull SMARTS:", smarts)
        
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            print("Reaction built: SUCCESS")
        except Exception as e:
            print("Reaction built: FAILED -", str(e))
            continue
        
        # Run reaction
        try:
            products = rxn.RunReactants((azide, alkyne))
            print("RunReactants(): SUCCESS")
            print("Number of product tuples:", len(products))
        except Exception as e:
            print("RunReactants(): FAILED -", str(e))
            continue
        
        if not products:
            print("No products generated")
            continue
        
        # Get first product
        product = products[0][0]
        
        print("\nProduct atom properties (AFTER RunReactants):")
        print("Index | Symbol | Charge | Degree | Aromatic")
        print("-" * 55)
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
            
            # Print charges after sanitization
            print("\nAtom charges AFTER sanitization:")
            print("Index | Symbol | Charge")
            print("-" * 30)
            for i in range(product.GetNumAtoms()):
                atom = product.GetAtomWithIdx(i)
                print("{:5} | {:6} | {:6}".format(
                    i,
                    atom.GetSymbol(),
                    atom.GetFormalCharge()
                ))
        except Exception as e:
            print("Sanitization: FAILED -", str(e))


def test_post_reaction_normalization():
    """
    Test if post-reaction charge normalization can fix the issue.
    """
    print("\n" + "=" * 80)
    print("POST-REACTION NORMALIZATION TEST")
    print("=" * 80)
    
    # Use current SMARTS
    smarts = "[N:1]=[N+:2]=[N-:3].C#C>>[C:4]1=[C:5]-[n:1]-[n+:2]=[n-:3]-1"
    
    print("\nSMARTS:", smarts)
    
    # Build reaction
    rxn = AllChem.ReactionFromSmarts(smarts)
    
    # Parse reactants
    azide = Chem.MolFromSmiles("[N:1]=[N+:2]=[N-:3]", sanitize=False)
    alkyne = Chem.MolFromSmiles("C#C", sanitize=False)
    
    # Run reaction
    products = rxn.RunReactants((azide, alkyne))
    product = products[0][0]
    
    print("\nBefore normalization:")
    print("Index | Symbol | Charge")
    print("-" * 30)
    for i in range(product.GetNumAtoms()):
        atom = product.GetAtomWithIdx(i)
        print("{:5} | {:6} | {:6}".format(
            i,
            atom.GetSymbol(),
            atom.GetFormalCharge()
        ))
    
    # Attempt charge normalization
    print("\nAttempting charge normalization...")
    for i in range(product.GetNumAtoms()):
        atom = product.GetAtomWithIdx(i)
        atom.SetFormalCharge(0)
    
    print("\nAfter manual charge reset:")
    print("Index | Symbol | Charge")
    print("-" * 30)
    for i in range(product.GetNumAtoms()):
        atom = product.GetAtomWithIdx(i)
        print("{:5} | {:6} | {:6}".format(
            i,
            atom.GetSymbol(),
            atom.GetFormalCharge()
        ))
    
    # Attempt sanitization after charge reset
    print("\nAttempting sanitization after charge reset...")
    try:
        Chem.SanitizeMol(product, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        print("Sanitization: SUCCESS")
    except Exception as e:
        print("Sanitization: FAILED -", str(e))


def test_standard_rdkit_patterns():
    """
    Test if there's a standard RDKit pattern that works.
    """
    print("\n" + "=" * 80)
    print("STANDARD RDKIT PATTERN INVESTIGATION")
    print("=" * 80)
    
    # Common variations found in RDKit documentation/examples
    patterns = [
        ("Neutral azide + neutral alkyne", "[N:1]=[N:2]=[N:3].[C:4]#[C:5]>>[C:4]1=[C:5]-[n:1]-[n:2]=[n:3]-1"),
        ("No atom maps in product", "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C]1=[C]-[n]-[n]=[n]-1"),
        ("Explicit bond orders only", "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1[C:5][N:1][N:2][N:3]1"),
    ]
    
    for pattern_name, smarts in patterns:
        print("\n" + "-" * 80)
        print("Testing:", pattern_name)
        print("SMARTS:", smarts)
        print("-" * 80)
        
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            print("Reaction built: SUCCESS")
            
            azide = Chem.MolFromSmiles("[N:1]=[N+:2]=[N-:3]", sanitize=False)
            alkyne = Chem.MolFromSmiles("C#C", sanitize=False)
            
            products = rxn.RunReactants((azide, alkyne))
            
            if not products:
                print("No products generated")
                continue
            
            product = products[0][0]
            
            print("Product charges:")
            for i in range(product.GetNumAtoms()):
                atom = product.GetAtomWithIdx(i)
                if atom.GetSymbol() == "N":
                    print("  N{}: charge={}".format(i, atom.GetFormalCharge()))
            
            try:
                Chem.SanitizeMol(product, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
                print("Sanitization: SUCCESS")
            except Exception as e:
                print("Sanitization: FAILED -", str(e))
                
        except Exception as e:
            print("FAILED -", str(e))


def main():
    print("\n")
    print("*" * 80)
    print("RDKIT CHARGE PRESERVATION INVESTIGATION")
    print("*" * 80)
    print("\n")
    
    # D) Minimal test
    test_minimal_azide_alkyne()
    
    # C) Post-reaction normalization
    test_post_reaction_normalization()
    
    # A) Standard RDKit patterns
    test_standard_rdkit_patterns()
    
    print("\n" + "=" * 80)
    print("INVESTIGATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
