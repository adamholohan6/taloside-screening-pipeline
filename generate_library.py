"""
Taloside Screening Pipeline - Physicochemical Descriptor Calculator
-------------------------------------------------------------------
Calculates ADMET properties (Molecular Weight, LogP, H-donors, H-acceptors, 
TOSA, Rotatable Bonds) for a library of taloside compounds using RDKit.

Requires: RDKit, Pandas
"""
import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)


def validate_smiles(smiles_string):
    """
    Validate a SMILES string.
    
    Args:
        smiles_string (str): SMILES notation string
        
    Returns:
        rdkit.Chem.Mol or None: RDKit molecule object if valid, None otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return mol
    except Exception as e:
        logger.warning(f"SMILES validation error: {e}")
        return None


def calculate_descriptors(mol):
    """
    Calculate physicochemical descriptors for a molecule.
    
    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object
        
    Returns:
        dict: Dictionary containing calculated descriptors
    """
    return {
        "Molecular_Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "H_Donors": Descriptors.NumHDonors(mol),
        "H_Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol)
    }


def process_compounds(compounds_dict):
    """
    Process a dictionary of compounds and calculate their descriptors.
    
    Args:
        compounds_dict (dict): Dictionary with compound names as keys and SMILES as values
        
    Returns:
        tuple: (pd.DataFrame, int) - DataFrame with results and count of failed compounds
    """
    results = []
    failed_compounds = []
    
    for name, smiles in compounds_dict.items():
        mol = validate_smiles(smiles)
        
        if mol is not None:
            descriptors = calculate_descriptors(mol)
            results.append({
                "Compound": name,
                "SMILES": smiles,
                **descriptors
            })
            logger.info(f"✓ Processed {name}")
        else:
            logger.warning(f"✗ Could not process {name}. Invalid SMILES string.")
            failed_compounds.append(name)
    
    if failed_compounds:
        logger.warning(f"Failed to process: {', '.join(failed_compounds)}")
    
    return pd.DataFrame(results), len(failed_compounds)


def load_compounds_from_dict():
    """
    Returns the library of taloside compounds.
    
    Returns:
        dict: Compound names mapped to SMILES strings
    """
    return {
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


def main():
    """Command-line interface for descriptor calculation."""
    parser = argparse.ArgumentParser(
        description="Calculate physicochemical descriptors for taloside compounds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_library.py
  python generate_library.py -o my_descriptors.csv
        """
    )
    
    parser.add_argument(
        "-o", "--output",
        default="taloside_descriptors.csv",
        help="Output CSV filename (default: taloside_descriptors.csv)"
    )
    
    args = parser.parse_args()
    
    try:
        logger.info("Loading taloside compound library...")
        compounds = load_compounds_from_dict()
        
        logger.info(f"Processing {len(compounds)} compounds...")
        df, failed_count = process_compounds(compounds)
        
        if df.empty:
            logger.error("No compounds were successfully processed.")
            sys.exit(1)
        
        # Save to CSV
        output_path = Path(args.output)
        df.to_csv(output_path, index=False)
        
        logger.info(f"\n✓ SUCCESS! Processed {len(df)} compounds.")
        logger.info(f"Data saved to: {output_path.absolute()}")
        
        if failed_count > 0:
            logger.warning(f"Note: {failed_count} compound(s) failed to process.")
        
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
