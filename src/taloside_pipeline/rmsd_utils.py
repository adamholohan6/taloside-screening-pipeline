from rdkit import Chem
from rdkit.Chem import rdMolAlign
import numpy as np

def get_pyranose_ring_atoms(mol: Chem.Mol):
    """Finds the C1, C2, C3, C4, C5, and O5 atoms of the pyranose ring."""
    # SMARTS for a 6-membered ring with one oxygen and five carbons
    ring_pat = Chem.MolFromSmarts("[C]1[C][C][C][C][O]1")
    matches = mol.GetSubstructMatches(ring_pat)
    if not matches:
        return None
    # For talosides, there should only be one such ring
    return matches[0]

def calculate_pyranose_rmsd(ref_mol: Chem.Mol, docked_mol: Chem.Mol):
    """Calculates RMSD between two molecules focusing only on the pyranose ring."""
    ref_atoms = get_pyranose_ring_atoms(ref_mol)
    dock_atoms = get_pyranose_ring_atoms(docked_mol)

    if not ref_atoms or not dock_atoms:
        return None

    # Create atom mapping for the ring
    # Since they are both the same scaffold core, we can assume the order matches
    # for the same SMARTS pattern match.
    atom_map = list(zip(dock_atoms, ref_atoms))

    # Calculate RMSD without alignment first (as they should be in the same coordinate frame)
    dock_coords = docked_mol.GetConformer().GetPositions()[list(dock_atoms)]
    ref_coords = ref_mol.GetConformer().GetPositions()[list(ref_atoms)]

    diff = dock_coords - ref_coords
    rmsd = np.sqrt((diff * diff).sum() / len(ref_atoms))

    return float(rmsd)
