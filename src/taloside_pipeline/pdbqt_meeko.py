"""Meeko-based PDBQT export helpers.

This module isolates the Meeko integration so the docking pipeline can switch
from manual PDBQT writing to a robust, Vina-compatible writer.
"""

from __future__ import annotations

from typing import List

from rdkit import Chem


def mol_to_pdbqt_string_meeko(mol: Chem.Mol, conf_id: int = 0, name: str = "LIG") -> str:
    """Generate a Vina-compatible ligand PDBQT using Meeko.

    Parameters
    ----------
    mol
        RDKit molecule with at least one 3D conformer.
    conf_id
        Which conformer to export.
    name
        Ligand name used in REMARK lines.
    """
    try:
        from meeko import MoleculePreparation
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "meeko is required for PDBQT export. Install with `pip install meeko`."
        ) from e

    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformer for PDBQT export")

    # Keep only the requested conformer. Avoid reusing Conformer objects
    # across different RDKit Mol instances (can trigger atom mismatch).
    active_conf = mol.GetConformer(conf_id)
    conf_copy = Chem.Conformer(active_conf)
    mol.RemoveAllConformers()
    mol.AddConformer(conf_copy, assignId=True)


    prep = MoleculePreparation()
    setups = prep.prepare(mol)

    if not setups:
        raise ValueError("Meeko could not prepare molecule for PDBQT export")

    parts: List[str] = []
    for setup in setups:
        # meeko API varies across versions. In this environment,
        # RDKitMoleculeSetup exposes `write_coord_string()`.
        if hasattr(setup, "write_pdbqt_string"):
            try:
                parts.append(setup.write_pdbqt_string(name=name))
            except TypeError:
                parts.append(setup.write_pdbqt_string())
        elif hasattr(setup, "write_pdbqt"):
            raise RuntimeError(
                "Your meeko version exposes `write_pdbqt` but not `write_pdbqt_string`. "
                "Adapting this to a string requires a file writer path; update pdbqt_meeko.py."
            )
        elif hasattr(setup, "write_coord_string"):
            # Fallback: coord string + meeko does not provide full PDBQT string in this build.
            # We return coord string only is NOT a valid PDBQT for Vina.
            raise RuntimeError(
                "Your meeko version does not support generating a full PDBQT text string 'write_pdbqt_string'. "
                "It only provides `write_coord_string()`, which is insufficient for Vina."
            )
        else:
            raise RuntimeError(
                f"Unsupported meeko setup object (missing PDBQT writer): {type(setup)}"
            )



    return "".join(parts)

