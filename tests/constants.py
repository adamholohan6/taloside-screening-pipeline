"""Shared test constants for taloside pipeline tests."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCAFFOLD = (
    "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)"
    "C4=C([N+]([O-])=O)C=CC=C4"
)

BUILDING_BLOCKS_MINI = [
    {"id": "BB-001-Ph", "smiles": "C#Cc1ccccc1"},
    {"id": "BB-002-4OMe", "smiles": "C#Cc1ccc(OC)cc1"},
]

BUILDING_BLOCKS_PHASE2 = [
    {"id": "BB-001-Ph", "smiles": "C#Cc1ccccc1"},
    {"id": "BB-002-4OMe", "smiles": "C#Cc1ccc(OC)cc1"},
    {"id": "BB-003-4Cl", "smiles": "C#Cc1ccc(Cl)cc1"},
    {"id": "BB-004-4F", "smiles": "C#Cc1ccc(F)cc1"},
    {"id": "BB-005-3Br", "smiles": "C#Cc1cc(Br)ccc1"},
    {"id": "BB-006-2NO2", "smiles": "C#Cc1ccccc1[N+](=O)[O-]"},
    {"id": "BB-007-Pyridine", "smiles": "C#Cc1ccccn1"},
    {"id": "BB-008-Furan", "smiles": "C#Cc1ccoc1"},
]

LIBRARY_COLUMNS = {
    "compound_id",
    "parent_scaffold_id",
    "building_block_id",
    "regioisomer",
    "product_smiles",
    "product_inchi",
    "product_inchikey",
    "molecular_weight",
    "h_donors",
    "h_acceptors",
    "logp",
    "tpsa",
    "rotatable_bonds",
    "generation_status",
}

PHASE2_EXPORT_FILES = [
    "01_all_generated_compounds.csv",
    "02_lipinski_passed.csv",
    "03_lipinski_failed.csv",
    "04_lipinski_clean_no_pains.csv",
    "05_lipinski_with_pains.csv",
    "06_pains_undetermined.csv",
    "07_lead_scored.csv",
]
