"""
Test script to verify the Phase 3 fix works with a single compound
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

# Import the FIXED version
from taloside_pipeline.phase3_docking_FIXED import DockingConfig, VinaDocking
import pandas as pd

# Load a single failing compound from Phase 2 output
df = pd.read_csv("phase2_output/07_lead_scored.csv")
single_compound = df.iloc[[0]]  # First compound

print("=" * 70)
print("TESTING PHASE 3 FIX - SINGLE COMPOUND")
print("=" * 70)
print("Compound ID:", single_compound.iloc[0]['compound_id'])
print("SMILES:", single_compound.iloc[0]['product_smiles'])
print("=" * 70)

# Configure Phase 3
config = DockingConfig(
    receptor_pdb=Path("data/docking/3ZSJ.pdb"),
    receptor_pdbqt=Path("data/docking/3ZSJ.pdbqt"),
    output_dir=Path("phase3_output"),
)

# Run ligand preparation only
docking = VinaDocking(config)
prepared_df = docking.prepare_ligands(single_compound)

print("\n" + "=" * 70)
print("RESULT")
print("=" * 70)
print("Generation status:", prepared_df.iloc[0]['generation_status'])
print("Ligand PDBQT:", prepared_df.iloc[0]['ligand_pdbqt'])
print("=" * 70)
