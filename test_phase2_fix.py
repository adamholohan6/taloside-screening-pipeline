"""
Test script to verify Phase 2 fix by regenerating products with aromatic notation preserved
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from taloside_pipeline.phase2_integration import run_phase2_pipeline

print("=" * 80)
print("TESTING PHASE 2 FIX - REGENERATING PRODUCTS")
print("=" * 80)
print("\nRunning Phase 2 pipeline with kekuleSmiles=False fix...")
print("This should generate products with aromatic triazole notation [n] instead of kekulized [N+]=[N-]")
print("=" * 80)

# Run Phase 2 pipeline
result = run_phase2_pipeline()

print("\n" + "=" * 80)
print("PHASE 2 COMPLETE")
print("=" * 80)
print("Generated", len(result['all']), "total products")
print("Check phase2_output/01_all_generated_compounds.csv for new SMILES")
print("Look for aromatic [n] notation in triazole rings instead of [N+]=[N-]")
print("=" * 80)
