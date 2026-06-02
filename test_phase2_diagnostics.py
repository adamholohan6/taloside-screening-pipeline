"""
Test script to regenerate Phase 2 library with diagnostic logging
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from taloside_pipeline.phase2_integration import run_phase2_pipeline

print("=" * 80)
print("REGENERATING PHASE 2 LIBRARY WITH DIAGNOSTIC LOGGING")
print("=" * 80)
print("\nThis will regenerate the triazole library with detailed logging to show:")
print("- Triazole ring atom properties before/after sanitization")
print("- SMILES export with and without kekuleSmiles parameter")
print("- When aromaticity is lost during processing")
print("=" * 80)

# Run Phase 2 pipeline
result = run_phase2_pipeline()

print("\n" + "=" * 80)
print("PHASE 2 COMPLETE")
print("=" * 80)
print("Check phase2_output/ for:")
print("- 01_all_generated_compounds.csv (new products)")
print("- library_output/glycolibrary_generator.log (diagnostic output)")
print("=" * 80)
