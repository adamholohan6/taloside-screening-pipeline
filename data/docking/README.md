# Phase 3 docking inputs (Galectin-3, PDB 3ZSJ)

Place prepared receptor files here before running Phase 3:

| File | Description |
|------|-------------|
| `3ZSJ.pdb` | Crystal structure (for lactose redock validation) |
| `3ZSJ.pdbqt` | Prepared receptor in AutoDock PDBQT format |

## Quick start

```bash
# After Phase 2
python -m taloside_pipeline.phase2_integration

# Phase 3 (requires Vina on PATH)
python -m taloside_pipeline.phase3_docking
```

Default binding site (CRD): center (10, 15, 5) Å, box 20 Å. Outputs go to `phase3_output/`, including `08_docking_results.csv`.
