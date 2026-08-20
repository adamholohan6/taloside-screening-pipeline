# Phase 3 docking inputs (Galectin-3, PDB 3ZSJ)

Place prepared receptor files here before running Phase 3:

| File | Description |
|------|-------------|
| `3ZSJ.pdb` | Crystal structure (for lactose redock validation) |
| `3ZSJ.pdbqt` | Prepared receptor in AutoDock PDBQT format |
| `3ZSJ_clean.pdbqt` | Apo receptor used for the clean Phase 3 rerun |
| `7RGX.pdb` | Crystal structure carrying the 57I ligand (pyranose-ring validation) |
| `7RGX_57I_docked.pdbqt` | Redocked 57I pose (tracked in this repository) |

`data/` is otherwise excluded from version control. The crystal structures are
public RCSB entries, so fetch them directly:

```bash
curl -o data/docking/3ZSJ.pdb https://files.rcsb.org/download/3ZSJ.pdb
curl -o data/docking/7RGX.pdb https://files.rcsb.org/download/7RGX.pdb
```

`7RGX_57I_docked.pdbqt` is the one exception: it is a Phase 3 output rather than
a public file, and it is committed so that the pyranose-ring RMSD reported in the
manuscript can be checked without installing AutoDock Vina or rerunning docking.

## Quick start

```bash
# After Phase 2
python -m taloside_pipeline.phase2_integration

# Phase 3 (requires Vina on PATH)
python -m taloside_pipeline.phase3_docking
```

Default binding site (CRD): center (10, 15, 5) Å, box 20 Å. Outputs go to `phase3_output/`, including `08_docking_results.csv`.

For the clean rerun, the receptor PDBQT was regenerated with `scripts/clean_receptor_pdbqt.py`, and docking outputs were written to `phase3_output_clean/`.

## Validation

Once the inputs above are in place, the 57I pyranose-ring RMSD reported in the
manuscript (§3.4) can be reproduced with:

```bash
python scripts/validation/compute_57i_pyranose_rmsd.py
```
