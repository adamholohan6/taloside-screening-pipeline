# Notebook Execution Summary

**Notebook:** `notebooks/Phase2_VirtualLibraryExpansion.ipynb`  
**Executed:** 2026-05-31  
**Method:** `py -3.14 -m jupyter nbconvert --to notebook --execute notebooks/Phase2_VirtualLibraryExpansion.ipynb --output Phase2_VirtualLibraryExpansion.ipynb --ExecutePreprocessor.timeout=600`

---

## Environment

| Component | Version |
|-----------|---------|
| Python | 3.14.x |
| Install command | `pip install -r requirements-dev.txt -e .` |
| Working directory | Repository root |

---

## Execution Result

| Metric | Value |
|--------|-------|
| Total code cells | 15 |
| Cells executed | 15 |
| Cells with errors | 0 |
| Output notebook size | ~238 KB (includes saved outputs and figures) |

All code cells have non-null `execution_count` values. The executed notebook overwrites the source file in place with outputs preserved.

---

## Pipeline Outputs Generated

### CSV files (`phase2_output/`)

| File | Rows | Status |
|------|------|--------|
| `01_all_generated_compounds.csv` | 16 | OK |
| `02_lipinski_passed.csv` | 14 | OK |
| `03_lipinski_failed.csv` | 2 | OK |
| `04_lipinski_clean_no_pains.csv` | 14 | OK |
| `05_lipinski_with_pains.csv` | 0 | OK (empty — header only) |
| `06_pains_undetermined.csv` | 0 | OK (empty — header only) |
| `07_lead_scored.csv` | 14 | OK |

### Figure files (`phase2_output/`)

| File | Description |
|------|-------------|
| `descriptor_distributions.png` | Histograms by regioisomer |
| `lipinski_breakdown.png` | Per-rule Lipinski pass counts |
| `lead_score_landscape.png` | MW vs LogP coloured by lead score |

---

## Key Workflow Metrics (from notebook summary cell)

| Metric | Value |
|--------|-------|
| Total enumerated | 16 |
| 1,4-CuAAC products | 8 |
| 1,5-RuAAC products | 8 |
| Lipinski pass rate | 87.5% (14/16) |
| Lipinski failed | 2 (2-nitro analogues, HBA threshold) |
| PAINS clean | 14 |
| PAINS flagged | 0 |
| PAINS undetermined | 0 |
| Top lead | `SCAF-001_BB-001-Ph_1` (score ≈ 0.873) |

These match the README Key Results and `python -m taloside_pipeline.phase2_integration` output.

---

## Import Resolution

The notebook resolves the repository root whether launched from:

- **Project root** — `src/` found directly; `sys.path` updated.
- **`notebooks/` directory** — parent `src/` detected automatically.

After `pip install -e .`, imports also work without manual `sys.path` manipulation, but the notebook retains path bootstrapping for Jupyter sessions that may not use the editable install kernel.

---

## Warnings Observed (non-blocking)

1. **Windows zmq/Proactor event loop** — RuntimeWarning during nbconvert kernel startup; execution completed normally.
2. **Empty PAINS CSV files** — pandas raises `EmptyDataError` when reading zero-row exports; expected for the current curated set.

---

## Verification Commands Used

```powershell
# Install
py -3.14 -m pip install -r requirements-dev.txt -e .

# Execute notebook
py -3.14 -m jupyter nbconvert --to notebook --execute notebooks/Phase2_VirtualLibraryExpansion.ipynb --output Phase2_VirtualLibraryExpansion.ipynb --ExecutePreprocessor.timeout=600

# Verify pipeline CLI (README quick-start)
py -3.14 -m taloside_pipeline.phase2_integration
```

Both the notebook and CLI completed successfully and wrote identical output file names to `phase2_output/`.
