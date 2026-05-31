# Notebook Alignment Report

**Date:** 2026-05-31  
**Scope:** Priority 1 fixes from `audit_report.md`

---

## 1. Duplicate Notebook Removed

| Action | Path |
|--------|------|
| **Deleted** | `src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb` |
| **Retained (canonical)** | `notebooks/Phase2_VirtualLibraryExpansion.ipynb` |

The duplicate was byte-identical to the notebooks copy. Only the `notebooks/` path is referenced in README and is the single source of truth.

---

## 2. README Installation & Documentation

### Changes made

- Added `pip install -e .` after runtime dependency installation.
- Documented dev/notebook setup via `requirements-dev.txt`.
- Added Windows note to use `python3` / `py -3` when `python` points to a legacy interpreter.
- Documented console entry point `taloside-phase2`.
- Expanded output table to include all seven CSV files (`01`–`07`).
- Noted that the notebook runs the same workflow and writes to `phase2_output/`.

### Verified commands (Python 3.14, after `pip install -e .`)

| Command | Result |
|---------|--------|
| `py -3.14 -m taloside_pipeline.phase2_integration` | **Success** — 16 compounds generated, 7 CSVs written |
| `py -3.14 -m jupyter nbconvert --execute notebooks/Phase2_VirtualLibraryExpansion.ipynb` | **Success** — all 15 code cells executed |

> **Note:** On Windows cp1252 consoles, `phase2_integration.py` logs containing Unicode arrows (`→`) may emit harmless `UnicodeEncodeError` logging warnings. Pipeline execution still completes successfully.

---

## 3. Notebook Alignment with Canonical Pipeline

The notebook was rewritten from scratch. All library-generation and filtering logic is imported from `taloside_pipeline`; duplicate inline implementations were removed.

### Removed (legacy notebook behaviour)

| Removed item | Reason |
|--------------|--------|
| `ReactionSMARTS.TRIAZOLE_FORMATION` | Deprecated; mixed unlabelled regioisomers |
| 20-building-block library | Did not match curated Phase 2 set |
| Inline `apply_lipinski_filter()` | Duplicated `phase2_integration.py` |
| Inline `rank_lead_compounds()` | Different formula from `compute_lead_scores()` |
| `GlycoLibraryGenerator` direct instantiation | Replaced by `generate_triazole_library()` |
| Separate `phase2_notebook_output/` directory | Unified to `phase2_output/` |
| Objective text “500–1000 compounds” | Replaced with accurate 16-compound workflow |
| Hard-coded 5-compound MSc subset | Replaced by `load_compounds_from_dict()` (10 compounds) |

### Added / aligned

| Step | Implementation |
|------|----------------|
| Library generation | `generate_triazole_library()` — CuAAC + RuAAC |
| Building blocks | Same 8 blocks as `run_phase2_pipeline()` |
| Configuration | Same `LibraryConfig` thresholds and `phase2_output/` |
| Lipinski filter | `apply_lipinski_filter(library_df, strict_mode=True)` |
| PAINS screening | `apply_pains_filter(lipinski_passed)` |
| Lead scoring | `compute_lead_scores(lipinski_passed)` + PAINS status merge |
| CSV export | All seven files (`01`–`07`) with identical names |
| MSc comparison | `load_compounds_from_dict()` from `descriptor_calculator.py` |

### Visualisation-only code (not pipeline duplicates)

Three plotting cells remain for exploratory analysis. They consume pipeline outputs and do not alter filtering logic:

- Descriptor histograms by regioisomer
- Lipinski rule pass counts (uses same threshold constants for labelling only)
- Lead-score MW vs LogP scatter plot

Figures are saved to `phase2_output/` as PNG files.

---

## 4. Dependency Updates

Added to `requirements-dev.txt`:

```
matplotlib>=3.5.0
seaborn>=0.12.0
jupyter>=1.0.0
ipykernel>=6.0.0
```

Runtime dependencies (`requirements.txt`) unchanged — notebook deps are dev-only.

---

## 5. Notebook Structure (post-alignment)

| Section | Content |
|---------|---------|
| 1. Setup & Imports | Standard + pipeline module imports |
| 2. Scaffold, Building Blocks & Configuration | Mirrors `phase2_integration.py` |
| 3. Run Canonical Pipeline Steps | generate → Lipinski → PAINS → score → export |
| 4. Descriptive Statistics | Summary tables |
| 5. Visualizations | Descriptor, Lipinski, lead-score plots |
| 6. Comparison | Original 10 MSc talosides vs generated library |
| 7. Summary | Funnel statistics matching pipeline log |
| Appendix | Environment versions for reproducibility |

**Cell count:** 23 total (15 code, 8 markdown) — reduced from the previous 40+ cell notebook with duplicated logic.

---

## 6. Expected Output Parity

After alignment, notebook and `python -m taloside_pipeline.phase2_integration` produce identical workflow results:

| Output | Expected rows |
|--------|---------------|
| `01_all_generated_compounds.csv` | 16 |
| `02_lipinski_passed.csv` | 14 |
| `03_lipinski_failed.csv` | 2 |
| `04_lipinski_clean_no_pains.csv` | 14 |
| `05_lipinski_with_pains.csv` | 0 (header only) |
| `06_pains_undetermined.csv` | 0 (header only) |
| `07_lead_scored.csv` | 14 |

Both paths write to `phase2_output/` at the repository root.

---

## 7. Files Modified

| File | Change |
|------|--------|
| `notebooks/Phase2_VirtualLibraryExpansion.ipynb` | Rewritten and executed with outputs |
| `src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb` | Deleted |
| `README.md` | Installation, outputs, entry points |
| `requirements-dev.txt` | Notebook dependencies |
| `notebook_alignment_report.md` | This report |
| `notebook_execution_summary.md` | Execution log |

---

*Priority 2+ items (tests, CI, coverage) intentionally deferred per task scope.*
