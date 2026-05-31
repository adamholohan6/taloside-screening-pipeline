# ChemRxiv Readiness Audit Report

**Repository:** taloside-screening-pipeline (Phase 2)  
**Audit date:** 2026-05-31  
**Scope:** Broken imports, stale tests, duplicate files, notebook execution, packaging/setup.py, README consistency, RDKit compatibility, unused files  
**Method:** Static code review, documentation cross-check, file-hash comparison, and limited runtime checks (`py -3.14` import test; full pytest/pipeline run blocked by missing dev dependencies in the active environment)

---

## Executive Summary

The core Python modules (`glycolibrary_generator.py`, `phase2_integration.py`, `descriptor_calculator.py`) are internally consistent and import cleanly when `src/` is on `PYTHONPATH`. However, the **notebook, README, and CONTRIBUTING docs describe a different workflow than the corrected pipeline**, and the README quick-start command **fails without an editable install**. These gaps are the primary blockers for ChemRxiv supplementary-material reproducibility.

---

## 1. Critical Issues

### C1. Notebook uses deprecated reaction SMARTS (contradicts corrected pipeline)

`notebooks/Phase2_VirtualLibraryExpansion.ipynb` instantiates `GlycoLibraryGenerator` with `ReactionSMARTS.TRIAZOLE_FORMATION["smarts"]`, which is explicitly marked **DEPRECATED** in `glycolibrary_generator.py` (produces unlabelled mixed regioisomers). The canonical pipeline uses `generate_triazole_library()` with separate CuAAC/RuAAC templates.

**Impact:** Notebook results will not match README/Zenodo claims (16 tagged regioisomers, 87.5% Lipinski pass rate) or `phase2_integration.py` outputs.

**Location:** Notebook cell “Instantiate generator” (~line 269 in `.ipynb` source)

---

### C2. README quick-start fails without package installation

README instructs:

```bash
python -m taloside_pipeline.phase2_integration
```

after only `pip install -r requirements.txt`. This **does not install the package** into the environment.

**Verified:** `py -3.14 -m taloside_pipeline.phase2_integration` → `ModuleNotFoundError: No module named 'taloside_pipeline'`.

**Working alternatives (not documented):**
- `pip install -e .` then `python -m taloside_pipeline.phase2_integration`
- `PYTHONPATH=src python -m taloside_pipeline.phase2_integration` (Windows: `$env:PYTHONPATH="src"`)

---

### C3. Notebook workflow diverges from published Phase 2 pipeline

| Aspect | `phase2_integration.py` / README | Notebook |
|--------|----------------------------------|----------|
| Building blocks | 8 curated blocks | 20 blocks |
| Regioisomer handling | Both 1,4-CuAAC and 1,5-RuAAC, tagged | Single deprecated SMARTS, unlabelled |
| Library size | 16 compounds (8 × 2) | ~20 (objective text claims 500–1000) |
| PAINS screening | Yes (`FilterCatalog`, exports `06_pains_undetermined.csv`) | **Not implemented** |
| Lead scoring | `compute_lead_scores()` (0.4 TPSA + 0.3 MW + 0.2 LogP + 0.1 RotBonds) | Inline `rank_lead_compounds()` (different weights/targets) |
| Lipinski filter | Module function returning 2 DataFrames | Inline copy returning 3 values |
| Original MSc set | 10 compounds in `descriptor_calculator.py` | Only 5 hard-coded in notebook |
| Output directory | `phase2_output/` | `phase2_notebook_output/` |

**Impact:** ChemRxiv readers cannot reproduce paper/README numbers from the notebook alone.

---

### C4. Duplicate notebook creates maintenance and citation risk

Identical file (SHA-256 `6B1FFC302B982E0D522C1E2A4707824C39A57A9AF46E9E2AC5D6066955DD2156`) exists at:

- `notebooks/Phase2_VirtualLibraryExpansion.ipynb` (referenced in README)
- `src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb` (not in README project tree)

**Impact:** Future edits may update one copy only; supplementary material path becomes ambiguous.

---

## 2. Major Issues

### M1. Stale / incomplete test suite

| Module | Approx. coverage (from `htmlcov/`) | Tests |
|--------|-------------------------------------|-------|
| `descriptor_calculator.py` | 67% | 6 tests |
| `glycolibrary_generator.py` | 72% | 5 tests |
| `phase2_integration.py` | **12%** (162/184 statements missed) | **None** |
| `library_generator.py` | 100% | 1 import alias test |

**Gaps:**
- No tests for `PAINSFilter`, `apply_lipinski_filter`, `apply_pains_filter`, `compute_lead_scores`, or `run_phase2_pipeline`
- No integration/smoke test for end-to-end CSV export
- `pytest.ini` defines `slow`, `integration`, `unit` markers but no test uses them
- Prior run achieved **49% total coverage**; `.pytest_cache/v/cache/lastfailed` is empty (no recorded failures, but coverage is thin on critical path)

**Note:** `pytest` is not installed in the default `py -3.14` environment; dev dependencies must be installed to run the suite.

---

### M2. Notebook dependencies missing from packaging files

The notebook imports **matplotlib**, **seaborn**, and implicitly requires **jupyter/ipykernel**. None appear in:

- `requirements.txt`
- `setup.py` `install_requires`
- `requirements-dev.txt`

**Impact:** Clean environment reproducibility fails for the primary visualization artifact cited in README.

---

### M3. CONTRIBUTING.md references removed / renamed files

| CONTRIBUTING says | Actual codebase |
|-------------------|-----------------|
| Edit `generate_library.py` → `load_compounds_from_dict()` | File does not exist; function is in `descriptor_calculator.py` |
| Run `python generate_library.py` | No such script; use `taloside-descriptors` entry point or `python -m taloside_pipeline.descriptor_calculator` |
| Python 3.6+ | README/setup require **Python 3.8+** |

---

### M4. README output documentation incomplete and partially stale

README table lists outputs `01`–`05` only. `phase2_integration.py` also exports:

- `06_pains_undetermined.csv`
- `07_lead_scored.csv`

README Key Results mention PAINS undetermined count but do not document file `06`. Lead scoring output `07` is undocumented.

README project structure omits:

- `data/input/`
- `notebooks/` duplicate path under `src/taloside_pipeline/`
- `LICENSE`, `CONTRIBUTING.md`

---

### M5. Notebook not executed / no pinned outputs

All code cells have `"execution_count": null` and empty `"outputs"`. The appendix records **Python 3.14.5** but provides no executed evidence that cells run cleanly on the stated stack.

**ChemRxiv concern:** Supplementary notebooks should either include executed outputs with pinned versions or a verified execution script/log.

---

### M6. No CI or reproducibility automation

No `.github/workflows/`, no `environment.yml`/`conda-lock`, no `CITATION.cff`, no pinned lockfile. Reproducibility relies entirely on minimum-version specifiers (`rdkit>=2022.09.1`, etc.).

---

### M7. Orphan reference data

`data/input/taloside_descriptors.csv` (10 compounds) is **never read by any module**. It appears to be example/reference output but is not wired into tests or docs. Descriptor values may drift from `descriptor_calculator.py` if SMILES or RDKit version changes (no test validates CSV against live calculation).

---

### M8. Windows / environment footgun for `python` command

On the audit machine, default `python` resolves to **MGLTools Python 2.5/2.7** (`C:\Program Files (x86)\MGLTools-1.5.7\python.exe`), which raises `SyntaxError` on type-hinted source. README uses generic `python` without version guidance.

---

## 3. Minor Issues

### m1. Unused / redundant imports (no runtime break, lint noise)

| File | Unused |
|------|--------|
| `glycolibrary_generator.py` | `field` (dataclasses), `inchi` |
| `phase2_integration.py` | `numpy as np`, `Descriptors`, `Lipinski`, `Crippen`, `Dict` |
| `phase2_integration.py` | `GlycoLibraryGenerator`, `ReactionSMARTS` imported but unused in module body |

---

### m2. `phase2_integration.py` standalone fallback import is fragile

```python
except ImportError:
    from glycolibrary_generator import ...
```

Only works when executed as a script with CWD inside `src/taloside_pipeline/`, not from project root. Relative package import path is the supported route.

---

### m3. `setup.py` metadata inconsistencies

- `Development Status :: 3 - Alpha` vs public **v1.0.0** and Zenodo DOI
- `license="MIT"` string form is deprecated in modern setuptools (prefer `license_files`)
- No `package_data` / `MANIFEST.in` to include notebooks or `data/input/` in sdist
- Console entry points `taloside-descriptors` and `taloside-phase2` defined but not mentioned in README

---

### m4. Descriptor calculation inconsistency across modules

- `glycolibrary_generator.py` applies `MolStandardize` before descriptors
- `descriptor_calculator.py` and notebook comparison cells use raw `Descriptors.*` / `Crippen.MolLogP` without standardization

Minor scientific inconsistency when comparing Phase 1 (10 compounds) vs Phase 2 library.

---

### m5. RDKit PAINS catalog loading is version-dependent

`PAINSFilter._get_catalog()` adds **only the first matching** catalog constant (`PAINS`, `PAINS_A`, or `FILTER_PAINS`). On builds where only `PAINS_A` exists, screening may be a subset of advertised “PAINS_A/B/C”. Defensive coding is present but behavior varies by RDKit build.

---

### m6. Test side effects

`test_generate_triazole_library_returns_expected_columns` writes to `Path("test_output")` in the project root (not `tmp_path`), potentially leaving artifacts.

---

### m7. Repository artifacts present in working tree

`htmlcov/` and `.pytest_cache/` exist locally. Both are `.gitignore`d (good) but should not be published in ChemRxiv zip if present.

---

### m8. Documentation placeholders

- README Contact: “ORCID (add your ORCID if you have one)”
- Notebook appendix cites git commit but does not capture one programmatically

---

### m9. Empty `src/__init__.py`

Contains only a docstring. Harmless with `package_dir={"": "src"}`, but adds no value.

---

## 4. Broken Imports Summary

| Category | Finding | Severity |
|----------|---------|----------|
| Python package imports (`taloside_pipeline.*`) | **OK** when `src/` on path or after `pip install -e .` | — |
| README `python -m taloside_pipeline.phase2_integration` | **Broken** without install/`PYTHONPATH` | Critical |
| Default Windows `python` (MGLTools 2.x) | **SyntaxError** on modern source | Major (environment) |
| `CONTRIBUTING.md` → `generate_library.py` | **Broken reference** (file missing) | Major (docs) |
| Notebook → `matplotlib`, `seaborn` | **Missing from requirements** | Major |
| `phase2_integration` script fallback import | **Fragile** outside package dir | Minor |

No broken intra-package relative imports were found in the Python source tree.

---

## 5. Duplicate Files

| File A | File B | Relationship |
|--------|--------|--------------|
| `notebooks/Phase2_VirtualLibraryExpansion.ipynb` | `src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb` | **Byte-identical duplicate** |

No other duplicate Python modules detected. `library_generator.py` is an intentional re-export wrapper, not a duplicate implementation.

---

## 6. Unused Files

| Path | Status |
|------|--------|
| `data/input/taloside_descriptors.csv` | Present; **not referenced** by code or tests |
| `src/__init__.py` | Minimal stub; **non-functional** |
| `src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb` | Duplicate of notebooks copy; **redundant** if `notebooks/` is canonical |
| `htmlcov/`, `.pytest_cache/` | Local test artifacts; should not ship |
| `library_generator.py` | **Used** (backward-compat imports/tests) — not unused |

---

## 7. RDKit Compatibility Notes

| Area | Assessment |
|------|------------|
| Minimum version `>=2022.09.1` | Appropriate for `MolStandardize`, `FilterCatalog`, InChIKey dedup |
| `MolStandardize` API usage | Standard for 2022.x+; may fail on very old RDKit (correctly documented) |
| Triazole sanitization (`SANITIZE_ALL ^ SANITIZE_KEKULIZE`) | Reasonable workaround for aromatic triazoles |
| InChIKey deduplication | Version-robust approach (good) |
| PAINS `FilterCatalog` enum probing | Good defensive pattern; catalog breadth varies by build (see m5) |
| Python 3.14 in notebook metadata | RDKit import succeeded on audit machine (`py -3.14`); not universally guaranteed—pin environment for ChemRxiv |
| `rdkit` PyPI package name in `setup.py` | Correct for pip installs; conda users may need separate instructions |

---

## 8. Recommended Fixes (Priority Order)

### Priority 1 — Reproducibility blockers (before ChemRxiv submission)

1. **Align notebook with `phase2_integration.py`**
   - Replace `TRIAZOLE_FORMATION` with `generate_triazole_library()`
   - Use the same 8 building blocks as the pipeline
   - Import and call `apply_lipinski_filter`, `apply_pains_filter`, `compute_lead_scores` from `phase2_integration` instead of inline duplicates
   - Execute all cells; save outputs; pin Python/RDKit versions in a new `environment.yml` or `requirements-lock.txt`

2. **Fix README installation and quick-start**
   - Add `pip install -e .` (or document `PYTHONPATH=src`)
   - Document all seven output CSV files (`01`–`07`)
   - Mention console scripts: `taloside-phase2`, `taloside-descriptors`
   - Specify `python3.8+` explicitly; warn against legacy `python` on Windows

3. **Remove duplicate notebook**
   - Keep single canonical copy under `notebooks/`; delete `src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb`

### Priority 2 — Test and documentation integrity

4. **Add `tests/test_phase2_integration.py`**
   - Smoke test: `run_phase2_pipeline()` produces 16 compounds, expected CSV columns, Lipinski pass rate ~87.5%
   - Unit tests for `PAINSFilter`, `compute_lead_scores`, undetermined PAINS routing

5. **Update `CONTRIBUTING.md`**
   - Point to `descriptor_calculator.py` for compound additions
   - Align Python version to 3.8+
   - Replace `generate_library.py` references with current entry points

6. **Extend dependency files**
   - Add notebook extras: `matplotlib`, `seaborn`, `jupyter`
   - Consider `[project.optional-dependencies] notebook = [...]` in future `pyproject.toml`

### Priority 3 — Packaging and publication polish

7. **Add CI** (GitHub Actions: install deps, run pytest, optionally execute notebook with `nbconvert --execute`)

8. **Wire or document `data/input/taloside_descriptors.csv`**
   - Either add test comparing CSV to `calculate_library_descriptors()`, or document as frozen reference output

9. **Clean metadata**
   - Set setuptools development status to stable/4 if v1.0.0 is final
   - Add `CITATION.cff` matching Zenodo DOI
   - Fill ORCID in README

10. **Minor code hygiene**
    - Remove unused imports
    - Use `tmp_path` in library generator test
    - Add `MANIFEST.in` or `package_data` for notebooks/data if distributing via PyPI/sdist

---

## Appendix: Files Audited

```
setup.py, requirements.txt, requirements-dev.txt, pytest.ini
README.md, CONTRIBUTING.md, LICENSE, .gitignore
src/taloside_pipeline/*.py
src/taloside_pipeline/Phase2_VirtualLibraryExpansion.ipynb
notebooks/Phase2_VirtualLibraryExpansion.ipynb
tests/*.py
data/input/taloside_descriptors.csv
htmlcov/ (coverage artifacts from prior run)
```

---

*End of audit report. No source files were modified during this audit.*
