# Coverage Report

**Generated:** 2026-05-31  
**Command:** `py -3.14 -m pytest --cov=src --cov-report=term-missing --cov-report=html`  
**Environment:** Python 3.14.5, pytest 9.0.3, pytest-cov 7.1.0, RDKit (via pip)

---

## Summary

| Metric | Before Priority 2 | After Priority 2 | Change |
|--------|-------------------|------------------|--------|
| **Total tests** | 11 | **40** | +29 |
| **Overall coverage** | 49% | **79%** | +30 pp |
| **`phase2_integration.py`** | 12% | **88%** | +76 pp |
| **`glycolibrary_generator.py`** | 72% | **74%** | +2 pp |
| **`descriptor_calculator.py`** | 67% | **67%** | — |

**Result:** 40 passed, 0 failed, 0 skipped (8.33 s)

---

## Coverage by Module

| Module | Statements | Missed | Coverage | Notes |
|--------|------------|--------|----------|-------|
| `src/__init__.py` | 0 | 0 | **100%** | Empty stub |
| `src/taloside_pipeline/__init__.py` | 7 | 0 | **100%** | Public API exports |
| `src/taloside_pipeline/library_generator.py` | 1 | 0 | **100%** | Compatibility re-export |
| `src/taloside_pipeline/phase2_integration.py` | 188 | 22 | **88%** | Lipinski, PAINS, lead score, pipeline |
| `src/taloside_pipeline/glycolibrary_generator.py` | 243 | 64 | **74%** | Library generation engine |
| `src/taloside_pipeline/descriptor_calculator.py` | 63 | 21 | **67%** | Phase 1 descriptor CLI |
| **TOTAL** | **502** | **107** | **79%** | |

HTML report: `htmlcov/index.html`  
JSON report: `coverage.json`

---

## Uncovered Lines (Priority Gaps)

### `phase2_integration.py` (22 lines missed)

| Lines | Function / context | Reason uncovered |
|-------|-------------------|------------------|
| 105, 109–114 | `PAINSFilter._get_catalog` | Catalog load failure path (requires broken RDKit install) |
| 139–141 | `PAINSFilter.screen_molecule` | Exception during catalog match |
| 241–242, 255–257 | `apply_pains_filter` | Sanitization fallbacks; PAINS-positive branch |
| 263–269 | `apply_pains_filter` | Generic exception handler |
| 376–377 | `run_phase2_pipeline` | Empty library early return |
| 478 | `if __name__ == "__main__"` | Not executed via import |

### `glycolibrary_generator.py` (64 lines missed)

Primarily error-handling branches, logging file handler setup, export helpers (`export_library`, `export_failed_products`), and edge cases in product sanitization/deduplication. Core generation path is covered by integration and SMARTS tests.

### `descriptor_calculator.py` (21 lines missed)

CLI entry point (`main()`, `if __name__`) and exception handler in `validate_smiles`. Functional API is covered by existing unit tests.

---

## Test Suite Breakdown

| File | Tests | Markers | Focus |
|------|-------|---------|-------|
| `test_smarts_validation.py` | 12 | `unit` | SMARTS parsing, regioisomer labels, reaction execution |
| `test_phase2_integration.py` | 11 | `unit`, `integration`, `slow` | Lipinski, PAINS, lead score, full workflow, pipeline smoke |
| `test_smoke.py` | 6 | `unit`, `smoke` | Package API, CLI entry points |
| `test_descriptor_calculator.py` | 6 | — | Phase 1 descriptors |
| `test_library_generator.py` | 5 | — | Library generator basics |

### Marker usage

```bash
pytest -m unit          # 32 tests
pytest -m integration   # 3 tests
pytest -m smoke         # 3 tests
pytest -m slow          # 1 test (run_phase2_pipeline_smoke)
pytest -m "not slow"    # 39 tests (fast CI subset)
```

---

## Coverage vs Audit Targets

| Audit recommendation | Status |
|---------------------|--------|
| SMARTS validation tests | **Done** — 12 tests in `test_smarts_validation.py` |
| Phase 2 integration tests | **Done** — 11 tests in `test_phase2_integration.py` |
| PAINSFilter unit tests | **Done** — clean, undetermined, invalid SMILES |
| `compute_lead_scores` tests | **Done** — formula, sorting, single-compound edge case |
| End-to-end smoke test | **Done** — `test_run_phase2_pipeline_smoke` + CLI smoke |
| `phase2_integration` >80% coverage | **Done** — 88% |

---

## Recommended Next Steps (Priority 3+)

1. Add test for PAINS-positive compound (mock catalog match) to cover lines 255–257.
2. Add test for catalog load failure via monkeypatch on `FilterCatalog.FilterCatalog`.
3. Add CLI test for `taloside-descriptors` / `descriptor_calculator.main()`.
4. Cover `export_library` / `export_failed_products` helpers if CSV export via generator is part of public API.

---

*Generated after Priority 2 test implementation.*
