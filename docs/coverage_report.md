# Test coverage report

**Generated:** 2026-08-20 (v2.2.0)
**Command:** `pytest --cov=src --cov-report=term-missing`
**Environment:** Python 3.14.5, RDKit 2026.03.2, Windows (win32)

This report is produced from an actual run. The previous version of this file
was dated 2026-05-31, predated the June 2026 bug fixes, and reported 40 tests.

---

## Summary

| Metric | Value |
|---|---|
| Tests collected | **53** |
| Tests passed | **53** |
| Tests failed | 0 |
| Overall statement coverage | **69%** |
| Statements | 1023 |
| Statements missed | 314 |

## Per-module coverage

| Module | Stmts | Miss | Cover |
|---|---:|---:|---:|
| `src/__init__.py` | 0 | 0 | 100% |
| `src/taloside_pipeline/__init__.py` | 8 | 0 | 100% |
| `src/taloside_pipeline/library_generator.py` | 1 | 0 | 100% |
| `src/taloside_pipeline/phase2_integration.py` | 211 | 24 | **89%** |
| `src/taloside_pipeline/glycolibrary_generator.py` | 249 | 72 | 71% |
| `src/taloside_pipeline/descriptor_calculator.py` | 65 | 21 | 68% |
| `src/taloside_pipeline/phase3_docking.py` | 489 | 197 | 60% |
| **TOTAL** | **1023** | **314** | **69%** |

## Test distribution by marker

| Marker | Count |
|---|---:|
| `unit` | 31 |
| `smoke` | 3 |
| `integration` | 3 |
| `slow` | 1 |

## Why `phase3_docking.py` sits lowest

At 60% it is the least-covered module, which is expected: most of its
uncovered lines are the AutoDock Vina and Open Babel subprocess paths. The unit
tests mock those subprocesses rather than invoking the external binaries, so
CI can run without Vina or Open Babel installed. The uncovered ranges are
concentrated in pose parsing, PDBQT preparation, and the docking loop itself.

## Superseded claims

| Claim | Where it appeared | Correct value |
|---|---|---|
| 40 tests | this file (2026-05-31 version) | 53 |
| 79% coverage | manuscript §2.10 | 69% |
| 88% of `phase2_integration.py` | manuscript §2.10 | 89% |

The README's figure of 53 tests was already correct.

## Reproducing

```bash
pip install -e .
pip install pytest pytest-cov
pytest --cov=src --cov-report=term-missing
```

CI runs the same suite on Python 3.10, 3.12 and 3.14 via
`.github/workflows/tests.yml` (excluding tests marked `slow`).
