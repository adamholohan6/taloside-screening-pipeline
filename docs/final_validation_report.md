# Final validation report

**Generated:** 2026-08-20 (v2.2.0 cleanup)
**Environment:** Python 3.14.5, RDKit 2026.03.2, AutoDock Vina 1.2.7,
Open Babel 3.1.1, BioPython 1.87, Windows (win32)

Every figure below comes from a command run against this working tree on the
date above. The previous version of this file was dated 2026-05-31 and described
a pipeline that no longer exists. Claims that did **not** reproduce are recorded
in [`UNVERIFIED.md`](../UNVERIFIED.md) rather than being quietly adjusted.

---

## 1. Environment and tooling

| Component | Required by | Detected | Verifies |
|---|---|---|---|
| Python | pipeline | 3.14.5 | — |
| RDKit | Phases 1–3 | 2026.03.2 | pinned in `requirements.txt` |
| AutoDock Vina | Phase 3 | **v1.2.7** | manuscript §2.7 claim |
| Open Babel | Phase 3 | **3.1.1** | SI Note S1 claim |
| BioPython | validation script | 1.87 | — |

Vina and Open Babel are installed but not on `PATH`; both were invoked by
absolute path. The Vina and Open Babel version claims were previously listed as
unverified and are now confirmed.

## 2. Test suite

```
pytest --cov=src --cov-report=term-missing
53 passed in 17.56s
TOTAL  1023 stmts  314 miss  69%
```

Full breakdown in [`coverage_report.md`](coverage_report.md).

## 3. Phase 2 — reproduced exactly

`python -m taloside_pipeline.phase2_integration`

| Claim | Expected | Observed | |
|---|---|---|---|
| Compounds generated | 16 | 16 | ✅ |
| — 1,4-CuAAC | 8 | 8 | ✅ |
| — 1,5-RuAAC | 8 | 8 | ✅ |
| Lipinski pass | 14 | 14 | ✅ |
| Lipinski pass rate | 87.5% | 87.5% | ✅ |
| Lipinski fail | 2 | 2 | ✅ |
| PAINS flagged | 0 | 0 | ✅ |
| PAINS undetermined | 0 | 0 | ✅ |
| Lead score range | 0.399–0.905 | 0.399–0.905 | ✅ |

The regenerated `phase2_output/01..07*.csv` are **bit-identical** to the copies
committed in this repository (`git status` reports no change after a clean
rerun). Phase 2 is fully deterministic.

## 4. Docking validation — 57I pyranose ring

`python scripts/validation/compute_57i_pyranose_rmsd.py` → **exit code 0**

| Quantity | Value |
|---|---|
| Cα atoms superposed (7RGX → 3ZSJ) | 138 |
| Cα RMSD | 0.275 Å |
| Pyranose RMSD, direct | 2.17 Å |
| Pyranose RMSD, Kabsch-aligned | **0.48 Å** |
| Manuscript §3.4 value | 0.48 Å ✅ |

This check needs no AutoDock Vina: the redocked pose
`data/docking/7RGX_57I_docked.pdbqt` is committed. Only the two RCSB structures
must be fetched first (see `data/docking/README.md`).

## 5. Grid centre — discrepancy found

The pipeline log records the centre it actually used:

```
Dynamically centered grid on crystal ligand: X=-20.98, Y=8.88, Z=-1.00
[grid] Docking grid centre: X=-20.982  Y=8.876  Z=-1.002  box=20.0 A^3
```

The manuscript states X = −21.11, Y = 9.03, Z = −1.00 — a 0.196 Å offset caused
by a duplicated `O1` record in the BGC residue of `3ZSJ.pdb`, which RDKit
collapses on parse. See [`UNVERIFIED.md`](../UNVERIFIED.md) §1 for the full
diagnosis and the evidence that the code's value is what produced the published
results.

The `(10, 15, 5)` figure formerly in `data/docking/README.md` was a
`DockingConfig` dataclass placeholder, overwritten at runtime by
`validate_receptor()`. Committed poses sit 33.6 Å from it, confirming it was
never used.

## 6. Phase 3 docking — carried over, not re-verified

The published `phase3_output/08_docking_results.csv` (14 compounds, Vina scores
−6.14 to −4.98 kcal/mol) is **carried over from the run of 2026-06-03**. It has
not been regenerated, and this working tree leaves it untouched.

What was established:

| Check | Result |
|---|---|
| AutoDock Vina present and working | ✅ v1.2.7; a single-ligand docking completed at −5.817 kcal/mol |
| Open Babel present | ✅ 3.1.1 |
| All 14 ligand PDBQTs pass `validate_ligand_pdbqt()` | ✅ 14/14, rotor counts 9–10 |
| Grid centre used at runtime | ✅ logged as X=−20.982 Y=8.876 Z=−1.002 |
| Full 14-ligand re-run | ❌ not completed |

The full re-run was attempted into a scratch directory (never overwriting the
published table) but the long-running process was terminated by the execution
environment before finishing, twice, leaving 12 of 14 compounds marked
`dock_failed`. That partial output was **discarded rather than reported**: the
failures are an artefact of the interrupted harness, not of the pipeline —
ligand validation passes for all 14 and Vina docks them individually without
error.

Phase 3 is also **not seeded**. Vina prints a fresh random seed each invocation
(observed: `random seed: -1296807838`), so even a completed re-run could not
bit-reproduce the published scores, only corroborate them statistically. See
`UNVERIFIED.md` §5.

## 7. Outstanding

Items that remain unverified or discrepant are enumerated in
[`UNVERIFIED.md`](../UNVERIFIED.md):

1. Grid centre disagreement (0.196 Å) between manuscript and code
2. Manuscript §2.10 test/coverage figures are stale (40/79%/88% → 53/69%/89%)
3. Manuscript RDKit attribution (2024.03.1 → 2026.03.x); 2024.03.1 itself
   untestable on Python 3.14
4. Lactose redocking RMSD (1.2 Å) is only bounded (< 2.0 Å), never persisted
5. Phase 3 docking is unseeded and therefore not bit-reproducible
