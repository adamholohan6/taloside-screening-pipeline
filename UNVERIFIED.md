# Unverified and discrepant claims

Claims that could not be reproduced from files in this repository, or that
reproduce to a *different* value than the one currently published. Nothing here
has been silently corrected; each entry records the original claim, what was
actually computed, and what is needed to settle it.

Last updated: 2026-08-20 (v2.2.0 cleanup)

---

## 1. Grid centre: manuscript disagrees with the code by 0.196 Å

| Source | Value | Basis |
|---|---|---|
| Manuscript §2.7 | X = −21.11, Y = 9.03, Z = −1.00 | Centroid of all 24 BGC + GAL **PDB records** |
| Code at runtime | X = −20.98, Y = 8.88, Z = −1.00 | Centroid of the 23 **RDKit-parsed atoms** |
| `data/docking/README.md` (before this cleanup) | (10, 15, 5) | `DockingConfig` dataclass default — never used |

**Cause.** The BGC residue in `3ZSJ.pdb` contains two records both named `O1`.
`extract_ligand_coords_from_pdb()` builds a PDB block from all 24 records and
parses it with `Chem.MolFromPDBBlock(...)`, which collapses the duplicate name,
yielding 23 atoms. The manuscript figure counts the duplicated oxygen twice; the
code does not.

**Which was used for the published results.** The code's value. Dynamic grid
centring was introduced in `6393ce0` (2026-06-03) and the outputs were
regenerated in `91cc7db` the same day; `extract_ligand_coords_from_pdb` has not
changed since. Independent confirmation: the committed poses in
`phase3_output/poses/` have a mean centroid 3.1 Å from the dynamic centre and
33.6 Å from the `(10, 15, 5)` dataclass default, so the default was demonstrably
not in play.

**Status.** Unresolved in the manuscript. The discrepancy is 0.196 Å inside a
20 Å box, so it is very unlikely to change any docking result, but §2.7 states a
centre that the code does not produce. Reproduce with:

```bash
python -c "import sys; sys.path.insert(0,'src'); from pathlib import Path; \
from taloside_pipeline.phase3_docking import extract_ligand_coords_from_pdb; \
c,_=extract_ligand_coords_from_pdb(Path('data/docking/3ZSJ.pdb')); print(c.mean(axis=0))"
```

---

## 2. Manuscript §2.10 software-validation figures are stale

| Claim (§2.10) | Actual (2026-08-20) |
|---|---|
| 40 pytest tests | **53 tests, all passing** |
| 79% overall coverage | **69%** |
| 88% of `phase2_integration.py` | **89%** |

Reproduce with `pytest --cov=src --cov-report=term-missing`.
Corrected numbers are in `docs/coverage_report.md`.

---

## 3. Manuscript RDKit version attribution was wrong

§2.2 attributes the charge-propagation failure mode to **RDKit 2024.03.1**. The
CHANGELOG and the comment at `glycolibrary_generator.py:161` both say
**2026.03.x**, and the latter is correct: running the pre-fix template under the
installed RDKit 2026.03.2 reproduces the documented error exactly
(`AtomValenceException: Explicit valence for atom # 4 N, 3, is greater than
permitted`), while the current template sanitises cleanly.

RDKit 2024.03.1 could not be tested directly: it publishes no distribution for
Python 3.14 (oldest available is 2025.3.4), so whether that specific version
*also* exhibits the behaviour is **untested**. The claim corrected in the code
and requirements is "the reported results were produced under 2026.03.2", which
is verified.

---

## 4. Lactose redocking RMSD (1.2 Å) is validated only indirectly

The manuscript reports a lactose redocking RMSD of 1.2 Å. This is not written to
any output file. `VinaDocking.validate_receptor()` computes it during the run and
aborts if it exceeds `rmsd_threshold_angstrom` (2.0 Å), so a successful pipeline
run establishes only that the value was **below 2.0 Å**, not that it was 1.2 Å.

To verify the exact figure, `validate_receptor()` would need to persist the
computed RMSD to the run log or results table. It currently does not.

---

## 5. Phase 3 docking is not seeded

`DockingConfig` sets `exhaustiveness = 8` but no `--seed` is passed to Vina
(`phase3_docking.py`, the argument list around line 586). RDKit conformer
embedding *is* seeded (`embed_ligand_3d(..., random_seed=42)`), but Vina's Monte
Carlo search is not, so Vina scores are not bit-reproducible between runs.

This means the published `08_docking_results.csv` cannot be reproduced exactly,
only statistically. Adding `--seed` to the Vina invocation would make Phase 3
deterministic; this has not been changed here because it would alter the
published numbers.

## 6. Phase 3 numbers are carried over, not re-verified

The 14 Vina scores in `phase3_output/08_docking_results.csv` (−6.14 to −4.98
kcal/mol) date from the run of 2026-06-03 and were **not** regenerated during
the v2.2.0 cleanup. A full re-run into a scratch directory was attempted twice
and was terminated by the execution environment before completing.

The surrounding conditions were verified: Vina 1.2.7 and Open Babel 3.1.1 are
installed and working, a single ligand docks successfully (−5.817 kcal/mol),
and all 14 ligand PDBQTs pass `validate_ligand_pdbqt()`. Only the end-to-end
run is outstanding. Because of §5 above it could in any case only corroborate,
never bit-reproduce, the published values.

To attempt it:

```bash
# with vina and obabel on PATH
python -m taloside_pipeline.phase3_docking
```

---

## Claims that DID verify

Recorded here so they are not re-investigated:

| Claim | Result |
|---|---|
| 16 compounds generated | ✅ exact |
| 14 Lipinski pass (87.5%) | ✅ exact |
| 0 PAINS hits | ✅ exact |
| Lead scores 0.399–0.905 | ✅ exact |
| 57I pyranose RMSD 0.48 Å (Kabsch) | ✅ exact, script exits 0 |
| Cα RMSD 0.275 Å over 138 atoms | ✅ exact |
| AutoDock Vina 1.2.7 | ✅ `vina --version` |
| Open Babel 3.1.1 | ✅ `obabel -V` |
| Zenodo concept DOI 10.5281/zenodo.20476421 | ✅ confirmed via API |
| Committed `phase2_output/*.csv` | ✅ bit-identical to a clean rerun |

---

## 7. Lactose redocking RMSD varies between runs (measured)

A partial re-run on 2026-08-20 completed the lactose validation redock and gave
**RMSD = 1.737 Å**, against the 1.2 Å reported in the manuscript. Both are below
the 2.0 Å threshold, so both pass validation.

This is the unseeded behaviour of §5 made concrete: Vina draws a fresh random
seed per invocation, so the lactose redock lands on a different pose each time
and the validation RMSD moves with it. The manuscript's 1.2 Å is therefore one
sample, not a reproducible constant. Adding `--seed` to the Vina invocation
would fix this; it has not been changed here because it would alter the
published numbers.
