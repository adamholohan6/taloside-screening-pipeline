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

## 6. Phase 3 re-run corroborates the published scores, but raw-Vina rank order is not stable

The full 14-compound Phase 3 docking was re-run on 2026-08-20 (Vina 1.2.7,
Open Babel 3.1.1, grid centre X=-20.98 Y=8.88 Z=-1.00, exhaustiveness 8) and
compared against the published `phase3_output/08_docking_results.csv` from
2026-06-03.

| Statistic | Value |
|---|---|
| Compounds docked | 14 / 14 |
| Mean absolute difference | **0.099 kcal/mol** |
| Maximum absolute difference | 0.326 kcal/mol |
| Pearson correlation | **0.908** |
| Published range | -6.138 to -4.976 (span 1.162) |
| Re-run range | -6.329 to -4.968 (span 1.361) |

**Combined-score ranking is stable.** The top five by combined score are
identical in both set *and* order, and rank 1
(`SCAF-001_BB-004-4F_CuAAC_1`) is unchanged. This is the ranking the manuscript's
conclusions rest on.

**Raw-Vina ranking is not stable.** The top five by raw Vina score differ
between runs. `SCAF-001_BB-007-Pyridine_CuAAC_1` moves from 2nd (-6.008) to 7th
(-5.682, a +0.326 shift) while `SCAF-001_BB-001-Ph_CuAAC_1` rises into the top
five (-5.728 -> -5.965). Any statement that depends on raw-Vina rank order
should be read as one sample, not a reproducible ordering.

Both observations follow from §5: Vina is unseeded. The published values are
corroborated to within ~0.1 kcal/mol on average, which is the right standard for
an unseeded search, but they are not bit-reproducible and the published table has
deliberately been left in place rather than replaced with a second sample.

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
