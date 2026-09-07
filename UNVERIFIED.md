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
only statistically.

`DockingConfig.seed` now exists for this: leave it at its `None` default and the
Vina command is byte-identical to the one that produced the published table; set
it to an integer and `--seed` is passed, making Phase 3 deterministic. The
default is deliberately unseeded so that the published numbers remain the output
of the pipeline's default configuration. The seed actually used is recorded in
`phase3_docking.log` on the `[seed]` line.

Verified on 2026-08-20 against Vina 1.2.7 using
`SCAF-001_BB-001-Ph_CuAAC_1.pdbqt` at the section 6 grid centre, exhaustiveness
8, two independent runs per condition:

| Condition | Run A (kcal/mol) | Run B (kcal/mol) | Identical |
|---|---|---|---|
| `seed=42` | -5.602, -5.554, -5.292 | -5.602, -5.554, -5.292 | yes |
| unseeded (default) | -5.895, -5.457, -5.142 | -5.919, -5.507, -5.153 | no |

So seeding does deliver bit-reproducible scores, and the unseeded default does
drift run to run -- here by up to 0.024 kcal/mol on the top pose, consistent
with the 0.099 kcal/mol mean absolute difference in section 6.

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

## 7. Lactose redocking RMSD does NOT match the manuscript, and stochasticity does not explain it

A re-run on 2026-08-20 gave **RMSD = 1.737 Å** against the 1.2 Å reported in the
manuscript (§2.7, §3.4, Figure S3). Both are below the 2.0 Å threshold, so both
pass validation.

An earlier version of this section attributed the gap to §5 (unseeded Vina).
**That explanation is wrong.** Measuring the distribution directly, five
independent unseeded runs:

| Run | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| RMSD (Å) | 1.726 | 1.715 | 1.740 | 1.723 | 1.731 |

Mean **1.727 Å**, SD **0.009 Å**, range 1.715-1.740 Å. With `seed=42` two runs
both give 1.716 Å exactly.

The run-to-run spread is therefore about 0.009 Å, roughly fifty times too small
to account for a 0.5 Å discrepancy. The redock RMSD is in practice a stable
quantity, and **the manuscript's 1.2 Å is not reproduced by the current
pipeline.** The cause is unidentified. Candidates not yet excluded: a different
receptor or ligand preparation than the one committed here, a different atom
subset (compare the pyranose-only 0.48 Å in the verified table), or a value
carried over from an earlier configuration. This should be resolved before
submission, since it is a number a reviewer can check directly.

## 8. The lactose validation reports a superposition RMSD, not a placement RMSD

`validate_receptor()` picks its RMSD metric by atom count:

* counts equal -> `compute_rmsd(..., align=False)`, an in-place RMSD in the
  binding-site frame, which is the correct redocking-validation metric
* counts differ -> `align_and_rmsd()`, which calls `rdMolAlign.AlignMol` and so
  **superposes the docked pose onto the crystal pose before measuring**

In practice the second branch always runs: the crystal ligand has 23 heavy
atoms, while the Open Babel PDBQT pose carries 31 (lactose's 8 hydroxyl
hydrogens are retained as polar H). Verified on 2026-08-20.

Superposition removes the translation and rotation that redocking validation
exists to test, so the reported figure is systematically optimistic: a pose in
the wrong place with the right internal geometry still scores well and still
passes the 2.0 Å threshold. What the fallback branch changes is not how the
pose is parsed but what is being measured.

### What the placement metric actually gives

Measured 2026-08-20, clean receptor, `seed=42`, all three returned poses:

| Pose | Vina score | Aligned RMSD (reported) | In-place RMSD | Centroid offset |
|---|---|---|---|---|
| 1 (top) | -4.589 | 1.716 Å | **7.854 Å** | 5.38 Å |
| 2 | -4.471 | 1.921 Å | 8.335 Å | 6.36 Å |
| 3 | -4.440 | **1.319 Å** | 8.732 Å | 7.04 Å |

Centroid offset is mapping-free and corroborates the in-place figures
independently of any substructure match: the crystal ligand's own maximum
internal span is 10.20 Å, so a 5.4 Å centroid displacement moves the pose more
than half its own length off the crystallographic site.

Two consequences:

1. **The redock does not reproduce the crystallographic binding mode.** The top
   pose sits 7.85 Å away in place. It passes the 2.0 Å threshold only because
   the reported metric superposes the pose onto the crystal first, which
   discards exactly the displacement being tested.
2. **The reported metric is anti-correlated with placement here.** Pose 3 is the
   furthest from the site (7.04 Å centroid offset) yet returns the *lowest*
   aligned RMSD of the three. Ranking poses by this number selects against
   correct placement.

Against the contaminated receptor the same run gives aligned 1.536 Å with an
8.65 Å centroid offset, so this is not an artefact of receptor preparation.

### Probable origin of the manuscript's 1.2 Å

Pose 3 returns **1.319 Å** by the reported metric. A value near 1.2 Å is
therefore reachable as the aligned RMSD of a non-top pose, or as the minimum
across returned poses, rather than the top pose's. This is consistent with the
manuscript's figure but is **inference, not a reproduction** -- the exact
configuration that produced 1.2 Å has not been recovered. What is established
is that 1.2 Å cannot be a placement RMSD, since every pose measured is above
7.8 Å in place.

Nothing here has been changed in the pipeline, because correcting the metric
would move a published validation number (§2.7, §3.4, Figure S3) and that is the
author's decision. But the claim the number is used to support -- that docking
reproduces the crystallographic lactose pose, and that the grid and receptor
preparation are therefore validated -- is not supported by the current
measurement, and this should be resolved before submission.

---

## 9. AI assistance, and what has and has not been audited

### Scope

Claude (Anthropic), used through Claude Code, contributed to this repository
between 3 June and 20 August 2026. Models used were Claude Opus 5, Claude
Sonnet 4.6, and Claude Haiku 4.5.

Twenty-two of the thirty-five commits reachable from `master` carry a
`Co-Authored-By: Claude` trailer (22 of 56 across all branches). All twenty-two
sit on the `master` lineage. The original Phase 1 and Phase 2 development, on
the `main` branch of 30-31 May 2026 -- a separate lineage that is not an
ancestor of `master` -- carries none.

By files touched, that work concentrated on:

| Area | File-touches |
|---|---|
| `src/` | 32 |
| `docs/` | 14 |
| `figures/` | 11 |
| `phase2_output/`, `manuscript/`, this file | 7 each |
| `tests/` | 6 |
| `scripts/` | 4 |

The most frequently modified individual files were `UNVERIFIED.md` (7),
`src/taloside_pipeline/phase3_docking.py` (5), `.gitignore` (5),
`requirements.txt` (4), `tests/test_phase3_docking.py` (3),
`src/taloside_pipeline/glycolibrary_generator.py` (3), and
`scripts/validation/compute_57i_pyranose_rmsd.py` (3).

### What this means for a reader

Substantial portions of the AI-assisted code, documentation and analysis were
**not independently audited at the time they were written**. Commit messages and
in-repository reports written during that period describe intent and results,
but do not constitute independent verification by the author.

This is not a hypothetical concern. The two files most directly implicated in
the validation defect recorded in §7 and §8 --
`src/taloside_pipeline/phase3_docking.py` and
`scripts/validation/compute_57i_pyranose_rmsd.py` -- are both among the
most-modified in this set.

The defect was also not caught by the review documents in this repository.
`UNVERIFIED.md`, `docs/audit_report.md`, `docs/final_validation_report.md` and
`manuscript/REVISION_NOTES.md` were each **created in AI-assisted commits dated
20 August 2026**. They are AI-generated review documents, not independent
verification by the author, and should be read as such -- including this
section.

### Status of the review

A review is underway. It is not complete, and no claim of comprehensive
verification is made here.

Verified to date by direct measurement, and recorded in §5-§8 of this document:
the seeding behaviour of Phase 3 docking; the run-to-run variance of the lactose
validation RMSD; the in-place placement of the docked lactose pose; and the
distinction between the aligned and direct 57I pyranose RMSD values.

Not yet verified: the remainder of the Phase 3 module, the figure-generation
code, the reports under `docs/`, and the numerical claims in the manuscript that
are not already itemised in this file.
