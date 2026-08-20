# Manuscript revision notes — v10 / SI v7

Every change made to `Taloside_Manuscript_Final.docx` in producing
`Taloside_Manuscript_v10.docx` and `Taloside_SI_v7.docx`, with before/after text
so the revision can be checked without diffing binaries.

**Date:** 2026-08-20
**Source:** `manuscript/Taloside_Manuscript_Final.docx` (210 paragraphs, 8 tables, 11 images)
**Outputs:**
- `manuscript/Taloside_Manuscript_v10.docx` — 123 paragraphs, 2 tables, 7 images, 13 pages
- `manuscript/Taloside_SI_v7.docx` — 88 paragraphs, 6 tables, 4 images, 9 pages

Scope was deliberately limited to factual corrections. The argumentative
structure, the limitations section and the scoring-sensitivity analysis are
untouched; no conclusion was strengthened, softened, or restructured.

---

## 1. File split

The source document contained the entire SI appended after §8 References,
restarting at "Supporting Information" / "# Contents". The two halves were split
at that boundary.

The appended SI was used as the basis for `Taloside_SI_v7.docx` rather than the
existing standalone `Taloside_SI_v6_Final.docx`, because **the appended copy is
newer**: it contains Figure S4 (Binding Pocket Pose Overlay), which v6 lacks.

| | Source | Main v10 | SI v7 |
|---|---:|---:|---:|
| Paragraphs | 210 | 123 | 88 |
| Tables | 8 | 2 | 6 |
| Images | 11 | 7 | 4 |

Counts reconcile exactly (123 + 88 = 211 including the split paragraph;
7 + 4 = 11 images; 2 + 6 = 8 tables).

### Figure numbering after the split

- **Main:** Figures 1, 2, 3, 4, 5 — complete and in order. §3.4 also carries a
  deliberate cross-reference labelled "Figure S4 (shown here for context)",
  retained from the source.
- **SI:** Figures S1, S2, S3, S4 — complete and in order, each appearing once as
  a section heading and once as a caption (the source's existing pattern).

Orphaned media left behind by the split were purged: the main file dropped from
3447 KB to 2638 KB and the SI from 3444 KB to 1645 KB, with all 7 and 4 images
respectively still rendering.

---

## 2. Text changes

### 2.1 §2.2 — RuAAC reaction SMARTS (the defect you flagged)

The printed RuAAC SMARTS was missing the negative charge on the azide terminal
nitrogen. Confirmed by character-by-character comparison against
`src/taloside_pipeline/glycolibrary_generator.py`; divergence begins at index 15.

**Before**
```
RuAAC:[N:1]=[N+:2]=[N:3].[C:4]#[C:5]>>[C:5]1=[C:4][N:1][N+0:2]=[N+0:3]1
```
**After**
```
RuAAC:[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4][N:1][N+0:2]=[N+0:3]1
```

**The CuAAC SMARTS was also diffed and is byte-identical to the source — no
change made.**

### 2.2 RDKit version attribution (three locations, not one)

The source attributed the charge-propagation failure mode to RDKit 2024.03.1.
This was tested rather than assumed: running the pre-fix template under the
installed RDKit 2026.03.2 reproduces the documented error exactly
(`AtomValenceException: Explicit valence for atom # 4 N, 3, is greater than
permitted`), and the CHANGELOG and code comments both say 2026.03.x.

RDKit 2024.03.1 could not be tested directly — it publishes no distribution for
Python 3.14 — so the corrected text states the version the results *were*
produced under, which is verified.

| Location | Before | After |
|---|---|---|
| §2.1 scaffold SMILES | `Validated scaffold SMILES (RDKit 2024.03.1):` | `… (RDKit 2026.03.2):` |
| §2.2 failure mode | `a failure mode in RDKit 2024.03.1.` | `a failure mode in RDKit 2026.03.2.` |
| §2.3 descriptors | `calculated using RDKit 2024.03.1 after MolStandardize` | `calculated using RDKit 2026.03.2 …` |

### 2.3 §2.7 — docking grid centre

**Before:** `Grid: 20 × 20 × 20 Å, centre X = −21.11, Y = 9.03, Z = −1.00 Å.`
**After:** `Grid: 20 × 20 × 20 Å, centre X = −20.98, Y = 8.88, Z = −1.00 Å.`

The published value was the centroid of all 24 BGC + GAL records in `3ZSJ.pdb`.
The BGC residue contains **two records both named `O1`**, and
`Chem.MolFromPDBBlock()` collapses the duplicate, so the pipeline actually
averages 23 atoms. The pipeline's own log for the reported run records:

```
Dynamically centered grid on crystal ligand: X=-20.98, Y=8.88, Z=-1.00
[grid] Docking grid centre: X=-20.982  Y=8.876  Z=-1.002  box=20.0 A^3
```

Evidence that this is the value that produced the published results: dynamic
centring was introduced in commit `6393ce0` (2026-06-03), the outputs were
regenerated in `91cc7db` the same day, `extract_ligand_coords_from_pdb` has not
changed since, and the committed poses lie 33.6 Å from the `(10, 15, 5)`
`DockingConfig` default versus 3.1 Å from the dynamic centre.

The offset is 0.196 Å inside a 20 Å box and is very unlikely to affect any
docking result, but the manuscript previously stated a centre the code does not
produce. Full diagnosis in `UNVERIFIED.md` §1.

### 2.4 §2.10 — software validation figures

**Before**
> The pipeline includes 40 pytest tests (79% overall coverage; 88% of
> phase2_integration.py).

**After**
> The pipeline includes 53 pytest tests (69% overall coverage; 89% of
> phase2_integration.py).

From an actual run: `pytest --cov=src --cov-report=term-missing` → 53 passed,
TOTAL 1023 statements / 314 missed / 69%, `phase2_integration.py` 211/24/89%.
Breakdown in `docs/coverage_report.md`.

### 2.5 §5 Conclusions — same figures, second location

Not listed in the brief, but the same stale numbers appear again here.

**Before:** `40 automated tests (79% coverage)`
**After:** `53 automated tests (69% coverage)`

### 2.6 Zenodo DOI (three locations)

The manuscript cited `10.5281/zenodo.20529842`, which is the **version-specific
DOI for v2.0.0** — a superseded release. Confirmed against the Zenodo API that
`10.5281/zenodo.20476421` is the concept DOI (`conceptrecid` 20476421), which
always resolves to the latest version.

| Location | Before | After |
|---|---|---|
| Abstract | `Code archive: …zenodo.20529842 (v2.0.0)` | `…zenodo.20476421 (concept DOI; this work archived as v2.2.0)` |
| §2.10 | `Zenodo: 10.5281/zenodo.20529842 (v2.0.0)` | `Zenodo: 10.5281/zenodo.20476421 (concept DOI; this work archived as v2.2.0)` |
| §6 | `Zenodo (v2.0.0): https://doi.org/10.5281/zenodo.20529842` | `Zenodo (concept DOI, always resolves to the latest release): https://doi.org/10.5281/zenodo.20476421` |

The DOIs were embedded as **hyperlinks**, so both the visible text and the
underlying relationship targets were updated. Verified: the only external Zenodo
target remaining in either file is `https://doi.org/10.5281/zenodo.20476421`.
Changing the display text alone would have left the links pointing at the
superseded record.

### 2.7 §6 — Code Availability release pointer (new paragraph)

Appended after the existing §6 entries:

> The exact code state underlying the results reported here is archived as
> release v2.2.0:
> https://github.com/adamholohan6/taloside-screening-pipeline/releases/tag/v2.2.0.
> Results were produced with RDKit 2026.03.2, AutoDock Vina 1.2.7 and Open Babel
> 3.1.1 under Python 3.14. Claims that could not be reproduced from the archive
> are listed in UNVERIFIED.md within it.

### 2.8 SI Note S1 — grid centre, and a missing value

The SI repeated the grid centre, and its z-coordinate was **absent entirely**.
This was not in the brief.

**Before**
```
Centre: x = -21.11 Å,  y = 9.03 Å,  z =  Å  (dynamically centred on BGC+GAL
crystal ligand centroid, 3ZSJ).
```
**After**
```
Centre: x = -20.98 Å,  y = 8.88 Å,  z = -1.00 Å  (dynamically centred on
BGC+GAL crystal ligand centroid, 3ZSJ).
```

The parenthetical was already correct about the mechanism; only the numbers were
wrong, and z was blank.

---

## 3. Deliberately NOT changed

| Claim | Why left alone |
|---|---|
| Lactose redocking RMSD = 1.2 Å (§2.7, §3.4, Figure S3) | **Not reproduced.** Five independent runs give 1.727 ± 0.009 Å; the value is stable, so run-to-run variation does not explain the 0.5 Å gap. Cause unidentified. Separately, the metric itself is a superposition RMSD rather than an in-place one, which makes it optimistic. Both pass the 2.0 Å threshold. Needs resolving before submission. See `UNVERIFIED.md` §7 and §8. |
| 57I pyranose RMSD = 0.48 Å (§2.7, §3.4) | Re-verified as correct — `scripts/validation/compute_57i_pyranose_rmsd.py` exits 0 reproducing 0.48 Å. No change needed. |
| Vina scores, combined scores, Tables 1/2/S2 rankings | Phase 3 has since been re-run in full, 14/14 (see `UNVERIFIED.md` §6): mean \|Δ\| 0.099 kcal/mol, r = 0.908, combined-score top five identical in set and order. Carried over unchanged rather than replaced with a second unseeded sample. Raw-Vina rank order is **not** stable between runs. |
| Cα RMSD 0.275 Å, exhaustiveness sensitivity (Table S4) | Cα RMSD re-verified as correct. Table S4 not re-run. |
| AutoDock Vina 1.2.7, Open Babel 3.1.1, AutoDockTools 1.5.7 | Vina and Open Babel confirmed by `--version`. AutoDockTools not independently checked. |
| Limitations section, scoring-sensitivity analysis, all argumentation | Out of scope by instruction. |

---

## 4. Verification performed

All 15 assertions below were checked programmatically against the output files:

- RuAAC SMARTS contains `[N-:3]`
- No occurrence of `2024.03.1` remains; `RDKit 2026.03.2` present
- Grid centre `−20.98` / `8.88` present; no `21.11` remains (main **and** SI)
- `53 pytest tests`, `69% overall coverage`, `89% of phase2_integration` present
- No `40 pytest tests` / `40 automated tests` remains
- `53 automated tests (69% coverage)` present in Conclusions
- Concept DOI `20476421` present; no `20529842` remains
- Only external Zenodo relationship target is the concept DOI
- v2.2.0 release pointer present
- SI grid centre reads `x = -20.98 … y = 8.88 … z = -1.00`

Both files were converted to PDF with LibreOffice to confirm they open and
paginate cleanly: main 13 pages / 7 images, SI 9 pages / 4 images.
