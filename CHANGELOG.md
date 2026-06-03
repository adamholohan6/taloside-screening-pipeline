# CHANGELOG

## Pipeline Correction and Revalidation (June 2026)

### What Was Broken

When the corrected two-template pipeline (separate CuAAC and RuAAC reaction
SMARTS) was first run under RDKit 2026.03.x, it produced **zero compounds**.
The symptom — every product silently dropped — masked three independent bugs
that together made the outputs unreliable even when earlier runs appeared to
succeed.

### Bug 1: Reactant Charges Propagating Into Products (SMARTS)

The triazole product templates used Kekulé uppercase atoms with no explicit
charge on ring N2/N3. Under RDKit 2026.03.x, formal charges from the azide
reactant (`[N+]`, `[N-]`) propagate into unspecified product atoms. The
inherited `[N-]` on a ring nitrogen produced a bond-order sum of 3, exceeding
RDKit's permitted valence of 2, so `SanitizeMol` raised an explicit-valence
error and every product was dropped. The fix was to pin N2 and N3 explicitly
neutral with `[N+0]` in both the CuAAC and RuAAC product templates. This is a
required defence, not a style choice: without it, neutrality is aspirational
rather than enforced.

### Bug 2: Both Alkyne Orientations Accepted (Geometry Filter)

The SMARTS pattern `[C:4]#[C:5]` matches a terminal alkyne in both
orientations. For each building block, each template therefore produced two
products — one with the correct triazole topology and one with the mirror
regiochemistry. The wrong-orientation product from the CuAAC run was
identical (same InChIKey) to the correct product from the RuAAC run, and vice
versa. Rather than catching this as a labelling error, the pipeline was
silently keeping both under the wrong label. A post-generation geometry filter
using two SMARTS patterns (`[#6][n]1[n][n][c]([c])[c]1` for 1,4;
`[#6][n]1[n][n][c][c]1[c]` for 1,5) now discards any product whose ring
topology does not match its template label. The result is 16 unique, correctly
labelled compounds — not 32 duplicates.

### Bug 3: Compound IDs Not Unique Across Regioisomers (Many-to-Many Merge)

`compound_id` was constructed from scaffold ID, building block ID, and a
per-building-block counter. Because the counter reset between the CuAAC and
RuAAC generator runs, both produced `SCAF-001_BB-001-Ph_1`. When the
lead-score merge joined on `compound_id`, each ID appeared twice on both
sides: 14 unique compounds became 56 rows in `07_lead_scored.csv`. The fix
embeds the regioisomer tag in the ID (`_CuAAC_1` / `_RuAAC_1`), making it
globally unique. Comments at both merge call sites document the assumption
and the composite-key upgrade path if it is ever violated again.

### Data Regeneration

All Phase 2 CSVs (`01` through `07`) were regenerated from the corrected
pipeline. The stale `08_docking_results.csv` — produced by the pre-fix
pipeline with old-format IDs — was deleted. Phase 3 was re-run against the
corrected `07_lead_scored.csv`, producing a new `08_docking_results.csv` in
which all 14 lead compounds dock successfully (Vina scores −4.98 to
−6.14 kcal/mol). The four docking failures in the old file (both regioisomers
of the 4-chlorophenyl and 3-bromophenyl series) are gone: the geometry-
corrected SMILES embed reliably into 3D conformers.

### Testing

`infer_triazole_regioisomer` — a post-hoc SMARTS inference function that
fired a mismatch warning on every compound but whose findings were always
overridden — was removed. It was detecting a real labelling problem (the
geometry bug above) but was not acting on it, making it dead code that looked
functional. Its stale test, which used non-standard reference SMILES, was
removed from `test_smarts_validation.py`. All 11 remaining SMARTS tests pass.
An Open Babel pre-flight check that used `--version` (exits non-zero on
OpenBabel 3.1.x) was corrected to `-V`.

### Current State

End-to-end: 16 geometrically verified compounds generated → 14 pass
carbohydrate-adjusted Lipinski filters → 14 PAINS-clean → 14 lead-scored →
**14/14 docked successfully** against the 3ZSJ Gal-3 CRD. Phase 2 and Phase 3
outputs now come from a single consistent pipeline version and carry matching
compound IDs throughout.
