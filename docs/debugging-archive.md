# Debugging scripts archive (June 2026 bug hunt)

Eighteen one-off investigation scripts were removed from the repository root in
the v2.2.0 cleanup. They were never part of the pipeline and were not collected
by pytest (`testpaths = tests`), but several were named `test_*.py`, which made
them look like part of the test suite.

They remain in git history. The last commit containing all of them is:

    b88caf073255aef4057e68bfbf2501ada0674cae

Recover any one of them with:

```bash
git show b88caf073255aef4057e68bfbf2501ada0674cae:<filename>
```

For example:

```bash
git show b88caf073255aef4057e68bfbf2501ada0674cae:rdkit_charge_investigation.py
git show b88caf073255aef4057e68bfbf2501ada0674cae:test_triazole_smarts.py
```

## Files removed

| File | Purpose |
|------|---------|
| `analyze_phase2_smarts.py` | SMARTS pattern inspection during the triazole bug hunt |
| `chemistry_investigation.py` | Triazole valence/aromaticity exploration |
| `diagnostic_evidence.py` | Evidence collection for the root-cause report |
| `investigate_scaffold.py` | Scaffold enumeration checks |
| `mapping_investigation.py` | Atom-map propagation through RunReactants |
| `proof_of_concept_unmapped.py` | Unmapped-template proof of concept |
| `rdkit_charge_investigation.py` | Formal-charge propagation into product templates |
| `validate_phase2_final.py` | Ad-hoc Phase 2 validation |
| `validate_phase2_products.py` | Ad-hoc product validation |
| `apply_patch.py` | One-off patch applier |
| `fix_backslash_newlines.py` | One-off text fixer |
| `test_phase2_diagnostics.py` | Ad-hoc Phase 2 diagnostics (not a pytest test) |
| `test_phase2_fix.py` | Ad-hoc fix verification (not a pytest test) |
| `test_phase3_fix.py` | Ad-hoc fix verification (not a pytest test) |
| `test_single_smiles.py` | Single-SMILES probe (not a pytest test) |
| `test_smiles_simple.py` | SMILES parsing probe (not a pytest test) |
| `test_triazole.py` | Triazole probe (not a pytest test) |
| `test_triazole_smarts.py` | Triazole SMARTS probe (not a pytest test) |

The maintained test suite lives in `tests/`.
