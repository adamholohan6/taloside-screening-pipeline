"""Integration and unit tests for phase2_integration module."""

import pandas as pd
import pytest
from rdkit import Chem

from taloside_pipeline.glycolibrary_generator import generate_triazole_library
from taloside_pipeline.phase2_integration import (
    PAINSFilter,
    apply_lipinski_filter,
    apply_pains_filter,
    compute_lead_scores,
    run_phase2_pipeline,
)

from tests.constants import LIBRARY_COLUMNS, PHASE2_EXPORT_FILES, SCAFFOLD


@pytest.fixture
def sample_library_df(building_blocks_phase2, phase2_config):
    df = generate_triazole_library(
        scaffold_smiles=SCAFFOLD,
        building_blocks=building_blocks_phase2,
        config=phase2_config,
    )
    assert not df.empty
    return df


@pytest.mark.unit
def test_apply_lipinski_filter_strict_mode(sample_library_df):
    passed, failed = apply_lipinski_filter(sample_library_df, strict_mode=True)
    total = len(passed) + len(failed)
    assert total == len(sample_library_df)
    assert len(passed) == 14
    assert len(failed) == 2
    assert (passed["molecular_weight"] <= 600).all()
    assert (passed["h_acceptors"] <= 12).all()


@pytest.mark.unit
def test_apply_lipinski_filter_classical_mode(sample_library_df):
    passed, failed = apply_lipinski_filter(sample_library_df, strict_mode=False)
    strict_passed, _ = apply_lipinski_filter(sample_library_df, strict_mode=True)
    assert len(passed) <= len(strict_passed)


@pytest.mark.unit
def test_compute_lead_scores_formula_and_sorting(sample_library_df):
    passed, _ = apply_lipinski_filter(sample_library_df, strict_mode=True)
    scored = compute_lead_scores(passed)

    assert "lead_score" in scored.columns
    assert scored["lead_score"].between(0, 1).all()
    assert scored["lead_score"].is_monotonic_decreasing
    assert list(scored.columns) == list(passed.columns) + ["lead_score"]


@pytest.mark.unit
def test_compute_lead_scores_single_compound():
    df = pd.DataFrame(
        [{
            "compound_id": "TEST-1",
            "tpsa": 100.0,
            "molecular_weight": 500.0,
            "logp": 1.0,
            "rotatable_bonds": 5,
        }]
    )
    scored = compute_lead_scores(df)
    assert len(scored) == 1
    assert scored.iloc[0]["lead_score"] == pytest.approx(1.0)


@pytest.mark.unit
def test_pains_filter_clean_molecule():
    mol = Chem.MolFromSmiles("CCO")
    is_pains, matches = PAINSFilter.screen_molecule(mol)
    assert is_pains is False
    assert matches == []


@pytest.mark.unit
def test_pains_filter_undetermined_when_catalog_unavailable(monkeypatch):
    monkeypatch.setattr(PAINSFilter, "_catalog_loaded", True)
    monkeypatch.setattr(PAINSFilter, "_catalog", None)

    mol = Chem.MolFromSmiles("CCO")
    is_pains, matches = PAINSFilter.screen_molecule(mol)
    assert is_pains is None
    assert matches == []


@pytest.mark.unit
def test_apply_pains_filter_routes_invalid_smiles_to_undetermined():
    df = pd.DataFrame([{"compound_id": "BAD-1", "product_smiles": "INVALID"}])
    clean, pains, undetermined = apply_pains_filter(df)

    assert len(clean) == 0
    assert len(pains) == 0
    assert len(undetermined) == 1
    assert undetermined.iloc[0]["pains_status"] == "UNDETERMINED - invalid SMILES"


@pytest.mark.unit
def test_apply_pains_filter_routes_catalog_unavailable(monkeypatch, sample_library_df):
    passed, _ = apply_lipinski_filter(sample_library_df, strict_mode=True)
    monkeypatch.setattr(PAINSFilter, "_catalog_loaded", True)
    monkeypatch.setattr(PAINSFilter, "_catalog", None)

    clean, pains, undetermined = apply_pains_filter(passed.head(2))
    assert len(clean) == 0
    assert len(pains) == 0
    assert len(undetermined) == 2
    assert (undetermined["pains_status"] == "UNDETERMINED - catalog unavailable").all()


@pytest.mark.integration
def test_apply_pains_filter_clean_library(sample_library_df):
    passed, _ = apply_lipinski_filter(sample_library_df, strict_mode=True)
    clean, pains, undetermined = apply_pains_filter(passed)

    assert len(clean) + len(pains) + len(undetermined) == len(passed)
    assert len(clean) == 14
    assert len(pains) == 0
    assert len(undetermined) == 0
    assert (clean["pains_status"] == "CLEAN").all()


@pytest.mark.integration
def test_full_phase2_workflow_components(sample_library_df, phase2_config):
    library_df = sample_library_df
    assert LIBRARY_COLUMNS.issubset(library_df.columns)
    assert len(library_df) == 16
    assert set(library_df["regioisomer"].unique()) == {"1,4-CuAAC", "1,5-RuAAC"}

    lipinski_passed, lipinski_failed = apply_lipinski_filter(library_df, strict_mode=True)
    assert len(lipinski_passed) == 14
    assert len(lipinski_failed) == 2

    clean_df, pains_df, undetermined_df = apply_pains_filter(lipinski_passed)
    assert len(clean_df) == 14

    scored_df = compute_lead_scores(lipinski_passed)
    assert len(scored_df) == 14
    assert scored_df.iloc[0]["lead_score"] >= scored_df.iloc[-1]["lead_score"]

    phase2_config.output_dir.mkdir(parents=True, exist_ok=True)
    exports = [
        (library_df, "01_all_generated_compounds.csv"),
        (lipinski_passed, "02_lipinski_passed.csv"),
        (lipinski_failed, "03_lipinski_failed.csv"),
        (clean_df, "04_lipinski_clean_no_pains.csv"),
        (pains_df, "05_lipinski_with_pains.csv"),
        (undetermined_df, "06_pains_undetermined.csv"),
        (scored_df, "07_lead_scored.csv"),
    ]
    for df_export, fname in exports:
        path = phase2_config.output_dir / fname
        df_export.to_csv(path, index=False)
        assert path.exists()


@pytest.mark.integration
@pytest.mark.slow
def test_run_phase2_pipeline_smoke(tmp_path, monkeypatch):
    """End-to-end smoke test for run_phase2_pipeline() CSV export."""
    monkeypatch.chdir(tmp_path)

    PAINSFilter._catalog_loaded = False
    PAINSFilter._catalog = None

    result = run_phase2_pipeline()
    assert result is not None

    output_dir = tmp_path / "phase2_output"
    for fname in PHASE2_EXPORT_FILES:
        assert (output_dir / fname).exists(), f"Missing export: {fname}"

    library_df = pd.read_csv(output_dir / "01_all_generated_compounds.csv")
    lipinski_passed = pd.read_csv(output_dir / "02_lipinski_passed.csv")
    lipinski_failed = pd.read_csv(output_dir / "03_lipinski_failed.csv")
    clean_df = pd.read_csv(output_dir / "04_lipinski_clean_no_pains.csv")
    scored_df = pd.read_csv(output_dir / "07_lead_scored.csv")

    assert len(library_df) == 16
    assert len(lipinski_passed) == 14
    assert len(lipinski_failed) == 2
    assert len(clean_df) == 14
    assert len(scored_df) == 14
    assert pytest.approx(len(lipinski_passed) / len(library_df), rel=0.01) == 0.875

    assert set(library_df["regioisomer"].unique()) == {"1,4-CuAAC", "1,5-RuAAC"}
    assert "lead_score" in scored_df.columns
    assert LIBRARY_COLUMNS.issubset(set(library_df.columns))

    assert lipinski_failed["building_block_id"].str.contains("2NO2").any()
