"""Smoke tests for package entry points and public API."""

import subprocess
import sys

import pytest

from tests.constants import ROOT, SCAFFOLD, BUILDING_BLOCKS_PHASE2, LIBRARY_COLUMNS, PHASE2_EXPORT_FILES


@pytest.mark.unit
def test_package_import_and_version():
    import taloside_pipeline

    assert taloside_pipeline.__version__ == "1.0.0"
    assert hasattr(taloside_pipeline, "run_phase2_pipeline")
    assert hasattr(taloside_pipeline, "generate_triazole_library")


@pytest.mark.unit
def test_public_api_exports():
    import taloside_pipeline

    expected = {
        "calculate_descriptors",
        "GlycoLibraryGenerator",
        "generate_triazole_library",
        "run_phase2_pipeline",
        "Phase2Integration",
    }
    assert expected.issubset(set(taloside_pipeline.__all__))


@pytest.mark.unit
def test_library_generator_compatibility_wrapper():
    from taloside_pipeline.library_generator import GlycoLibraryGenerator, ReactionSMARTS

    assert GlycoLibraryGenerator is not None
    assert ReactionSMARTS.TRIAZOLE_1_4_CuAAC["regioisomer"] == "1,4-CuAAC"


@pytest.mark.smoke
def test_module_main_phase2_integration_runnable():
    """python -m taloside_pipeline.phase2_integration must be invokable."""
    result = subprocess.run(
        [sys.executable, "-m", "taloside_pipeline.phase2_integration"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    # Allow success even if Windows console logging emits Unicode warnings
    assert result.returncode == 0
    for fname in PHASE2_EXPORT_FILES:
        assert (ROOT / "phase2_output" / fname).exists()


@pytest.mark.smoke
def test_descriptor_calculator_entry_point(tmp_path):
    from taloside_pipeline.descriptor_calculator import calculate_library_descriptors

    out = tmp_path / "smoke_descriptors.csv"
    df, failed = calculate_library_descriptors(out)
    assert failed == 0
    assert len(df) == 10
    assert out.exists()


@pytest.mark.smoke
def test_triazole_library_quick_smoke():
    from taloside_pipeline.glycolibrary_generator import LibraryConfig, generate_triazole_library

    df = generate_triazole_library(
        scaffold_smiles=SCAFFOLD,
        building_blocks=BUILDING_BLOCKS_PHASE2,
        config=LibraryConfig(max_products=500),
    )
    assert len(df) == 16
    assert LIBRARY_COLUMNS.issubset(df.columns)
