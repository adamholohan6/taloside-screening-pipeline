"""Shared fixtures for taloside pipeline tests."""

from pathlib import Path
import sys

import pytest

from tests.constants import (
    BUILDING_BLOCKS_MINI,
    BUILDING_BLOCKS_PHASE2,
    SCAFFOLD,
)

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@pytest.fixture
def scaffold():
    return SCAFFOLD


@pytest.fixture
def building_blocks_mini():
    return BUILDING_BLOCKS_MINI.copy()


@pytest.fixture
def building_blocks_phase2():
    return BUILDING_BLOCKS_PHASE2.copy()


@pytest.fixture
def phase2_config(tmp_path):
    from taloside_pipeline.glycolibrary_generator import LibraryConfig

    return LibraryConfig(
        max_products=500,
        min_product_mw=250.0,
        max_product_mw=800.0,
        include_stereoisomers=True,
        filter_hypervalent=True,
        output_dir=tmp_path / "phase2_output",
    )
