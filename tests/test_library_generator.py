"""Tests for the taloside library generator."""

from pathlib import Path

from rdkit import Chem

from taloside_pipeline.glycolibrary_generator import (
    GlycoLibraryGenerator,
    LibraryConfig,
    ReactionSMARTS,
    generate_triazole_library,
)


SCAFFOLD = (
    "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)"
    "C4=C([N+]([O-])=O)C=CC=C4"
)
BUILDING_BLOCKS = [
    {'id': 'BB-001-Ph', 'smiles': 'C#Cc1ccccc1'},
    {'id': 'BB-002-4OMe', 'smiles': 'C#Cc1ccc(OC)cc1'},
]


def test_library_config_defaults():
    config = LibraryConfig()
    assert config.max_products == 50000
    assert config.include_stereoisomers is True
    assert config.sanitize_products is True
    assert config.filter_hypervalent is True
    assert config.min_product_mw == 150.0
    assert config.max_product_mw == 1000.0


def test_compatibility_wrapper_imports():
    # Import through the legacy module path to ensure the wrapper works.
    from taloside_pipeline.library_generator import LibraryConfig as LegacyConfig
    assert LegacyConfig is LibraryConfig


def test_triazole_regioisomer_constants_present():
    assert ReactionSMARTS.TRIAZOLE_1_4_CuAAC["regioisomer"] == "1,4-CuAAC"
    assert ReactionSMARTS.TRIAZOLE_1_5_RuAAC["regioisomer"] == "1,5-RuAAC"


def test_generate_triazole_library_returns_expected_columns(tmp_path):
    df = generate_triazole_library(
        scaffold_smiles=SCAFFOLD,
        building_blocks=BUILDING_BLOCKS,
        config=LibraryConfig(max_products=10, output_dir=tmp_path),
    )
    assert not df.empty
    expected = {
        'compound_id', 'parent_scaffold_id', 'building_block_id', 'regioisomer',
        'product_smiles', 'product_inchi', 'product_inchikey', 'molecular_weight',
        'h_donors', 'h_acceptors', 'logp', 'tpsa', 'rotatable_bonds', 'generation_status',
    }
    assert expected.issubset(df.columns)
    assert set(df['regioisomer'].unique()) <= {'1,4-CuAAC', '1,5-RuAAC'}


def test_generator_initializes_with_valid_inputs():
    generator = GlycoLibraryGenerator(
        scaffold_smiles=SCAFFOLD,
        building_blocks=BUILDING_BLOCKS,
        reaction_smarts=ReactionSMARTS.TRIAZOLE_1_4_CuAAC['smarts'],
        regioisomer_label='1,4-CuAAC',
        config=LibraryConfig(max_products=5),
    )
    assert generator.scaffold_mol is not None
    assert len(generator.building_blocks) == 2
