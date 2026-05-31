"""SMARTS validation tests for reaction templates."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from taloside_pipeline.glycolibrary_generator import (
    GlycoLibraryGenerator,
    LibraryConfig,
    ReactionSMARTS,
    infer_triazole_regioisomer,
)

from tests.constants import BUILDING_BLOCKS_MINI, SCAFFOLD


# All ReactionSMARTS dict entries that carry a 'smarts' key
TRIAZOLE_SMARTS = (
    ReactionSMARTS.TRIAZOLE_1_4_CuAAC,
    ReactionSMARTS.TRIAZOLE_1_5_RuAAC,
    ReactionSMARTS.TRIAZOLE_FORMATION,
)
OTHER_SMARTS = (
    ReactionSMARTS.AMIDE_COUPLING,
    ReactionSMARTS.SUZUKI_COUPLING,
)
ALL_REACTION_DEFS = TRIAZOLE_SMARTS + OTHER_SMARTS


@pytest.mark.unit
@pytest.mark.parametrize("rxn_def", ALL_REACTION_DEFS, ids=lambda d: d["name"])
def test_reaction_smarts_parses(rxn_def):
    """Each declared SMARTS string must parse into a valid RDKit reaction."""
    rxn = AllChem.ReactionFromSmarts(rxn_def["smarts"])
    assert rxn is not None
    assert rxn.GetNumReactantTemplates() >= 1
    assert rxn.GetNumProductTemplates() >= 1


@pytest.mark.unit
def test_cuaac_and_ruaac_smarts_are_distinct():
    """CuAAC and RuAAC templates must differ (regiochemistry correction)."""
    cuaac = ReactionSMARTS.TRIAZOLE_1_4_CuAAC["smarts"]
    ruaac = ReactionSMARTS.TRIAZOLE_1_5_RuAAC["smarts"]
    assert cuaac != ruaac


@pytest.mark.unit
def test_deprecated_triazole_formation_differs_from_canonical():
    """Legacy SMARTS must not silently alias the corrected templates."""
    deprecated = ReactionSMARTS.TRIAZOLE_FORMATION["smarts"]
    assert deprecated != ReactionSMARTS.TRIAZOLE_1_4_CuAAC["smarts"]
    assert deprecated != ReactionSMARTS.TRIAZOLE_1_5_RuAAC["smarts"]
    assert ReactionSMARTS.TRIAZOLE_FORMATION["regioisomer"] == "UNLABELLED"


@pytest.mark.unit
@pytest.mark.parametrize(
    "rxn_def,expected_regio",
    [
        (ReactionSMARTS.TRIAZOLE_1_4_CuAAC, "1,4-CuAAC"),
        (ReactionSMARTS.TRIAZOLE_1_5_RuAAC, "1,5-RuAAC"),
    ],
)
def test_triazole_templates_have_regioisomer_labels(rxn_def, expected_regio):
    assert rxn_def["regioisomer"] == expected_regio


@pytest.mark.unit
def test_cuaac_reaction_runs_on_scaffold_and_alkyne():
    """CuAAC SMARTS must produce at least one product from the taloside scaffold."""
    scaffold_mol = Chem.MolFromSmiles(SCAFFOLD)
    alkyne_mol = Chem.MolFromSmiles(BUILDING_BLOCKS_MINI[0]["smiles"])
    assert scaffold_mol is not None
    assert alkyne_mol is not None

    rxn = AllChem.ReactionFromSmarts(ReactionSMARTS.TRIAZOLE_1_4_CuAAC["smarts"])
    products = rxn.RunReactants((scaffold_mol, alkyne_mol))
    assert len(products) >= 1


@pytest.mark.unit
def test_infer_triazole_regioisomer_patterns():
    """Regioisomer inference heuristics must label known substructures."""
    assert infer_triazole_regioisomer("N2[N+]=[N-]C=C2", fallback="") == "1,4-CuAAC"
    assert infer_triazole_regioisomer("N2C=C([N-]=[N+]2", fallback="") == "1,5-RuAAC"
    assert infer_triazole_regioisomer("unknown_smiles", fallback="fallback") == "fallback"


@pytest.mark.unit
def test_glyco_library_generator_rejects_invalid_smarts():
    with pytest.raises(ValueError, match="Invalid reaction SMARTS"):
        GlycoLibraryGenerator(
            scaffold_smiles=SCAFFOLD,
            building_blocks=BUILDING_BLOCKS_MINI,
            reaction_smarts="not_a_valid_reaction>>",
            regioisomer_label="1,4-CuAAC",
            config=LibraryConfig(max_products=5),
        )
