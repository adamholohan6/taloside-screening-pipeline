"""Unit tests for phase3_docking (no Vina executable required)."""

from pathlib import Path
from unittest import mock

import pandas as pd
import pytest
from rdkit import Chem
from taloside_pipeline.phase3_docking import (
    DockingConfig,
    assert_clean_receptor_pdbqt,
    VinaDocking,
    compute_rmsd,
    embed_ligand_3d,
    extract_ligand_coords_from_pdb,
    minmax_normalize,
    mol_to_pdbqt_string,
    parse_vina_affinity,
    validate_ligand_pdbqt,
)


@pytest.fixture
def docking_config(tmp_path):
    receptor_pdb = tmp_path / "3ZSJ.pdb"
    receptor_pdbqt = tmp_path / "3ZSJ.pdbqt"
    receptor_pdbqt.write_text("REMARK receptor stub\n", encoding="utf-8")
    return DockingConfig(
        receptor_pdb=receptor_pdb,
        receptor_pdbqt=receptor_pdbqt,
        output_dir=tmp_path / "phase3_output",
        lead_csv=tmp_path / "07_lead_scored.csv",
    )


@pytest.fixture
def clean_receptor_pdbqt():
    return Path("data/docking/3ZSJ_clean.pdbqt")


@pytest.fixture
def mini_lead_df():
    return pd.DataFrame(
        {
            "compound_id": ["TEST-001", "TEST-BAD"],
            "product_smiles": [
                "c1ccccc1",
                "not_a_smiles",
            ],
            "lead_score": [0.8, 0.2],
            "generation_status": ["success", "success"],
        }
    )


@pytest.mark.unit
def test_parse_vina_affinity():
    stdout = """
    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -8.42          0          0
       2       -7.91      1.234      2.345
    """
    assert parse_vina_affinity(stdout) == pytest.approx(-8.42)


@pytest.mark.unit
def test_minmax_normalize_invert():
    s = pd.Series([-9.0, -7.0, -8.0])
    norm = minmax_normalize(s, higher_is_better=False)
    assert norm.iloc[0] == pytest.approx(1.0)  # most negative → best
    assert norm.iloc[1] == pytest.approx(0.0)


@pytest.mark.unit
def test_embed_and_pdbqt_benzene():
    mol = Chem.MolFromSmiles("c1ccccc1")
    embedded, status = embed_ligand_3d(mol)
    assert status == "success"
    assert embedded is not None
    pdbqt = mol_to_pdbqt_string(embedded, conf_id=0, name="benzene")
    assert "ROOT" in pdbqt
    assert "TORSDOF" in pdbqt
    assert "ATOM" in pdbqt


@pytest.mark.unit
def test_clean_receptor_is_clean(clean_receptor_pdbqt):
    assert_clean_receptor_pdbqt(clean_receptor_pdbqt)


@pytest.mark.unit
def test_contaminated_receptor_is_rejected():
    receptor = Path("data/docking/3ZSJ.pdbqt")
    with pytest.raises(AssertionError, match="forbidden residues"):
        assert_clean_receptor_pdbqt(receptor)


@pytest.mark.unit
def test_ligand_pdbqt_has_charges_and_torsions(tmp_path):
    """Open Babel must produce a PDBQT with non-zero Gasteiger charges and TORSDOF."""
    # n-propylbenzene: 2 rotatable bonds, simple structure, exercises the full
    # embed → obabel → validate path without depending on pipeline output files.
    mol = Chem.MolFromSmiles("CCCc1ccccc1")
    embedded, status = embed_ligand_3d(mol)
    assert status == "success"
    pdbqt = tmp_path / "propylbenzene.pdbqt"
    mol_to_pdbqt_string(embedded, conf_id=0, name="LIG")   # smoke-test the string form
    from taloside_pipeline.phase3_docking import mol_to_pdbqt_file_openbabel
    mol_to_pdbqt_file_openbabel(embedded, pdbqt, conf_id=0, name="LIG")
    validate_ligand_pdbqt(pdbqt, expected_rotatable_bonds=2)


@pytest.mark.unit
def test_all_zero_charge_pdbqt_is_rejected(tmp_path):
    """validate_ligand_pdbqt must reject a PDBQT whose charges are all zero.

    This reproduces the pre-Open Babel manual writer output where Gasteiger
    charges were never applied.  The stub is written inline so the test is
    self-contained and does not depend on leftover files from old pipeline runs.
    """
    pdbqt = tmp_path / "zero_charges.pdbqt"
    pdbqt.write_text(
        "ROOT\n"
        "ATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00"
        "     0.000 C \n"
        "ENDROOT\n"
        "TORSDOF 0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="all-zero charges"):
        validate_ligand_pdbqt(pdbqt)


@pytest.mark.unit
def test_validation_extracts_full_bgc_gal():
    coords, mol = extract_ligand_coords_from_pdb(Path("data/docking/3ZSJ.pdb"))
    assert coords.shape[0] == 23
    assert mol.GetNumAtoms() == 23


@pytest.mark.unit
def test_prepare_ligands_marks_failures(docking_config, mini_lead_df):
    docker = VinaDocking(docking_config)
    prepared = docker.prepare_ligands(mini_lead_df)
    assert prepared.loc[0, "generation_status"] == "success"
    assert prepared.loc[1, "generation_status"] == "prep_failed"
    assert Path(prepared.loc[0, "ligand_pdbqt"]).exists()


@pytest.mark.unit
def test_merge_combined_score_formula(docking_config):
    docker = VinaDocking(docking_config)
    docking_df = pd.DataFrame(
        {
            "compound_id": ["A", "B"],
            "vina_score": [-9.0, -7.0],
            "docking_status": ["dock_success", "dock_success"],
        }
    )
    scored_df = pd.DataFrame(
        {"compound_id": ["A", "B"], "lead_score": [0.2, 0.8]}
    )
    merged = docker.merge_with_lead_scores(docking_df, scored_df)
    assert "combined_score" in merged.columns
    assert merged.loc[merged["compound_id"] == "B", "combined_score"].iloc[0] > merged.loc[
        merged["compound_id"] == "A", "combined_score"
    ].iloc[0]
    out = docking_config.output_dir / "08_docking_results.csv"
    assert out.exists()


@pytest.mark.unit
def test_run_docking_mock_vina(docking_config, mini_lead_df):
    docker = VinaDocking(docking_config)
    prepared = docker.prepare_ligands(mini_lead_df.iloc[[0]])
    vina_stdout = "   1       -7.50          0          0\n"

    with mock.patch("subprocess.run") as mock_run:
        mock_run.return_value = mock.Mock(
            returncode=0,
            stdout=vina_stdout,
            stderr="",
        )
        (docking_config.output_dir / "poses" / "TEST-001_out.pdbqt").parent.mkdir(
            parents=True, exist_ok=True
        )
        docked = docker.run_docking(prepared)

    assert docked.iloc[0]["docking_status"] == "dock_success"
    assert docked.iloc[0]["vina_score"] == pytest.approx(-7.5)


@pytest.mark.unit
def test_compute_rmsd_identical():
    coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    import numpy as np

    arr = np.array(coords, dtype=float)
    assert compute_rmsd(arr, arr) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.unit
def test_validate_receptor_asserts_high_rmsd(docking_config):
    """Validation must raise when redock RMSD exceeds the 2.0 Å threshold."""
    docker = VinaDocking(docking_config)
    docker.config.rmsd_threshold_angstrom = 2.0
    docker.config.receptor_pdb.write_text("REMARK stub\n", encoding="utf-8")

    import numpy as np

    crystal = np.zeros((4, 3))
    docked = crystal + np.array([5.0, 0.0, 0.0])

    with mock.patch(
        "taloside_pipeline.phase3_docking.extract_ligand_coords_from_pdb",
        return_value=(crystal, Chem.MolFromSmiles("O")),
    ):
        with mock.patch(
            "taloside_pipeline.phase3_docking.embed_ligand_3d",
            return_value=(Chem.MolFromSmiles("O"), "success"),
        ):
            with mock.patch(
                "taloside_pipeline.phase3_docking.mol_to_pdbqt_file_openbabel",
                return_value=Path("phase3_output/validation/lactose_redock.pdbqt"),
            ):
                with mock.patch.object(docker, "_build_vina_command", return_value=["vina"]):
                    with mock.patch("subprocess.run") as mock_run:
                        mock_run.return_value = mock.Mock(
                            returncode=0,
                            stdout="   1       -6.0          0          0\n",
                            stderr="",
                        )
                        lig_dir = docker.config.output_dir / "validation"
                        lig_dir.mkdir(parents=True, exist_ok=True)
                        pose_out = lig_dir / "lactose_redock_out.pdbqt"
                        pose_out.write_text("ROOT\nENDROOT\nTORSDOF 0\n", encoding="utf-8")
                        (lig_dir / "lactose_redock.log").write_text("", encoding="utf-8")

                        with mock.patch(
                            "taloside_pipeline.phase3_docking.parse_pdbqt_coords",
                            return_value=docked,
                        ):
                            with pytest.raises(AssertionError, match="RMSD"):
                                docker.validate_receptor()
