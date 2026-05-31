"""Tests for descriptor_calculator module."""

import pytest
from rdkit import Chem
from taloside_pipeline.descriptor_calculator import (
    validate_smiles,
    calculate_descriptors,
    process_compounds,
    load_compounds_from_dict,
    calculate_library_descriptors,
)


def test_validate_smiles_valid():
    mol = validate_smiles("CCO")
    assert mol is not None
    assert isinstance(mol, Chem.Mol)


def test_validate_smiles_invalid():
    assert validate_smiles("INVALID_SMILES") is None


def test_calculate_descriptors():
    mol = validate_smiles("CCO")
    descriptors = calculate_descriptors(mol)
    assert descriptors["Molecular_Weight"] == pytest.approx(46.07, rel=0.1)
    assert descriptors["H_Donors"] == 1
    assert descriptors["H_Acceptors"] == 1
    assert descriptors["TPSA"] > 0


def test_process_compounds():
    compounds = {"Ethanol": "CCO", "Methanol": "CO", "Invalid": "INVALID"}
    df, failed_count = process_compounds(compounds)
    assert len(df) == 2
    assert failed_count == 1
    assert {"Ethanol", "Methanol"}.issubset(set(df["Compound"]))


def test_load_compounds_from_dict():
    compounds = load_compounds_from_dict()
    assert isinstance(compounds, dict)
    assert len(compounds) == 10
    assert "Talo-1" in compounds
    assert "Talo-10" in compounds


def test_calculate_library_descriptors_smoke(tmp_path):
    out = tmp_path / "taloside_descriptors.csv"
    df, failed = calculate_library_descriptors(out)
    assert failed == 0
    assert out.exists()
    assert len(df) == 10
