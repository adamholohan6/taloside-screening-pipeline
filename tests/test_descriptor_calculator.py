"""Tests for descriptor_calculator module."""
import pytest
from rdkit import Chem
from taloside_pipeline.descriptor_calculator import (
    validate_smiles,
    calculate_descriptors,
    process_compounds,
    load_compounds_from_dict,
)


def test_validate_smiles_valid():
    """Test validation of valid SMILES string."""
    smiles = "CCO"  # ethanol
    mol = validate_smiles(smiles)
    assert mol is not None
    assert isinstance(mol, Chem.Mol)


def test_validate_smiles_invalid():
    """Test validation of invalid SMILES string."""
    smiles = "INVALID_SMILES"
    mol = validate_smiles(smiles)
    assert mol is None


def test_calculate_descriptors():
    """Test descriptor calculation for a simple molecule."""
    smiles = "CCO"  # ethanol
    mol = validate_smiles(smiles)
    descriptors = calculate_descriptors(mol)
    
    assert "Molecular_Weight" in descriptors
    assert "LogP" in descriptors
    assert "H_Donors" in descriptors
    assert "H_Acceptors" in descriptors
    assert "TPSA" in descriptors
    assert "Rotatable_Bonds" in descriptors
    
    # Check some expected values for ethanol
    assert descriptors["Molecular_Weight"] == pytest.approx(46.07, rel=0.1)
    assert descriptors["H_Donors"] == 1
    assert descriptors["H_Acceptors"] == 1


def test_process_compounds():
    """Test processing a dictionary of compounds."""
    compounds = {
        "Ethanol": "CCO",
        "Methanol": "CO",
        "Invalid": "INVALID"
    }
    
    df, failed_count = process_compounds(compounds)
    
    assert len(df) == 2  # Only valid compounds
    assert failed_count == 1
    assert "Ethanol" in df["Compound"].values
    assert "Methanol" in df["Compound"].values


def test_load_compounds_from_dict():
    """Test loading the default compound library."""
    compounds = load_compounds_from_dict()
    
    assert isinstance(compounds, dict)
    assert len(compounds) == 10  # 10 taloside compounds
    assert "Talo-1" in compounds
    assert "Talo-10" in compounds
