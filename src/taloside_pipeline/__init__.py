"""
Taloside Screening Pipeline

A Python-based computational chemistry pipeline for calculating physicochemical 
descriptors of taloside compounds and generating virtual libraries.
"""

__version__ = "0.2.0"
__author__ = "Adam Holohan"

from .descriptor_calculator import calculate_descriptors, validate_smiles
from .library_generator import GlycoLibraryGenerator, LibraryConfig
from .phase2_integration import run_phase2_pipeline

__all__ = [
    "calculate_descriptors",
    "validate_smiles", 
    "GlycoLibraryGenerator",
    "LibraryConfig",
    "Phase2Integration",
]
