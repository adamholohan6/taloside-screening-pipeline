"""Taloside Screening Pipeline.

Utilities for descriptor calculation, combinatorial taloside library generation,
and the Phase 2 filtering workflow.
"""

__version__ = "1.0.0"
__author__ = "Adam Holohan"

from .descriptor_calculator import (
    calculate_descriptors,
    calculate_library_descriptors,
    load_compounds_from_dict,
    process_compounds,
    validate_smiles,
)
from .glycolibrary_generator import (
    GlycoLibraryGenerator,
    LibraryConfig,
    ReactionSMARTS,
    configure_logging,
    generate_triazole_library,
)
from .phase2_integration import run_phase2_pipeline
from .phase3_docking import DockingConfig, VinaDocking, run_phase3_pipeline

# Backwards-compatible alias retained for older scripts.
Phase2Integration = run_phase2_pipeline

__all__ = [
    "__version__",
    "__author__",
    "calculate_descriptors",
    "calculate_library_descriptors",
    "load_compounds_from_dict",
    "process_compounds",
    "validate_smiles",
    "GlycoLibraryGenerator",
    "LibraryConfig",
    "ReactionSMARTS",
    "configure_logging",
    "generate_triazole_library",
    "run_phase2_pipeline",
    "Phase2Integration",
    "DockingConfig",
    "VinaDocking",
    "run_phase3_pipeline",
]
