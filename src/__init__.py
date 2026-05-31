"""Taloside Screening Pipeline - Package initialization."""

__version__ = "0.2.0"
__author__ = "Adam Holohan"
__email__ = "adamholohan6@gmail.com"

from .generate_library import (
    validate_smiles,
    calculate_descriptors,
    process_compounds,
    TalosideLibraryGenerator,
)
from .glycolibrary_generator import (
    GlycoLibraryGenerator,
    LibraryConfig,
    ReactionSMARTS,
    configure_logging,
)
from .phase2_integration import (
    PAINSFilter,
    apply_lipinski_filter,
    apply_pains_filter,
    run_phase2_pipeline,
)

__all__ = [
    "validate_smiles",
    "calculate_descriptors",
    "process_compounds",
    "TalosideLibraryGenerator",
    "GlycoLibraryGenerator",
    "LibraryConfig",
    "ReactionSMARTS",
    "configure_logging",
    "PAINSFilter",
    "apply_lipinski_filter",
    "apply_pains_filter",
    "run_phase2_pipeline",
]
