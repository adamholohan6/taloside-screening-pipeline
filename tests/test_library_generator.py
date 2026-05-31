"""Tests for library_generator module."""
import pytest
from pathlib import Path
from taloside_pipeline.library_generator import LibraryConfig, GlycoLibraryGenerator


def test_library_config_defaults():
    """Test LibraryConfig default values."""
    config = LibraryConfig()
    
    assert config.max_products == 50000
    assert config.include_stereoisomers == True
    assert config.sanitize_products == True
    assert config.filter_hypervalent == True
    assert config.min_product_mw == 150.0
    assert config.max_product_mw == 1000.0


def test_library_config_custom():
    """Test LibraryConfig with custom values."""
    config = LibraryConfig(
        max_products=1000,
        include_stereoisomers=False,
        min_product_mw=200.0
    )
    
    assert config.max_products == 1000
    assert config.include_stereoisomers == False
    assert config.min_product_mw == 200.0


def test_library_generator_initialization():
    """Test GlycoLibraryGenerator initialization."""
    config = LibraryConfig()
    generator = GlycoLibraryGenerator(config)
    
    assert generator.config == config


def test_library_generator_output_dir():
    """Test that output directory can be set."""
    config = LibraryConfig(output_dir=Path("test_output"))
    assert config.output_dir == Path("test_output")
