#!/usr/bin/env python
"""Setup configuration for taloside-screening-pipeline package."""

import re

from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README
long_description = Path("README.md").read_text(encoding="utf-8")

# Single source of truth for the version: src/taloside_pipeline/__init__.py.
# Parsed rather than imported so that building does not require rdkit.
_init = Path("src/taloside_pipeline/__init__.py").read_text(encoding="utf-8")
_match = re.search(r"^__version__ = ['\"]([^'\"]+)['\"]", _init, re.M)
if _match is None:
    raise RuntimeError("Cannot find __version__ in src/taloside_pipeline/__init__.py")
version = _match.group(1)

setup(
    name="taloside-screening-pipeline",
    version=version,
    author="Adam Holohan",
    author_email="adamholohan6@gmail.com",
    description="Taloside virtual library generation, drug-likeness filtering, and PAINS screening pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adamholohan6/taloside-screening-pipeline",
    project_urls={
        "Bug Tracker": "https://github.com/adamholohan6/taloside-screening-pipeline/issues",
        "Documentation": "https://github.com/adamholohan6/taloside-screening-pipeline/blob/main/README.md",
    },
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            "taloside-descriptors=taloside_pipeline.descriptor_calculator:main",
            "taloside-phase2=taloside_pipeline.phase2_integration:run_phase2_pipeline",
            "taloside-phase3=taloside_pipeline.phase3_docking:run_phase3_pipeline",
        ],
    },
    python_requires=">=3.10",
    install_requires=[
        "rdkit==2026.03.2",
        "pandas>=1.3.0",
        "numpy>=1.19.0",
        "biopython>=1.80",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "pytest-xdist>=2.5.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.910",
            "isort>=5.10.0",
            "pyyaml>=6.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
    keywords="chemistry rdkit cheminformatics drug-discovery admet descriptors",
)
