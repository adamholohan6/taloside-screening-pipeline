#!/usr/bin/env python
"""Setup configuration for taloside-screening-pipeline package."""

from setuptools import setup, find_packages
from pathlib import Path

# Read long description from README
long_description = Path("README.md").read_text(encoding="utf-8")

setup(
    name="taloside-screening-pipeline",
    version="1.0.0",
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
        ],
    },
    python_requires=">=3.8",
    install_requires=[
        "rdkit>=2022.09.1",
        "pandas>=1.3.0",
        "numpy>=1.19.0",
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
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="chemistry rdkit cheminformatics drug-discovery admet descriptors",
)
