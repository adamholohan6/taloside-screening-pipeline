"""Compatibility wrapper for legacy imports.

The corrected implementation lives in :mod:`taloside_pipeline.glycolibrary_generator`.
This module preserves older imports that referenced ``library_generator``.
"""

from .glycolibrary_generator import *  # noqa: F401,F403
