"""
Phase 2 Integration: Virtual Library Generation → Lipinski Filtering → PAINS Detection
========================================================================================

Workflow:
  1. Generate combinatorial library using GlycoLibraryGenerator (both regioisomers)
  2. Standardise structures and calculate physicochemical descriptors
  3. Apply Lipinski's Rule of 5 with carbohydrate-adjusted thresholds
  4. Screen for PAINS using RDKit's validated FilterCatalog
  5. Compute lead score (defined formula, not implicit)
  6. Export filtered + flagged populations for downstream analysis

PAINS Fallback Policy (reviewer correction):
  The original pipeline treated all compounds as PAINS-clean when the
  FilterCatalog failed to load ('conservative' was a misnomer - permissive
  is correct). Compounds with PAINS status that cannot be determined are now
  exported to 06_pains_undetermined.csv rather than silently passed as clean.

Regiochemistry (reviewer correction):
  Products are generated separately for CuAAC (1,4-triazole) and RuAAC
  (1,5-triazole) and tagged accordingly. The original single SMARTS produced
  unlabelled mixed regioisomers.

Lead Score (reviewer correction):
  The lead score is now explicitly defined (formula documented in Section 2.6
  of the manuscript) rather than implicitly computed in the notebook only.

References:
  - Lipinski et al. (1997). Adv. Drug Deliv. Rev.
  - Baell & Holloway (2010). J. Med. Chem.
  - RDKit FilterCatalog (validated PAINS substructures)
  - Veber et al. (2002). J. Med. Chem. (TPSA threshold)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, FilterCatalog

try:
    from .glycolibrary_generator import (
        GlycoLibraryGenerator,
        LibraryConfig,
        ReactionSMARTS,
        configure_logging,
        generate_triazole_library,
    )
except ImportError:  # pragma: no cover - allows running as a standalone script
    from glycolibrary_generator import (  # type: ignore
        GlycoLibraryGenerator,
        LibraryConfig,
        ReactionSMARTS,
        configure_logging,
        generate_triazole_library,
    )


logger = configure_logging(log_file=Path("phase2_integration.log"))


# ============================================================================
# PAINS FILTER — CORRECTED FALLBACK
# ============================================================================

class PAINSFilter:
    """
    PAINS detection using RDKit's FilterCatalog (PAINS_A/B/C substructure sets).

    Fallback policy (reviewer correction):
      If the catalog cannot be loaded, screen_molecule() returns (None, [])
      to signal UNDETERMINED status. Callers must route these to a separate
      CSV rather than treating them as clean.

    Limitations:
      PAINS_A/B/C does not cover all assay interference mechanisms.
      Colloidal aggregation, autofluorescence, redox cycling, and metal
      chelation are outside the scope of this filter and should be assessed
      separately before proceeding to fluorescence polarisation or ITC assays.
    """

    _catalog: Optional[FilterCatalog.FilterCatalog] = None
    _catalog_loaded: bool = False

    @classmethod
    def reset(cls) -> None:
        """
        CRITICAL FIX: Reset module-level state for test isolation.
        
        Call this in test setUp() and tearDown() to avoid test pollution
        from stateful _catalog and _catalog_loaded variables.
        
        Example:
            def setUp(self):
                PAINSFilter.reset()
            
            def tearDown(self):
                PAINSFilter.reset()
        """
        cls._catalog = None
        cls._catalog_loaded = False
        logger.debug("PAINSFilter state reset")

    @classmethod
    def _get_catalog(cls) -> Optional[FilterCatalog.FilterCatalog]:
        if cls._catalog_loaded:
            return cls._catalog

        logger.info("Initialising RDKit PAINS catalog...")
        try:
            params = FilterCatalog.FilterCatalogParams()
            # Handle different RDKit versions
            pains_attr = None
            catalog_enum = getattr(FilterCatalog.FilterCatalogParams, "FilterCatalogs", None)
            if catalog_enum is not None:
                for attr_name in ("PAINS", "PAINS_A", "FILTER_PAINS"):
                    if hasattr(catalog_enum, attr_name):
                        pains_attr = getattr(catalog_enum, attr_name)
                        break
            if pains_attr is None:
                raise AttributeError("No PAINS constant found in FilterCatalogParams.FilterCatalogs")
            params.AddCatalog(pains_attr)
            cls._catalog = FilterCatalog.FilterCatalog(params)
            logger.info(f"  Loaded PAINS catalog ({cls._catalog.GetNumEntries()} entries)")
        except Exception as e:
            logger.warning(
                f"Could not load PAINS catalog: {e}. "
                "Compounds will be flagged as PAINS status UNDETERMINED."
            )
            cls._catalog = None

        cls._catalog_loaded = True
        return cls._catalog

    @classmethod
    def screen_molecule(
        cls, mol: Chem.Mol
    ) -> Tuple[Optional[bool], List[str]]:
        """
        Screen molecule for PAINS substructures.

        Returns:
            (True, [matched patterns])  — PAINS detected
            (False, [])                 — clean
            (None, [])                  — catalog unavailable; status UNDETERMINED
        """
        catalog = cls._get_catalog()
        if catalog is None:
            return None, []   # UNDETERMINED — do NOT treat as clean

        try:
            entries = catalog.GetMatches(mol)
            matches = [e.GetDescription() for e in entries]
            return len(matches) > 0, matches
        except Exception as e:
            logger.debug(f"PAINS screening error: {e}")
            return None, []  # UNDETERMINED on error


# ============================================================================
# LIPINSKI FILTER
# ============================================================================

def apply_lipinski_filter(
    df: pd.DataFrame,
    strict_mode: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply Lipinski's Rule of 5 with carbohydrate-adjusted thresholds.

    Threshold rationale:
      Glycomimetics frequently exceed the standard Ro5 MW (500 Da) and HBA (10)
      limits due to their polyhydroxylated cores. Adjusted thresholds follow
      precedent for carbohydrate-based drugs (MW ≤600, HBD ≤6, HBA ≤12).
      Note: LogP ≤4 is more conservative than the standard ≤5.

    TPSA note:
      TPSA is not part of Lipinski's Ro5; the Veber criterion (≤140 Å²) is a
      separate heuristic for passive membrane permeation. All taloside compounds
      have TPSA > 140 Å². This is acceptable for an extracellular target
      (Gal-3 is secreted) but precludes claims of oral bioavailability.
    """
    if strict_mode:
        mw_threshold  = 600
        logp_threshold = 4.0
        hbd_threshold  = 6
        hba_threshold  = 12
    else:
        mw_threshold  = 500
        logp_threshold = 5.0
        hbd_threshold  = 5
        hba_threshold  = 10

    mask = (
        (df['molecular_weight'] <= mw_threshold) &
        (df['logp']             <= logp_threshold) &
        (df['h_donors']         <= hbd_threshold) &
        (df['h_acceptors']      <= hba_threshold)
    )
    passed = df[mask].copy()
    failed = df[~mask].copy()

    logger.info(f"[OK] Lipinski filtering (carbohydrate-adjusted thresholds):")
    logger.info(f"  MW ≤{mw_threshold} | LogP ≤{logp_threshold} | HBD ≤{hbd_threshold} | HBA ≤{hba_threshold}")
    logger.info(f"  Passed: {len(passed)}  Failed: {len(failed)}  "
                f"Pass rate: {100*len(passed)/max(len(df),1):.1f}%")

    return passed, failed


# ============================================================================
# PAINS SCREENING
# ============================================================================

def apply_pains_filter(
    df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Screen compounds for PAINS substructures with ROBUST empty DataFrame handling.

    CRITICAL FIX: Check for empty input DataFrame and preserve schema in output.

    Returns three DataFrames:
      clean_df       — no PAINS alerts
      pains_df       — PAINS alerts detected
      undetermined_df — catalog unavailable or screening error (NOT treated as clean)

    The 'undetermined' category corrects the original fallback which silently
    passed all compounds when the catalog failed to load.
    """
    if df.empty:
        logger.warning("apply_pains_filter() received empty DataFrame")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    clean_records        = []
    pains_records        = []
    undetermined_records = []

    logger.info("Screening library for PAINS (RDKit FilterCatalog PAINS_A/B/C)...")

    for idx, row in df.iterrows():
        smiles = row.get('product_smiles', '')
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                row_copy = row.copy()
                row_copy['pains_flag']     = None
                row_copy['pains_status']   = 'UNDETERMINED - invalid SMILES'
                row_copy['pains_patterns'] = ''
                undetermined_records.append(row_copy)
                continue

            try:
                Chem.SanitizeMol(
                    mol,
                    Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES ^ Chem.SANITIZE_VALENCE
                )
            except Exception:
                Chem.FastFindRings(mol)

            try:
                mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass

            is_pains, matched = PAINSFilter.screen_molecule(mol)

            record = row.copy()
            record['pains_patterns'] = '; '.join(matched) if matched else ''

            if is_pains is None:
                # Catalog unavailable — UNDETERMINED, not clean
                record['pains_flag']   = None
                record['pains_status'] = 'UNDETERMINED - catalog unavailable'
                undetermined_records.append(record)
            elif is_pains:
                record['pains_flag']   = True
                record['pains_status'] = 'PAINS'
                pains_records.append(record)
            else:
                record['pains_flag']   = False
                record['pains_status'] = 'CLEAN'
                clean_records.append(record)

        except Exception as e:
            logger.warning(f"Failed to screen compound at index {idx}: {e}")
            row_copy = row.copy()
            row_copy['pains_flag']     = None
            row_copy['pains_status']   = f'UNDETERMINED - error: {e}'
            row_copy['pains_patterns'] = ''
            undetermined_records.append(row_copy)

    # CRITICAL FIX: Check if lists are non-empty before creating DataFrames
    # This prevents crashes from empty pd.concat operations
    if clean_records:
        clean_df = pd.DataFrame(clean_records)
    else:
        # Preserve schema even when empty
        clean_df = pd.DataFrame(columns=list(df.columns) + ['pains_flag', 'pains_status', 'pains_patterns'])

    if pains_records:
        pains_df = pd.DataFrame(pains_records)
    else:
        pains_df = pd.DataFrame(columns=list(df.columns) + ['pains_flag', 'pains_status', 'pains_patterns'])

    if undetermined_records:
        undetermined_df = pd.DataFrame(undetermined_records)
    else:
        undetermined_df = pd.DataFrame(columns=list(df.columns) + ['pains_flag', 'pains_status', 'pains_patterns'])

    logger.info(f"[OK] PAINS screening complete:")
    logger.info(f"  Clean (no PAINS):  {len(clean_df)}")
    logger.info(f"  Flagged (PAINS):   {len(pains_df)}")
    logger.info(f"  UNDETERMINED:      {len(undetermined_df)}  ← see 06_pains_undetermined.csv")

    return clean_df, pains_df, undetermined_df


# ============================================================================
# LEAD SCORING
# ============================================================================

def compute_lead_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute a composite lead score for compound ranking.

    Formula (defined here explicitly; previously only implicit in notebook):
      lead_score = 0.4 × (1 − norm_TPSA)
                 + 0.3 × (1 − norm_MW)
                 + 0.2 × (1 − norm_LogP)
                 + 0.1 × (1 − norm_RotBonds)

    where norm(x) = (x − min) / (max − min) over the library.

    Rationale: higher weight on TPSA and MW because these two descriptors
    most strongly correlate with permeability and synthetic tractability for
    extracellular glycomimetic targets. Lower LogP is preferred for solubility.
    RotBonds contributes least because conformational flexibility is less
    critical for a rigid sugar scaffold.

    Note: this score is a heuristic for within-library ranking only. It does
    not incorporate target affinity. Docking scores should supersede it.
    """
    df = df.copy()
    for col, norm_col in [
        ('tpsa',            'norm_tpsa'),
        ('molecular_weight','norm_mw'),
        ('logp',            'norm_logp'),
        ('rotatable_bonds', 'norm_rot'),
    ]:
        mn, mx = df[col].min(), df[col].max()
        df[norm_col] = (df[col] - mn) / (mx - mn) if mx > mn else 0.0

    df['lead_score'] = (
        0.4 * (1 - df['norm_tpsa']) +
        0.3 * (1 - df['norm_mw'])   +
        0.2 * (1 - df['norm_logp']) +
        0.1 * (1 - df['norm_rot'])
    )
    df = df.drop(columns=['norm_tpsa','norm_mw','norm_logp','norm_rot'])
    return df.sort_values('lead_score', ascending=False)


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def run_phase2_pipeline():
    """
    End-to-end Phase 2 pipeline:
      Generation (both regioisomers, tagged) → Lipinski → PAINS → Lead scoring → Export
    """

    SCAFFOLD = (
        "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)"
        "C4=C([N+]([O-])=O)C=CC=C4"
    )

    BUILDING_BLOCKS = [
        {'id': 'BB-001-Ph',      'smiles': 'C#Cc1ccccc1'},
        {'id': 'BB-002-4OMe',    'smiles': 'C#Cc1ccc(OC)cc1'},
        {'id': 'BB-003-4Cl',     'smiles': 'C#Cc1ccc(Cl)cc1'},
        {'id': 'BB-004-4F',      'smiles': 'C#Cc1ccc(F)cc1'},
        {'id': 'BB-005-3Br',     'smiles': 'C#Cc1cc(Br)ccc1'},
        {'id': 'BB-006-2NO2',    'smiles': 'C#Cc1ccccc1[N+](=O)[O-]'},
        {'id': 'BB-007-Pyridine','smiles': 'C#Cc1ccccn1'},
        {'id': 'BB-008-Furan',   'smiles': 'C#Cc1ccoc1'},
    ]

    config = LibraryConfig(
        max_products=500,
        min_product_mw=250.0,
        max_product_mw=800.0,
        include_stereoisomers=True,
        filter_hypervalent=True,
        output_dir=Path("phase2_output")
    )

    logger.info("=" * 70)
    logger.info("PHASE 2: VIRTUAL LIBRARY EXPANSION (corrected pipeline)")
    logger.info("=" * 70)

    # ── Step 1: Generate both regioisomers, tagged ────────────────────────
    library_df = generate_triazole_library(
        scaffold_smiles=SCAFFOLD,
        building_blocks=BUILDING_BLOCKS,
        config=config,
        logger_instance=logger
    )

    if library_df.empty:
        logger.error("No library generated — check scaffold SMILES and building blocks.")
        return

    logger.info(f"\nTotal enumerated products: {len(library_df)}")
    regio_counts = library_df['regioisomer'].value_counts()
    for label, n in regio_counts.items():
        logger.info(f"  {label}: {n} compounds")

    # ── Step 2: Lipinski filtering ────────────────────────────────────────
    logger.info("\n" + "=" * 70)
    logger.info("PHASE 2b: LIPINSKI FILTERING (carbohydrate-adjusted)")
    logger.info("=" * 70)
    lipinski_passed, lipinski_failed = apply_lipinski_filter(library_df, strict_mode=True)

    # ── Step 3: PAINS screening ───────────────────────────────────────────
    logger.info("\n" + "=" * 70)
    logger.info("PHASE 2c: PAINS SCREENING (PAINS_A/B/C)")
    logger.info("=" * 70)
    clean_df, pains_df, undetermined_df = apply_pains_filter(lipinski_passed)

    # ── Step 4: Lead scoring (on clean or all Lipinski-passed) ─────────────
    logger.info("\n" + "=" * 70)
    logger.info("PHASE 2d: LEAD SCORING")
    logger.info("=" * 70)

    # Score over all Lipinski-passed so ranking spans the full set
    scored_df = compute_lead_scores(lipinski_passed)
    # Re-attach PAINS status
    if 'pains_status' in clean_df.columns:
        scored_df = scored_df.merge(
            lipinski_passed[['compound_id']].assign(
                pains_status=lambda x: x['compound_id'].map(
                    {r['compound_id']: r.get('pains_status', '') 
                     for _, r in pd.concat([clean_df, pains_df, undetermined_df]).iterrows()}
                ) if not all(df.empty for df in [clean_df, pains_df, undetermined_df]) else ''
            ),
            on='compound_id', how='left'
        )

    logger.info(f"Lead score range: {scored_df['lead_score'].min():.3f} – "
                f"{scored_df['lead_score'].max():.3f}")
    logger.info(f"Top 5 compounds:")
    for _, row in scored_df.head(5).iterrows():
        logger.info(f"  {row['compound_id']} | {row.get('regioisomer','')} | "
                    f"score={row['lead_score']:.3f} | MW={row['molecular_weight']:.1f} | "
                    f"LogP={row['logp']:.2f} | TPSA={row['tpsa']:.1f}")

    # ── Step 5: Export ────────────────────────────────────────────────────
    logger.info("\n" + "=" * 70)
    logger.info("EXPORTING RESULTS")
    logger.info("=" * 70)

    config.output_dir.mkdir(parents=True, exist_ok=True)

    exports = [
        (library_df,      "01_all_generated_compounds.csv"),
        (lipinski_passed, "02_lipinski_passed.csv"),
        (lipinski_failed, "03_lipinski_failed.csv"),
        (clean_df,        "04_lipinski_clean_no_pains.csv"),
        (pains_df,        "05_lipinski_with_pains.csv"),
        (undetermined_df, "06_pains_undetermined.csv"),   # NEW — replaces silent pass
        (scored_df,       "07_lead_scored.csv"),
    ]

    for df_export, fname in exports:
        path = config.output_dir / fname
        df_export.to_csv(path, index=False)
        logger.info(f"[OK] {fname} ({len(df_export)} compounds)")

    # ── Summary ──────────────────────────────────────────────────────────
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY STATISTICS")
    logger.info("=" * 70)
    total = len(library_df)
    logger.info(f"\nGeneration Funnel:")
    logger.info(f"  Total enumerated:         {total:5} (100.0%)")
    logger.info(f"  → 1,4-CuAAC products:     "
                f"{len(library_df[library_df['regioisomer']=='1,4-CuAAC']):5}  (synthesisable by standard click)")
    logger.info(f"  → 1,5-RuAAC products:     "
                f"{len(library_df[library_df['regioisomer']=='1,5-RuAAC']):5}  (requires Ru catalysis)")
    logger.info(f"  → Lipinski passed:        {len(lipinski_passed):5} ({100*len(lipinski_passed)/max(total,1):.1f}%)")
    logger.info(f"     → PAINS clean:         {len(clean_df):5}")
    logger.info(f"     → PAINS flagged:       {len(pains_df):5}")
    logger.info(f"     → PAINS undetermined:  {len(undetermined_df):5}  ← inspect manually")
    logger.info(f"  → Lipinski failed:        {len(lipinski_failed):5}")

    logger.info("\n" + "=" * 70)
    logger.info("[OK] PHASE 2 PIPELINE COMPLETE")
    logger.info("=" * 70)

    return {
        'all':           library_df,
        'lipinski_passed': lipinski_passed,
        'lipinski_failed': lipinski_failed,
        'clean':         clean_df,
        'pains':         pains_df,
        'undetermined':  undetermined_df,
        'scored':        scored_df,
    }


if __name__ == "__main__":
    run_phase2_pipeline()
