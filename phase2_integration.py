"""
Phase 2 Integration: Virtual Library Generation → Lipinski Filtering → PAINS Detection

Workflow:
  1. Generate combinatorial library using GlycoLibraryGenerator
  2. Calculate physicochemical descriptors
  3. Apply Lipinski's Rule of 5 filters
  4. Screen for PAINS (Pan-Assay Interference Compounds) substructures
  5. Export filtered + flagged populations for Phase 3 analysis

PAINS Rationale:
  - PAINS are frequent hitters in high-throughput screening assays
  - Not inherently "bad" but likely to show assay artifacts
  - Important to flag and investigate separately for academic drug discovery

Reference:
  - Lipinski et al. (1997). "Drug-like properties and causes of poor solubility 
    and permeability." Adv. Drug Deliv. Rev.
  - Baell & Holloway (2010). "New Substructure Filters for Removal of 
    Pan Assay Interference Compounds (PAINS)." J. Med. Chem.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
import logging

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen

from glycolibrary_generator import (
    GlycoLibraryGenerator,
    LibraryConfig,
    ReactionSMARTS,
    configure_logging
)


logger = configure_logging(log_file=Path("phase2_integration.log"))


# ============================================================================
# PAINS FILTER DEFINITIONS
# ============================================================================

class PAINSFilter:
    """
    Pan-Assay Interference Compounds (PAINS) detection.
    
    Reference: Baell & Holloway (2010) - J. Med. Chem.
    
    PAINS are frequent hitters in HTS assays due to:
    - Michael acceptors (electrophilic reactivity)
    - Chelating groups (metal coordination)
    - Redox-active moieties (assay interference)
    - Covalent binders (non-specific protein modification)
    
    For carbohydrate drugs like taloside derivatives, some PAINS
    may be less problematic (e.g., sugars are polar, less likely to aggregate)
    but should still be tracked for experimental validation.
    """
    
    # SMARTS patterns for common PAINS
    SMARTS_PATTERNS = {
        "michael_acceptor_A": "[#6,#7;!$(*#*)]=!@[#7,#8]",  # Generic Michael acceptor
        "michael_acceptor_B": "[#6,#7]-[#6]=[#6][#6,#7,#8,#16]",  # Conjugated
        
        "aromatic_nitrogen_oxide": "[n+;H0]([*])[O-]",  # Aromatic N-oxide
        
        "quaternary_sulfonium": "[S+]([!$(*=*)],~[*])(~[*])[!$(*#*)]",  # Sulfonium
        
        "phenolic_ester": "[#6][#6]=[#6][#8][#6](=[#8])[#6,#7,#8]",  # Phenyl esters
        
        "thioesters": "[#6][#6]=[#6][#16][#6](=[#8])[#6,#7,#8]",  # Thio esters
        
        "biaryl_diazo": "[#6]~[#6][#7]=[#7]",  # Diazo compounds
        
        "aldehyde": "[#6;!$(C=O)][#6]=O",  # Aldehydes (excluding carbonyls)
        
        "unstable_nitro": "[#6][#7](=[O])[O]",  # Nitro groups on aliphatic C
        
        "alkyl_halides": "[#6][Cl,Br,I]",  # Alkyl halides (SN2 reactive)
        
        "hydrazine": "[#7][#7]",  # Hydrazines and azo
        
        "thiocyanate": "[#6][#16][#6]#[#7]",  # Thiocyanates
        
        "chromone": "[#6]1=[#6][#6](=[#8])[#8][#6]1",  # Chromones
        
        "anthracene": "[#6]1=[#6][#6]2=[#6][#6]=[#6][#6]2=[#6][#6]1",  # Anthracenes
    }
    
    @classmethod
    def compile_patterns(cls) -> Dict[str, Chem.Mol]:
        """
        Compile SMARTS patterns to RDKit molecule patterns.
        
        Returns:
            Dict mapping pattern name to compiled RDKit pattern
        """
        compiled = {}
        for name, smarts in cls.SMARTS_PATTERNS.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    compiled[name] = pattern
                else:
                    logger.warning(f"Failed to compile SMARTS: {name}")
            except Exception as e:
                logger.warning(f"Error compiling {name}: {e}")
        return compiled
    
    @classmethod
    def screen_molecule(cls, mol: Chem.Mol, compiled_patterns: Dict) -> Tuple[bool, List[str]]:
        """
        Screen molecule against PAINS patterns.
        
        Args:
            mol: RDKit molecule to screen
            compiled_patterns: Pre-compiled SMARTS patterns
            
        Returns:
            Tuple (is_pains, matched_patterns)
        """
        matched_patterns = []
        
        for pattern_name, pattern in compiled_patterns.items():
            if mol.HasSubstructMatch(pattern):
                matched_patterns.append(pattern_name)
        
        is_pains = len(matched_patterns) > 0
        return is_pains, matched_patterns


# ============================================================================
# LIPINSKI FILTER
# ============================================================================

def apply_lipinski_filter(
    df: pd.DataFrame,
    strict_mode: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Apply Lipinski's Rule of 5 filters to generated library.
    
    Lipinski's Rules for oral bioavailability:
      - MW ≤ 500 Da (≤600 for carbohydrate drugs)
      - LogP ≤ 5 (≤4 for carbohydrate drugs)
      - H-bond donors ≤ 5 (≤6 for carbohydrates)
      - H-bond acceptors ≤ 10 (≤12 for carbohydrates)
      - No sulfhydryl or phosphate groups (simplified here)
    
    Args:
        df: DataFrame with descriptor columns
        strict_mode: If True, apply stricter thresholds for carbohydrate drugs
        
    Returns:
        Tuple (passed_df, failed_df)
    """
    
    if strict_mode:
        # Carbohydrate-adjusted thresholds
        mw_threshold = 600
        logp_threshold = 4.0
        hbd_threshold = 6
        hba_threshold = 12
    else:
        # Classical Lipinski
        mw_threshold = 500
        logp_threshold = 5.0
        hbd_threshold = 5
        hba_threshold = 10
    
    passed = df[
        (df['molecular_weight'] <= mw_threshold) &
        (df['logp'] <= logp_threshold) &
        (df['h_donors'] <= hbd_threshold) &
        (df['h_acceptors'] <= hba_threshold)
    ].copy()
    
    failed = df[
        ~df.index.isin(passed.index)
    ].copy()
    
    logger.info(f"[OK] Lipinski filtering complete:")
    logger.info(f"  Passed: {len(passed)} compounds")
    logger.info(f"  Failed: {len(failed)} compounds")
    logger.info(f"  Pass rate: {100*len(passed)/len(df):.1f}%")
    
    return passed, failed


def apply_pains_filter(
    df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Screen library for PAINS substructures.
    
    Args:
        df: DataFrame with 'product_smiles' column
        
    Returns:
        Tuple (clean_df, pains_df, failed_df)
        - clean_df: Compounds with no PAINS flags
        - pains_df: Compounds containing PAINS substructures
        - failed_df: Compounds where screening failed
    """
    
    logger.info("Compiling PAINS patterns...")
    compiled_patterns = PAINSFilter.compile_patterns()
    logger.info(f"  Loaded {len(compiled_patterns)} PAINS patterns")
    
    clean_records = []
    pains_records = []
    failed_records = []
    
    logger.info("Screening library for PAINS...")
    
    for idx, row in df.iterrows():
        smiles = row['product_smiles']
        
        try:
            # Parse molecule without sanitization to avoid valence errors
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                failed_records.append(idx)
                continue
            
            # Compute rings and basic properties
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
            
            # Try to sanitize, but ignore valence errors for triazole rings
            try:
                Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
            except:
                # If sanitization fails, continue anyway - ring info is already set
                pass
            
            is_pains, matched_patterns = PAINSFilter.screen_molecule(mol, compiled_patterns)
            
            record = row.copy()
            record['pains_flag'] = is_pains
            record['pains_patterns'] = '; '.join(matched_patterns) if matched_patterns else 'None'
            
            if is_pains:
                pains_records.append(record)
            else:
                clean_records.append(record)
        
        except Exception as e:
            logger.warning(f"Failed to screen {idx}: {e}")
            failed_records.append(idx)
    
    clean_df = pd.DataFrame(clean_records)
    pains_df = pd.DataFrame(pains_records)
    failed_df = df.loc[failed_records] if failed_records else pd.DataFrame()
    
    logger.info(f"[OK] PAINS screening complete:")
    logger.info(f"  Clean (no PAINS): {len(clean_df)} compounds")
    logger.info(f"  Flagged (PAINS): {len(pains_df)} compounds")
    logger.info(f"  Failed screening: {len(failed_df)} compounds")
    
    return clean_df, pains_df, failed_df


# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def run_phase2_pipeline():
    """
    End-to-end Phase 2: Generation → Lipinski → PAINS Screening
    """
    
    # Load taloside scaffold from MSc work
    SCAFFOLD = "O=C(O[C@H]1[C@@H](OCN=[N+]=[N-])[C@@H](O)[C@@H](CO)O[C@H]1OC)C4=C([N+]([O-])=O)C=CC=C4"
    
    # Building blocks: Aromatic alkynes (terminal alkyne for click chemistry)
    BUILDING_BLOCKS = [
        {'id': 'BB-001-Ph', 'smiles': 'C#Cc1ccccc1'},
        {'id': 'BB-002-4OMe', 'smiles': 'C#Cc1ccc(OC)cc1'},
        {'id': 'BB-003-4Cl', 'smiles': 'C#Cc1ccc(Cl)cc1'},
        {'id': 'BB-004-4F', 'smiles': 'C#Cc1ccc(F)cc1'},
        {'id': 'BB-005-3Br', 'smiles': 'C#Cc1cc(Br)ccc1'},
        {'id': 'BB-006-2NO2', 'smiles': 'C#Cc1ccccc1[N+](=O)[O-]'},
        {'id': 'BB-007-Pyridine', 'smiles': 'C#Cc1ccccn1'},
        {'id': 'BB-008-Furan', 'smiles': 'C#Cc1ccoc1'},
    ]
    
    config = LibraryConfig(
        max_products=500,
        min_product_mw=250.0,
        max_product_mw=800.0,
        include_stereoisomers=True,
        output_dir=Path("phase2_output")
    )
    
    try:
        # =====================================================================
        # STEP 1: GENERATE LIBRARY
        # =====================================================================
        logger.info("=" * 70)
        logger.info("PHASE 2: VIRTUAL LIBRARY EXPANSION")
        logger.info("=" * 70)
        
        generator = GlycoLibraryGenerator(
            scaffold_smiles=SCAFFOLD,
            building_blocks=BUILDING_BLOCKS,
            reaction_smarts=ReactionSMARTS.TRIAZOLE_FORMATION["smarts"],
            config=config,
            logger_instance=logger
        )
        
        library_df = generator.generate_library()
        
        if library_df.empty:
            logger.error("No library generated!")
            return
        
        # =====================================================================
        # STEP 2: LIPINSKI FILTERING
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 2b: LIPINSKI FILTERING")
        logger.info("=" * 70)
        
        lipinski_passed, lipinski_failed = apply_lipinski_filter(
            library_df,
            strict_mode=True
        )
        
        # =====================================================================
        # STEP 3: PAINS SCREENING (on Lipinski-passed compounds)
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("PHASE 2c: PAINS SUBSTRUCTURE SCREENING")
        logger.info("=" * 70)
        
        clean_compounds, pains_compounds, failed_pains = apply_pains_filter(
            lipinski_passed
        )
        
        # =====================================================================
        # EXPORT RESULTS
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("EXPORTING RESULTS")
        logger.info("=" * 70)
        
        config.output_dir.mkdir(parents=True, exist_ok=True)
        
        # All generated compounds
        all_path = config.output_dir / "01_all_generated_compounds.csv"
        library_df.to_csv(all_path, index=False)
        logger.info(f"[OK] {all_path.name} ({len(library_df)} compounds)")
        
        # Lipinski passed
        lip_pass_path = config.output_dir / "02_lipinski_passed.csv"
        lipinski_passed.to_csv(lip_pass_path, index=False)
        logger.info(f"[OK] {lip_pass_path.name} ({len(lipinski_passed)} compounds)")
        
        # Lipinski failed
        lip_fail_path = config.output_dir / "03_lipinski_failed.csv"
        lipinski_failed.to_csv(lip_fail_path, index=False)
        logger.info(f"[OK] {lip_fail_path.name} ({len(lipinski_failed)} compounds)")
        
        # PAINS clean (high-confidence leads)
        clean_path = config.output_dir / "04_lipinski_clean_no_pains.csv"
        clean_compounds.to_csv(clean_path, index=False)
        logger.info(f"[OK] {clean_path.name} ({len(clean_compounds)} compounds - HIGH PRIORITY)")
        
        # PAINS flagged (experimental investigation)
        pains_path = config.output_dir / "05_lipinski_with_pains.csv"
        pains_compounds.to_csv(pains_path, index=False)
        logger.info(f"[OK] {pains_path.name} ({len(pains_compounds)} compounds - FLAG FOR ASSAY VALIDATION)")
        
        # Failed products from generator
        generator.export_failed_products()
        
        # =====================================================================
        # SUMMARY STATISTICS
        # =====================================================================
        logger.info("\n" + "=" * 70)
        logger.info("SUMMARY STATISTICS")
        logger.info("=" * 70)
        
        logger.info(f"\nGeneration Funnel:")
        logger.info(f"  Total generated:        {len(library_df):5} (100.0%)")
        logger.info(f"  -> Lipinski passed:      {len(lipinski_passed):5} ({100*len(lipinski_passed)/len(library_df):5.1f}%)")
        logger.info(f"  -> No PAINS (Clean):     {len(clean_compounds):5} ({100*len(clean_compounds)/len(library_df):5.1f}%)")
        logger.info(f"\n  Lipinski failed:        {len(lipinski_failed):5} ({100*len(lipinski_failed)/len(library_df):5.1f}%)")
        logger.info(f"  PAINS flagged:          {len(pains_compounds):5} ({100*len(pains_compounds)/len(lipinski_passed):5.1f}% of Lipinski-passed)")
        
        logger.info(f"\nTop PAINS Patterns (if any):")
        if not pains_compounds.empty and 'pains_patterns' in pains_compounds.columns:
            pains_counts = {}
            for patterns_str in pains_compounds['pains_patterns']:
                if patterns_str != 'None':
                    for pattern in patterns_str.split('; '):
                        pains_counts[pattern] = pains_counts.get(pattern, 0) + 1
            
            for pattern, count in sorted(pains_counts.items(), key=lambda x: x[1], reverse=True)[:5]:
                logger.info(f"  - {pattern}: {count} compound(s)")
        
        logger.info("\n" + "=" * 70)
        logger.info("[OK] PHASE 2 PIPELINE COMPLETE")
        logger.info("=" * 70)
        
        return {
            'all': library_df,
            'lipinski_passed': lipinski_passed,
            'lipinski_failed': lipinski_failed,
            'clean': clean_compounds,
            'pains': pains_compounds,
        }
    
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    run_phase2_pipeline()
