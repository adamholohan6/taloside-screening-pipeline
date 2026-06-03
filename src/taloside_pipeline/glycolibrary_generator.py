"""
Taloside Screening Pipeline - Phase 2: Virtual Library Expansion Engine
========================================================================

Generates combinatorial chemical libraries from core scaffolds and building block fragments
using RDKit's Chemical Reactions framework (SMARTS-based).

Features:
  - Combinatorial product generation from scaffold + building blocks
  - Stereochemistry preservation and validation
  - Hypervalency and valence error detection
  - Duplicate and chimeric compound filtering via InChIKey (version-robust)
  - Regioisomer tagging (1,4-CuAAC vs 1,5-RuAAC) for triazole products
  - CSV export with traceable lineage (scaffold ID, building block ID)
  - Comprehensive error logging to dead-letter queue

CHANGELOG (reviewer-driven corrections):
  - BUG FIX: Replaced single ambiguous triazole SMARTS with two separate templates
    (CuAAC → 1,4-product; RuAAC → 1,5-product). Original single SMARTS produced
    both regioisomers without labelling them, misrepresenting synthetic accessibility.
  - BUG FIX: Deduplication now uses stdInChIKey (via Chem.MolToInchiKey) instead of
    canonical SMILES. Canonical SMILES with kekulisation skipped can differ for
    identical molecules across RDKit versions.
  - BUG FIX: Products are tagged with 'regioisomer' column ('1,4-CuAAC' or '1,5-RuAAC').
  - BUG FIX: MolStandardize cleanup applied before descriptor calculation.

Author: Adam Holohan
License: MIT
Requires: RDKit >= 2022.09.1, Pandas >= 1.3.0
"""

import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, inchi
from rdkit.Chem.MolStandardize import rdMolStandardize


# ============================================================================
# CONFIGURATION & LOGGING
# ============================================================================

@dataclass
class LibraryConfig:
    """Configuration container for library generation parameters."""
    max_products: int = 50000
    include_stereoisomers: bool = True
    sanitize_products: bool = True
    filter_hypervalent: bool = True
    min_product_mw: float = 150.0
    max_product_mw: float = 1000.0
    output_dir: Path = Path("library_output")


def configure_logging(log_file: Optional[Path] = None) -> logging.Logger:
    logger = logging.getLogger("GlycoLibraryGenerator")
    logger.setLevel(logging.DEBUG)
    if logger.handlers:
        return logger  # avoid duplicate handlers on reimport

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(
        logging.Formatter('[%(levelname)-8s] %(name)s: %(message)s')
    )
    logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(
            logging.Formatter(
                '[%(asctime)s] [%(levelname)-8s] %(name)s: %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
        )
        logger.addHandler(file_handler)

    return logger


logger = configure_logging()


# ============================================================================
# CORE VALIDATION UTILITIES
# ============================================================================

def validate_smiles(smiles: str) -> Optional[Chem.Mol]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        return mol
    except Exception as e:
        logger.debug(f"SMILES validation error for '{smiles}': {e}")
        return None


def standardize_mol(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """
    Apply MolStandardize cleanup (Cleanup → Normalize → Reionize → Uncharger)
    and canonicalize tautomers before descriptor calculation.

    This ensures consistent cLogP and TPSA values regardless of the input
    representation, correcting a weakness in the original pipeline where
    descriptors were calculated on unsanitised/non-standardised structures.
    """
    try:
        clean = rdMolStandardize.Cleanup(mol)
        clean = rdMolStandardize.Normalize(clean)
        clean = rdMolStandardize.Reionize(clean)
        uncharger = rdMolStandardize.Uncharger()
        clean = uncharger.uncharge(clean)
        te = rdMolStandardize.TautomerEnumerator()
        clean = te.Canonicalize(clean)
        return clean
    except Exception as e:
        logger.debug(f"Standardization failed: {e}")
        return mol  # fall back to unstandardised rather than dropping compound


def mol_to_inchikey(mol: Chem.Mol) -> Optional[str]:
    """
    Compute InChIKey for deduplication. More robust than canonical SMILES
    across RDKit versions and kekulisation settings.
    """
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None



# ============================================================================
# REACTION SMARTS DEFINITIONS
# ============================================================================

class ReactionSMARTS:
    """
    Container for standard reaction SMARTS patterns for carbohydrate chemistry.

    TRIAZOLE REGIOCHEMISTRY NOTE (reviewer correction):
      Copper-catalysed azide-alkyne cycloaddition (CuAAC) exclusively produces
      1,4-disubstituted 1,2,3-triazoles. Ruthenium-catalysed AAC (RuAAC) gives
      1,5-disubstituted products. The original pipeline used a single SMARTS that
      generated both regioisomers without labelling them, overstating synthetic
      accessibility under standard click conditions. Two separate templates are
      now defined and used; products are tagged accordingly.
    """

    # CuAAC: exclusively 1,4-disubstituted triazole
    # Product template pins ring N2 and N3 to NEUTRAL via explicit [N+0] markers.
    # Rationale: formal charges from the azide reactant ([N+], [N-]) propagate
    # into product atoms whose charge is left unspecified in the SMARTS (observed
    # under RDKit 2026.03.x). The inherited [N-] then has bond-order sum = 3,
    # which exceeds RDKit's permitted valence of 2 for N-, so SanitizeMol raises
    # "Explicit valence for atom N, 3, is greater than permitted" and every
    # product is dropped. The [N+0] markers are a required defense: they force
    # neutrality at the template level rather than relying on it (the previous
    # "all-Kekule with neutral atoms" comment was aspirational, not enforced).
    # With neutral ring nitrogens, SanitizeMol perceives aromaticity from the
    # electron count and produces clean, round-trippable canonical SMILES.
    TRIAZOLE_1_4_CuAAC = {
        "name": "triazole_1_4_CuAAC",
        "regioisomer": "1,4-CuAAC",
        "smarts": "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5][N:1][N+0:2]=[N+0:3]1",
        "description": (
            "CuAAC: 1,3-dipolar azide-alkyne cycloaddition giving the "
            "1,4-disubstituted 1,2,3-triazole exclusively. "
            "C4 of alkyne bonds to C-4 of triazole (adjacent to N-1). "
            "Ring N2/N3 are pinned neutral with explicit [N+0] markers to stop "
            "reactant azide charges propagating into the product; aromaticity is "
            "then perceived by SanitizeMol from electron count."
        )
    }

    # RuAAC: exclusively 1,5-disubstituted triazole
    # Same design as CuAAC: ring N2/N3 pinned neutral via explicit [N+0] markers
    # to prevent reactant azide charge propagation into the product template.
    TRIAZOLE_1_5_RuAAC = {
        "name": "triazole_1_5_RuAAC",
        "regioisomer": "1,5-RuAAC",
        "smarts": "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:5]1=[C:4][N:1][N+0:2]=[N+0:3]1",
        "description": (
            "RuAAC: 1,3-dipolar azide-alkyne cycloaddition giving the "
            "1,5-disubstituted 1,2,3-triazole exclusively under Ru catalysis. "
            "C5 of alkyne bonds to C-5 of triazole (adjacent to N-1). "
            "Ring N2/N3 are pinned neutral with explicit [N+0] markers, same "
            "required defense against reactant charge propagation as CuAAC."
        )
    }

    # Convenience tuple: both regioisomers
    TRIAZOLE_BOTH = (TRIAZOLE_1_4_CuAAC, TRIAZOLE_1_5_RuAAC)

    # Legacy alias retained for backwards compatibility - now deprecated.
    # DO NOT USE for new work; use TRIAZOLE_1_4_CuAAC and TRIAZOLE_1_5_RuAAC.
    # Will be removed in v3.0.
    TRIAZOLE_FORMATION = {
        "name": "triazole_formation_DEPRECATED",
        "regioisomer": "UNLABELLED",
        "smarts": "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5][N:1][N:2]=[N:3]1",
        "description": (
            "DEPRECATED - produces unlabelled regioisomers. "
            "Use TRIAZOLE_1_4_CuAAC and TRIAZOLE_1_5_RuAAC instead."
        )
    }

    # Other chemistry (unchanged)
    AMIDE_COUPLING = {
        "name": "amide_coupling",
        "regioisomer": None,
        "smarts": "[N:1].[C:2](=O)[O:3]>>[N:1]C(=O)[C:2]",
        "description": "Amide bond formation between primary amine and carboxylic acid."
    }

    SUZUKI_COUPLING = {
        "name": "suzuki_coupling",
        "regioisomer": None,
        "smarts": "[c:1][Br,I].[c:2]B(O)O>>[c:1][c:2]",
        "description": "Suzuki-Miyaura cross-coupling: aryl halide + aryl boronic acid."
    }


# ============================================================================
# MAIN GLYCOLIBRARY GENERATOR CLASS
# ============================================================================

class GlycoLibraryGenerator:
    """
    High-throughput combinatorial library generation engine for glycomimetic scaffolds.

    Generates products from a core scaffold + building blocks using SMARTS-encoded
    reaction templates. For triazole libraries, pass both TRIAZOLE_1_4_CuAAC and
    TRIAZOLE_1_5_RuAAC to enumerate both regioisomers with synthetic-route tags.
    """

    def __init__(
        self,
        scaffold_smiles: str,
        building_blocks: List[Dict[str, str]],
        reaction_smarts: str,
        regioisomer_label: str = "",
        config: LibraryConfig = None,
        logger_instance: logging.Logger = None
    ):
        """
        Args:
            scaffold_smiles:    SMILES for core scaffold
            building_blocks:    List of {'id': str, 'smiles': str}
            reaction_smarts:    SMARTS for the chemical reaction
            regioisomer_label:  Tag for products (e.g. '1,4-CuAAC', '1,5-RuAAC', '')
            config:             LibraryConfig (uses defaults if None)
            logger_instance:    Logger (uses module logger if None)
        """
        self.logger = logger_instance or logger
        self.config = config or LibraryConfig()
        self.regioisomer_label = regioisomer_label

        # Validate scaffold
        self.scaffold_smiles = scaffold_smiles
        self.scaffold_mol = validate_smiles(scaffold_smiles)
        if self.scaffold_mol is None:
            raise ValueError(f"Invalid scaffold SMILES: {scaffold_smiles}")
        self.logger.info(
            f"[OK] Scaffold loaded ({self.scaffold_mol.GetNumAtoms()} atoms)"
        )

        # Validate building blocks
        self.building_blocks = []
        invalid_count = 0
        for bb in building_blocks:
            bb_mol = validate_smiles(bb.get('smiles', ''))
            if bb_mol is None:
                self.logger.warning(f"[X] Invalid BB SMILES (ID: {bb.get('id','?')})")
                invalid_count += 1
                continue
            self.building_blocks.append({**bb, 'mol': bb_mol})

        if invalid_count:
            self.logger.warning(
                f"Skipped {invalid_count} invalid building block(s). "
                f"Proceeding with {len(self.building_blocks)} valid blocks."
            )
        if not self.building_blocks:
            raise ValueError("No valid building blocks provided.")
        self.logger.info(f"[OK] Loaded {len(self.building_blocks)} building blocks")

        # Parse reaction SMARTS
        self.reaction_smarts = reaction_smarts
        self.rxn = self._parse_reaction_smarts(reaction_smarts)
        if self.rxn is None:
            raise ValueError(f"Invalid reaction SMARTS: {reaction_smarts}")

        self.products: List[Dict] = []
        self.failed_products: List[Dict] = []

    def _parse_reaction_smarts(self, smarts: str) -> Optional[object]:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None or rxn.GetNumReactantTemplates() == 0:
                return None
            return rxn
        except Exception as e:
            self.logger.error(f"Failed to parse reaction SMARTS: {e}")
            return None

    def _apply_reaction(
        self, scaffold_mol: Chem.Mol, bb_mol: Chem.Mol
    ) -> List[Tuple[Chem.Mol, ...]]:
        try:
            return self.rxn.RunReactants((scaffold_mol, bb_mol))
        except Exception as e:
            self.logger.debug(f"Reaction failed: {e}")
            return []

    def _sanitize_and_validate_product(
        self, product_mol: Chem.Mol, scaffold_id: str, bb_id: str
    ) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        Sanitize and validate a product molecule from RunReactants().

        With the corrected all-Kekule product SMARTS (uppercase atoms, no formal
        charges on ring nitrogens), standard SanitizeMol correctly perceives
        triazole aromaticity from the electron count and produces clean,
        round-trippable SMILES.

        Note: the aggressive ring-nitrogen repair loop that previously preceded
        SanitizeMol has been removed. It was only needed to mask broken product
        SMARTS that included aromatic atoms or formal charges in the product
        template. With clean SMARTS, no pre-sanitization repair is needed.
        """
        try:
            Chem.SanitizeMol(product_mol)
        except Exception as e:
            return False, None, f"Sanitization failed: {e}"

        mol_wt = Descriptors.MolWt(product_mol)
        if not (self.config.min_product_mw <= mol_wt <= self.config.max_product_mw):
            return False, None, (
                f"MW out of range ({mol_wt:.2f}; "
                f"expected {self.config.min_product_mw}–{self.config.max_product_mw})"
            )

        try:
            product_smiles = Chem.MolToSmiles(
                product_mol, isomericSmiles=self.config.include_stereoisomers
            )
        except Exception as e:
            return False, None, f"Failed to generate SMILES: {e}"

        return True, product_smiles, None

    def generate_library(self) -> pd.DataFrame:
        """
        Execute combinatorial library generation.

        Returns a DataFrame with columns including 'regioisomer' tagging each
        product as '1,4-CuAAC', '1,5-RuAAC', or '' (for non-triazole reactions).
        Deduplication uses InChIKey for version-robust uniqueness.
        """
        self.logger.info("=" * 70)
        self.logger.info("INITIATING COMBINATORIAL LIBRARY GENERATION")
        self.logger.info(f"Regioisomer label: '{self.regioisomer_label}'")
        self.logger.info("=" * 70)

        seen_inchikeys: Set[str] = set()
        product_count = 0

        # Pre-compile triazole geometry-filter patterns.
        #
        # [C:4]#[C:5] in the reaction SMARTS matches a terminal alkyne in BOTH
        # orientations (terminal-C as :4 or as :5), so RunReactants produces two
        # products per building block even though each catalyst gives only one
        # regioisomer in practice:
        #
        #   CuAAC → 1,4-disubstituted: aryl lands on C4 (two atoms from N1)
        #   RuAAC → 1,5-disubstituted: aryl lands on C5 (adjacent to N1)
        #
        # The patterns below match the N1-bearing-scaffold path around the ring
        # to the aryl-bearing carbon:
        #
        #   _pat_14: [#6]-N1-N2-N3-C4(aryl)-C5 (aryl on C4 → 1,4)
        #   _pat_15: [#6]-N1-N2-N3-C4-C5(aryl) (aryl on C5 → 1,5)
        #
        # Products whose geometry doesn't match the template label are routed to
        # failed_products with reason 'geometry_mismatch'. This makes the library
        # contents honest: every compound labelled 1,4-CuAAC actually has 1,4
        # topology, and every 1,5-RuAAC compound actually has 1,5 topology.
        _pat_14 = Chem.MolFromSmarts("[#6][n]1[n][n][c]([c])[c]1")  # aryl on C4
        _pat_15 = Chem.MolFromSmarts("[#6][n]1[n][n][c][c]1[c]")    # aryl on C5
        _do_geometry_filter = (
            self.regioisomer_label in ("1,4-CuAAC", "1,5-RuAAC")
            and _pat_14 is not None
            and _pat_15 is not None
        )

        for bb_idx, bb in enumerate(self.building_blocks, 1):
            if product_count >= self.config.max_products:
                self.logger.warning(f"Reached product limit ({self.config.max_products}).")
                break

            bb_id = bb['id']
            bb_mol = bb['mol']
            self.logger.info(f"[{bb_idx}/{len(self.building_blocks)}] Processing {bb_id}...")

            reaction_products = self._apply_reaction(self.scaffold_mol, bb_mol)
            if not reaction_products:
                self.logger.warning(f"  ✗ No products from {bb_id}.")
                self.failed_products.append({
                    'building_block_id': bb_id,
                    'failure_reason': 'No products from reaction'
                })
                continue

            bb_product_count = 0
            for prod_idx, prod_tuple in enumerate(reaction_products):
                if product_count >= self.config.max_products:
                    break
                for product_mol in prod_tuple:
                    if product_count >= self.config.max_products:
                        break

                    is_valid, product_smiles, failure_reason =                         self._sanitize_and_validate_product(product_mol, "SCAF-001", bb_id)

                    if not is_valid:
                        self.failed_products.append({
                            'building_block_id': bb_id,
                            'failure_reason': failure_reason
                        })
                        continue

                    # Standardise before InChIKey and descriptors.
                    # Use the RDKit molecule directly rather than round-tripping
                    # through SMILES, because some triazole representations can
                    # become non-parsable after kekulisation/aromaticity handling.
                    try:
                        std_mol = standardize_mol(Chem.Mol(product_mol))
                    except Exception:
                        std_mol = Chem.Mol(product_mol)

                    if std_mol is None:
                        self.failed_products.append({
                            'building_block_id': bb_id,
                            'failure_reason': 'Standardisation produced None'
                        })
                        continue

                    try:
                        Chem.AssignStereochemistry(std_mol, cleanIt=True, force=True)
                    except Exception:
                        pass

                    # Deduplicate by InChIKey (not canonical SMILES)
                    ikey = mol_to_inchikey(std_mol)
                    if ikey is None:
                        ikey = product_smiles  # fallback
                    if ikey in seen_inchikeys:
                        self.logger.debug(f"  ℹ Duplicate InChIKey from {bb_id}; skipping.")
                        continue
                    seen_inchikeys.add(ikey)

                    # Descriptors on standardised molecule
                    try:
                        mol_wt       = Descriptors.MolWt(std_mol)
                        h_donors     = Lipinski.NumHDonors(std_mol)
                        h_acceptors  = Lipinski.NumHAcceptors(std_mol)
                        logp         = Crippen.MolLogP(std_mol)
                        tpsa         = Descriptors.TPSA(std_mol)
                        rot_bonds    = Lipinski.NumRotatableBonds(std_mol)
                        product_inchi = Chem.MolToInchi(std_mol)
                        std_smiles   = Chem.MolToSmiles(std_mol, isomericSmiles=True)
                    except Exception as e:
                        self.logger.warning(f"Descriptor calculation failed for {bb_id}: {e}")
                        continue

                    # Geometry filter: discard products whose triazole ring topology
                    # doesn't match the template's expected regiochemistry.
                    if _do_geometry_filter:
                        _check = Chem.MolFromSmiles(std_smiles)
                        if _check is not None:
                            _is_14 = _check.HasSubstructMatch(_pat_14)
                            _is_15 = _check.HasSubstructMatch(_pat_15)
                            _keep = (
                                (self.regioisomer_label == "1,4-CuAAC" and _is_14) or
                                (self.regioisomer_label == "1,5-RuAAC" and _is_15)
                            )
                            if not _keep:
                                _actual = (
                                    "1,5" if _is_15 else
                                    "1,4" if _is_14 else
                                    "unknown"
                                )
                                self.failed_products.append({
                                    'building_block_id': bb_id,
                                    'failure_reason': (
                                        f"geometry_mismatch: expected "
                                        f"{self.regioisomer_label}, product has "
                                        f"{_actual} topology (wrong alkyne orientation)"
                                    )
                                })
                                continue

                    # Build a regioisomer-unique compound ID.
                    # Each building block can produce more than one product per
                    # reaction run (e.g. both alkyne orientations match [C:4]#[C:5]),
                    # so we need BOTH a regioisomer tag AND a per-product counter to
                    # guarantee globally unique IDs. Without the regioisomer tag, the
                    # CuAAC and RuAAC runs both emit SCAF-001_<bb>_1 and the
                    # lead-score merge produces a many-to-many Cartesian explosion.
                    _regio_suffix_map = {
                        "1,4-CuAAC": "CuAAC",
                        "1,5-RuAAC": "RuAAC",
                    }
                    _regio_tag = _regio_suffix_map.get(self.regioisomer_label)
                    if _regio_tag is not None:
                        # triazole reactions: <bb>_<CuAAC|RuAAC>_<n>
                        compound_id = f"SCAF-001_{bb_id}_{_regio_tag}_{bb_product_count + 1}"
                    else:
                        # non-triazole reactions: preserve original numeric-only format
                        compound_id = f"SCAF-001_{bb_id}_{bb_product_count + 1}"
                    # The regioisomer label comes from the reaction template used to
                    # produce this compound. Each GlycoLibraryGenerator instance
                    # runs exactly one template (CuAAC or RuAAC), so the label is
                    # authoritative at the generator level.
                    regioisomer = self.regioisomer_label

                    self.products.append({
                        'compound_id':        compound_id,
                        'parent_scaffold_id': 'SCAF-001',
                        'building_block_id':  bb_id,
                        'regioisomer':        regioisomer,
                        'product_smiles':     std_smiles,
                        'product_inchi':      product_inchi,
                        'product_inchikey':   ikey,
                        'molecular_weight':   mol_wt,
                        'h_donors':           h_donors,
                        'h_acceptors':        h_acceptors,
                        'logp':               logp,
                        'tpsa':               tpsa,
                        'rotatable_bonds':    rot_bonds,
                        'generation_status':  'success'
                    })
                    product_count += 1
                    bb_product_count += 1

            self.logger.info(f"  [OK] Generated {bb_product_count} product(s) from {bb_id}")

        self.logger.info("=" * 70)
        self.logger.info(f"LIBRARY GENERATION COMPLETE")
        self.logger.info(f"Total products: {product_count}")
        self.logger.info(f"Unique products (by InChIKey): {len(seen_inchikeys)}")
        self.logger.info(f"Failed products: {len(self.failed_products)}")
        self.logger.info("=" * 70)

        return pd.DataFrame(self.products)

    def export_library(self, df: pd.DataFrame, filename: str = "generated_library.csv") -> Path:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.output_dir / filename
        df.to_csv(output_path, index=False)
        self.logger.info(f"[OK] Library exported to: {output_path.absolute()}")
        return output_path

    def export_failed_products(self, filename: str = "failed_products.csv") -> Path:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.output_dir / filename
        if self.failed_products:
            pd.DataFrame(self.failed_products).to_csv(output_path, index=False)
            self.logger.info(f"[OK] Dead-letter log exported to: {output_path.absolute()}")
        else:
            self.logger.info("No failed products to export.")
        return output_path


def generate_triazole_library(
    scaffold_smiles: str,
    building_blocks: List[Dict[str, str]],
    config: LibraryConfig = None,
    logger_instance: logging.Logger = None
) -> pd.DataFrame:
    """
    Convenience function: generate both CuAAC (1,4) and RuAAC (1,5) triazole
    products and return a single combined, tagged DataFrame.

    This is the recommended entry point for triazole library generation,
    replacing direct use of TRIAZOLE_FORMATION which was the source of the
    regiochemistry labelling bug.
    """
    config = config or LibraryConfig()
    frames = []
    for rxn_def in ReactionSMARTS.TRIAZOLE_BOTH:
        gen = GlycoLibraryGenerator(
            scaffold_smiles=scaffold_smiles,
            building_blocks=building_blocks,
            reaction_smarts=rxn_def["smarts"],
            regioisomer_label=rxn_def["regioisomer"],
            config=config,
            logger_instance=logger_instance
        )
        df = gen.generate_library()
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


if __name__ == "__main__":
    # 1. Define a dummy core scaffold with an azide group
    test_scaffold = "CC(N=[N+]=[N-])C1=CC=CC=C1" 
    
    # 2. Define a couple of test building blocks (alkynes)
    test_bbs = [
        {'id': 'BB-001', 'smiles': 'C#CC1=CC=CC=C1'},  # Phenylacetylene
        {'id': 'BB-002', 'smiles': 'C#CCO'}            # Propargyl alcohol
    ]
    
    # 3. Run the generator
    print("Starting test run...")
    df = generate_triazole_library(
        scaffold_smiles=test_scaffold,
        building_blocks=test_bbs
    )
    
    # 4. Print the results
    print("\n--- TEST RESULTS ---")
    if not df.empty:
        print(df[['compound_id', 'regioisomer', 'product_smiles']])
        
        # ADD THIS LINE TO SAVE:
        df.to_csv("my_generated_library.csv", index=False)
        print("\n[SUCCESS] Library saved to 'my_generated_library.csv'")
    else:
        print("No products generated. Check RDKit valency errors.")
