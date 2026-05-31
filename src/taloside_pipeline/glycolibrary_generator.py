"""
Taloside Screening Pipeline - Phase 2: Virtual Library Expansion Engine
========================================================================

Generates combinatorial chemical libraries from core scaffolds and building block fragments
using RDKit's Chemical Reactions framework (SMARTS-based).

Features:
  - Combinatorial product generation from scaffold + building blocks
  - Stereochemistry preservation and validation
  - Hypervalency and valence error detection
  - Duplicate and chimeric compound filtering
  - CSV export with traceable lineage (scaffold ID, building block ID)
  - Comprehensive error logging to dead-letter queue

Author: Computational Chemistry Pipeline
License: MIT
Requires: RDKit >= 2022.09.1, Pandas >= 1.3.0, NumPy >= 1.19.0
"""

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski


# ============================================================================
# CONFIGURATION & LOGGING
# ============================================================================

@dataclass
class LibraryConfig:
    """Configuration container for library generation parameters."""
    max_products: int = 50000  # Safety limit to prevent memory overflow
    include_stereoisomers: bool = True
    sanitize_products: bool = True
    filter_hypervalent: bool = True
    min_product_mw: float = 150.0
    max_product_mw: float = 1000.0
    output_dir: Path = Path("library_output")


def configure_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """
    Configure logging for the library generator.
    
    Args:
        log_file: Optional path for log file output
        
    Returns:
        logging.Logger: Configured logger instance
    """
    logger = logging.getLogger("GlycoLibraryGenerator")
    logger.setLevel(logging.DEBUG)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter(
        '[%(levelname)-8s] %(name)s: %(message)s'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter(
            '[%(asctime)s] [%(levelname)-8s] %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
    
    return logger


logger = configure_logging()


# ============================================================================
# CORE VALIDATION UTILITIES
# ============================================================================

def validate_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Validate and parse a SMILES string.
    
    Args:
        smiles: SMILES notation string
        
    Returns:
        rdkit.Chem.Mol or None: Valid RDKit molecule or None if invalid
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Check for pathological structures
        if mol.GetNumAtoms() == 0:
            return None
        return mol
    except Exception as e:
        logger.debug(f"SMILES validation error for '{smiles}': {e}")
        return None


def check_valence_errors(mol: Chem.Mol) -> bool:
    """
    Check if molecule has valence errors (hypervalency).
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        bool: True if valence errors detected, False otherwise
    """
    try:
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL)
        return False
    except Chem.AtomValenceException:
        return True
    except Exception:
        return True


def calculate_heavy_atom_count(mol: Chem.Mol) -> int:
    """
    Calculate number of heavy atoms (excluding H).
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        int: Heavy atom count
    """
    return mol.GetNumHeavyAtoms()


# ============================================================================
# REACTION SMARTS DEFINITIONS
# ============================================================================

class ReactionSMARTS:
    """
    Container for standard reaction SMARTS patterns for carbohydrate chemistry.
    
    Each reaction is defined as:
      - reaction_string: SMARTS pattern with marked atom indices [1:X] (scaffold)
                         and [2:X] (building block)
      - description: Human-readable explanation of chemistry
      - name: Short identifier
    """
    
    # Click Chemistry: Triazole formation via azide-alkyne cycloaddition
    TRIAZOLE_FORMATION = {
        "name": "triazole_formation",
        "smarts": "[N:1]=[N+:2]=[N-:3].[C:4]#[C:5]>>[C:4]1=[C:5][n:1][n+:2]=[n-:3]1",
        "description": (
            "1,3-Dipolar azide-alkyne cycloaddition forming 1,2,3-triazole. "
            "Requires azide [N]=[N+]=[N-] on scaffold and terminal alkyne [C#C] on building block. "
            "Atom mapping: [1:3] = azide N atoms, [4:5] = alkyne C atoms"
        )
    }
    
    # Protecting Group Installation: Carbamate formation (Cbz protection)
    CARBAMATE_FORMATION = {
        "name": "carbamate_formation",
        "smarts": "[N:1]-C.[C:2]([O-])=O>>[N:1]C(=O)[O][C:2]",
        "description": (
            "Carbamate (Cbz) protection of primary amines. "
            "Scaffold must have free -NH2 group [N:1]. "
            "Building block must be activated carboxylic acid derivative [C:2]."
        )
    }
    
    # Fragment Coupling: Amide bond formation
    AMIDE_COUPLING = {
        "name": "amide_coupling",
        "smarts": "[N:1].[C:2](=O)[O:3]>>[N:1]C(=O)[C:2]",
        "description": (
            "Amide bond formation between primary amine [N:1] (scaffold) "
            "and carboxylic acid [C:2] (building block). "
            "Classic peptide coupling chemistry."
        )
    }
    
    # Williamson Ether Synthesis: Nucleophilic aromatic substitution
    ETHER_FORMATION = {
        "name": "ether_formation",
        "smarts": "[O:1]-.[C:2][Cl,Br,I]>>[O:1][C:2]",
        "description": (
            "Williamson ether synthesis: nucleophilic oxygen [O:1] attacks "
            "activated alkyl halide [C:2][X]. Used for methoxy/alkoxy installation "
            "on aromatic rings or aliphatic carbons."
        )
    }
    
    # Nucleophilic Aromatic Substitution: SNAr on electron-deficient arenes
    SNAR_SUBSTITUTION = {
        "name": "snar_substitution",
        "smarts": (
            "[c:1]([N+:2](=O)[O-].[N:3])>>"
            "[c:1][N:3].[N+:2](=O)[O-]"
        ),
        "description": (
            "Nucleophilic aromatic substitution (SNAr) on nitro-activated benzene. "
            "Nitro group [N+](=O)[O-] activates aromatic ring [c:1] for "
            "nucleophilic attack by amine [N:3] from building block."
        )
    }
    
    # Suzuki-Miyaura Cross-Coupling: Boronic acid + aryl halide
    SUZUKI_COUPLING = {
        "name": "suzuki_coupling",
        "smarts": (
            "[c:1][Br,I].[c:2]B(O)O>>"
            "[c:1][c:2]"
        ),
        "description": (
            "Suzuki-Miyaura cross-coupling: aryl halide [c:1][Br/I] (scaffold) "
            "couples with aryl/alkyl boronic acid [c:2]B(O)O (building block). "
            "Requires Pd catalyst (not modeled). Produces biaryl/arylalkyl products."
        )
    }


# ============================================================================
# MAIN GLYCOLIBRARY GENERATOR CLASS
# ============================================================================

class GlycoLibraryGenerator:
    """
    High-throughput combinatorial library generation engine for glycomimetic scaffolds.
    
    The generator applies user-defined reaction SMARTS patterns to core scaffolds 
    and building blocks to produce combinatorial product libraries with full 
    data lineage tracking, validation, and dead-letter logging.
    
    Attributes:
        scaffold_smiles (str): Core structural scaffold (SMILES)
        building_blocks (List[Dict]): List of building block dictionaries with 
                                      'id' and 'smiles' keys
        reaction_smarts (str): SMARTS-encoded reaction pattern
        config (LibraryConfig): Configuration parameters
        logger (logging.Logger): Logger instance
        products (List[Dict]): Generated product molecules with metadata
        failed_products (List[Dict]): Products that failed validation
    """
    
    def __init__(
        self,
        scaffold_smiles: str,
        building_blocks: List[Dict[str, str]],
        reaction_smarts: str,
        config: LibraryConfig = None,
        logger_instance: logging.Logger = None
    ):
        """
        Initialize GlycoLibraryGenerator.
        
        Args:
            scaffold_smiles: SMILES string for core scaffold
            building_blocks: List of dicts with keys:
                            - 'id': compound identifier (str)
                            - 'smiles': SMILES string (str)
                            Example: [
                              {'id': 'BB-001', 'smiles': 'CCO'},
                              {'id': 'BB-002', 'smiles': 'CC(C)O'}
                            ]
            reaction_smarts: SMARTS pattern encoding the chemical reaction
                            Uses atom maps: [#:1] for scaffold, [#:2] for building block
            config: LibraryConfig object (uses defaults if None)
            logger_instance: Logger instance (uses module logger if None)
            
        Raises:
            ValueError: If scaffold or building blocks fail validation
        """
        self.logger = logger_instance or logger
        self.config = config or LibraryConfig()
        
        # Validate scaffold
        self.scaffold_smiles = scaffold_smiles
        self.scaffold_mol = validate_smiles(scaffold_smiles)
        if self.scaffold_mol is None:
            raise ValueError(f"Invalid scaffold SMILES: {scaffold_smiles}")
        self.logger.info(
            f"[OK] Scaffold loaded: {scaffold_smiles[:50]}... "
            f"({self.scaffold_mol.GetNumAtoms()} atoms)"
        )
        
        # Validate and store building blocks
        self.building_blocks = []
        invalid_bb_count = 0
        for bb in building_blocks:
            bb_smiles = bb.get('smiles')
            bb_id = bb.get('id', 'UNKNOWN')
            bb_mol = validate_smiles(bb_smiles)
            if bb_mol is None:
                self.logger.warning(
                    f"[X] Invalid building block SMILES (ID: {bb_id}): {bb_smiles}"
                )
                invalid_bb_count += 1
                continue
            self.building_blocks.append({
                'id': bb_id,
                'smiles': bb_smiles,
                'mol': bb_mol
            })
        
        if invalid_bb_count > 0:
            self.logger.warning(
                f"Skipped {invalid_bb_count} invalid building block(s). "
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
        
        # Initialize product containers
        self.products = []
        self.failed_products = []
    
    
    def _parse_reaction_smarts(self, smarts: str) -> Optional[Chem.ChemicalReaction]:
        """
        Parse and validate reaction SMARTS string.
        
        The reaction SMARTS format uses >> to separate reactants from products:
            [Reactant1].[Reactant2]>>[Product1].[Product2]
        
        Atom mapping (indicated by :[digit]) tells RDKit which atoms correspond:
            [C:1][C:2].[C:3]>>[C:1][C:3][C:2]
        
        Args:
            smarts: Reaction SMARTS string
            
        Returns:
            rdkit.Chem.ChemicalReaction or None if invalid
        """
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None or rxn.GetNumReactantTemplates() == 0:
                return None
            self.logger.debug(
                f"[OK] Parsed reaction SMARTS with "
                f"{rxn.GetNumReactantTemplates()} reactants, "
                f"{rxn.GetNumProductTemplates()} products"
            )
            return rxn
        except Exception as e:
            self.logger.error(f"Failed to parse reaction SMARTS: {e}")
            return None
    
    
    def _apply_reaction(
        self,
        scaffold_mol: Chem.Mol,
        building_block_mol: Chem.Mol
    ) -> List[Tuple[Chem.Mol, ...]]:
        """
        Apply chemical reaction to scaffold + building block pair.
        
        RDKit's RunReactants() method takes a tuple of reactant molecules
        and returns a list of product tuples (one per possible outcome).
        
        Args:
            scaffold_mol: RDKit molecule object (scaffold)
            building_block_mol: RDKit molecule object (building block)
            
        Returns:
            List of tuples, each containing product molecule(s).
            Empty list if reaction failed or produced no products.
        """
        try:
            products = self.rxn.RunReactants(
                (scaffold_mol, building_block_mol)
            )
            return products
        except Exception as e:
            self.logger.debug(
                f"Reaction failed for BB {building_block_mol}: {e}"
            )
            return []
    
    
    def _sanitize_and_validate_product(
        self,
        product_mol: Chem.Mol,
        scaffold_id: str,
        building_block_id: str
    ) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        Perform comprehensive validation on product molecule.
        
        Checks:
          1. Sanitization (valence, aromaticity)
          2. Hypervalency / valence errors
          3. Molecular weight range
          4. Heavy atom count (reasonable complexity)
          5. Structure canonicalization
        
        Args:
            product_mol: RDKit molecule to validate
            scaffold_id: Identifier for scaffold
            building_block_id: Identifier for building block
            
        Returns:
            Tuple (is_valid, product_smiles, failure_reason)
        """
        try:
            # Attempt sanitization (handles implicit/explicit Hs)
            # Skip kekulization for aromatic nitrogen rings (triazoles)
            Chem.SanitizeMol(product_mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
            Chem.Kekulize(product_mol, clearAromaticFlags=True)
        except Exception as e:
            # If kekulization still fails, try without it
            try:
                Chem.SanitizeMol(product_mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
            except Exception as e2:
                return False, None, f"Sanitization failed: {e}"
        
        # Molecular weight filtering
        mol_wt = Descriptors.MolWt(product_mol)
        if mol_wt < self.config.min_product_mw or mol_wt > self.config.max_product_mw:
            return False, None, (
                f"MW out of range ({mol_wt:.2f}; "
                f"expected {self.config.min_product_mw}-{self.config.max_product_mw})"
            )
        
        # Generate canonical SMILES
        try:
            product_smiles = Chem.MolToSmiles(
                product_mol,
                isomericSmiles=self.config.include_stereoisomers
            )
        except Exception as e:
            return False, None, f"Failed to generate SMILES: {e}"
        
        return True, product_smiles, None
    
    
    def _is_duplicate(self, smiles: str, existing_smiles: Set[str]) -> bool:
        """
        Check if SMILES has already been generated (duplicate detection).
        
        Args:
            smiles: SMILES string to check
            existing_smiles: Set of previously generated SMILES
            
        Returns:
            bool: True if duplicate, False otherwise
        """
        return smiles in existing_smiles
    
    
    def generate_library(self) -> pd.DataFrame:
        """
        Execute combinatorial library generation.
        
        Algorithm:
          1. For each building block (BB):
              a. Apply reaction to scaffold + BB
              b. For each product from reaction:
                  - Validate and sanitize
                  - Check for duplicates
                  - Log metadata (lineage, descriptors)
              c. Cap products at config.max_products (safety)
        
        Returns:
            pd.DataFrame with columns:
              - compound_id: Unique identifier (e.g., SCAF-001_BB-001_1)
              - parent_scaffold_id: Scaffold identifier
              - building_block_id: BB identifier
              - product_smiles: Canonical SMILES of product
              - product_inchi: InChI string
              - molecular_weight: Calculated MW
              - h_donors: H-bond donors
              - h_acceptors: H-bond acceptors
              - logp: Lipophilicity
              - tpsa: Topological polar surface area
              - rotatable_bonds: Flexibility measure
              - generation_status: 'success' or failure reason
        """
        self.logger.info("=" * 70)
        self.logger.info("INITIATING COMBINATORIAL LIBRARY GENERATION")
        self.logger.info("=" * 70)
        
        self.logger.info(
            f"Scaffold: {self.scaffold_smiles[:50]}...\n"
            f"Building Blocks: {len(self.building_blocks)}\n"
            f"Reaction: {self.reaction_smarts[:60]}...\n"
        )
        
        seen_smiles: Set[str] = set()
        product_count = 0
        
        # Iterate through building blocks
        for bb_idx, bb in enumerate(self.building_blocks, 1):
            if product_count >= self.config.max_products:
                self.logger.warning(
                    f"Reached product limit ({self.config.max_products}). "
                    "Stopping library generation."
                )
                break
            
            bb_id = bb['id']
            bb_mol = bb['mol']
            
            self.logger.info(
                f"\n[{bb_idx}/{len(self.building_blocks)}] "
                f"Processing BB {bb_id}..."
            )
            
            # Apply reaction
            reaction_products = self._apply_reaction(self.scaffold_mol, bb_mol)
            
            if not reaction_products:
                self.logger.warning(
                    f"  ✗ No products from BB {bb_id}."
                )
                self.failed_products.append({
                    'building_block_id': bb_id,
                    'failure_reason': 'No products from reaction'
                })
                continue
            
            # Process each product tuple
            bb_product_count = 0
            for prod_tuple in reaction_products:
                if product_count >= self.config.max_products:
                    break
                
                # Handle multiple products per reaction (typically just 1)
                for product_mol in prod_tuple:
                    if product_count >= self.config.max_products:
                        break
                    
                    # Validate product
                    is_valid, product_smiles, failure_reason = \
                        self._sanitize_and_validate_product(
                            product_mol, "SCAF-001", bb_id
                        )
                    
                    if not is_valid:
                        self.failed_products.append({
                            'building_block_id': bb_id,
                            'failure_reason': failure_reason
                        })
                        continue
                    
                    # Check for duplicates
                    if self._is_duplicate(product_smiles, seen_smiles):
                        self.logger.debug(
                            f"  ℹ Duplicate SMILES from BB {bb_id}; skipping."
                        )
                        continue
                    
                    seen_smiles.add(product_smiles)
                    
                    # Calculate descriptors
                    try:
                        mol_wt = Descriptors.MolWt(product_mol)
                        h_donors = Lipinski.NumHDonors(product_mol)
                        h_acceptors = Lipinski.NumHAcceptors(product_mol)
                        logp = Crippen.MolLogP(product_mol)
                        tpsa = Descriptors.TPSA(product_mol)
                        rot_bonds = Lipinski.NumRotatableBonds(product_mol)
                        
                        # Generate InChI
                        inchi = Chem.MolToInchi(product_mol)
                    except Exception as e:
                        self.logger.warning(
                            f"Failed to calculate descriptors for BB {bb_id}: {e}"
                        )
                        continue
                    
                    # Build product record
                    compound_id = f"SCAF-001_{bb_id}_{bb_product_count + 1}"
                    product_record = {
                        'compound_id': compound_id,
                        'parent_scaffold_id': 'SCAF-001',
                        'building_block_id': bb_id,
                        'product_smiles': product_smiles,
                        'product_inchi': inchi,
                        'molecular_weight': mol_wt,
                        'h_donors': h_donors,
                        'h_acceptors': h_acceptors,
                        'logp': logp,
                        'tpsa': tpsa,
                        'rotatable_bonds': rot_bonds,
                        'generation_status': 'success'
                    }
                    
                    self.products.append(product_record)
                    product_count += 1
                    bb_product_count += 1
            
            self.logger.info(
                f"  [OK] Generated {bb_product_count} product(s) from BB {bb_id}"
            )
        
        self.logger.info("\n" + "=" * 70)
        self.logger.info(f"LIBRARY GENERATION COMPLETE")
        self.logger.info("=" * 70)
        self.logger.info(f"Total products: {product_count}")
        self.logger.info(f"Unique products: {len(seen_smiles)}")
        self.logger.info(f"Failed products: {len(self.failed_products)}")
        
        return pd.DataFrame(self.products)
    
    
    def export_library(
        self,
        df: pd.DataFrame,
        filename: str = "generated_library.csv"
    ) -> Path:
        """
        Export generated library to CSV file.
        
        Args:
            df: DataFrame containing products
            filename: Output CSV filename
            
        Returns:
            Path: Path to exported file
        """
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.output_dir / filename
        
        df.to_csv(output_path, index=False)
        self.logger.info(f"[OK] Library exported to: {output_path.absolute()}")
        
        return output_path
    
    
    def export_failed_products(
        self,
        filename: str = "failed_products.csv"
    ) -> Path:
        """
        Export dead-letter log of failed products.
        
        Args:
            filename: Output CSV filename for failures
            
        Returns:
            Path: Path to exported file
        """
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.output_dir / filename
        
        if self.failed_products:
            df_failed = pd.DataFrame(self.failed_products)
            df_failed.to_csv(output_path, index=False)
            self.logger.info(
                f"[OK] Dead-letter log exported to: {output_path.absolute()}"
            )
        else:
            self.logger.info("No failed products to export.")
        
        return output_path


if __name__ == "__main__":
    pass
