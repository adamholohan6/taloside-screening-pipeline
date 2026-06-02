"""
Phase 3: Molecular docking against Galectin-3 (PDB: 3ZSJ) via AutoDock Vina
============================================================================

Workflow:
  1. Load Phase 2 lead list (07_lead_scored.csv)
  2. Standardise SMILES, embed 3D conformers, write ligand PDBQT (no Open Babel)
  3. Dock each ligand with Vina (subprocess)
  4. Merge docking affinities with Phase 2 lead scores -> 08_docking_results.csv
  5. Optional receptor validation: redock lactose vs crystal pose (RMSD < 2.0 A)

Requires:
  - AutoDock Vina on PATH (or path set in DockingConfig.vina_executable)
  - Prepared receptor PDBQT (3ZSJ) at paths given in DockingConfig
  - RDKit >= 2022.09.1
"""

from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

try:
    from .library_generator import configure_logging
except ImportError:
    from library_generator import configure_logging  # type: ignore

try:
    from .pdbqt_meeko import mol_to_pdbqt_string_meeko
except ImportError:
    mol_to_pdbqt_string_meeko = None  # type: ignore


try:
    from .descriptor_calculator import validate_smiles
except ImportError:
    from descriptor_calculator import validate_smiles  # type: ignore


def standardize_mol(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Light standardisation: strip salts (largest fragment), neutralise charges."""
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        largest = rdMolStandardize.FragmentParent(mol)
        uncharge = rdMolStandardize.Uncharger()
        std = uncharge.uncharge(largest)
        Chem.SanitizeMol(std)
        return std
    except Exception:
        return None


# beta-D-galactopyranosyl-(1->4)-D-glucopyranose (lactose), CRD ligand in 3ZSJ
LACTOSE_SMILES = (
    "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@H](O)[C@@H]2O)"
    "[C@H](O)[C@@H](O)[C@@H]1O"
)

# HETATM residue names commonly used for lactose / galactosides in 3ZSJ-like structures
LIGAND_RESNAMES = frozenset(
    {"LAC", "LCT", "BGC", "GAL", "GLC", "LNT", "NAG", "BMA", "FUL"}
)

_VINA_AFFINITY_RE = re.compile(
    r"^\s*1\s+(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)",
    re.MULTILINE,
)


@dataclass
class DockingConfig:
    """Configuration for AutoDock Vina docking against Galectin-3 (3ZSJ)."""

    receptor_pdb: Path
    receptor_pdbqt: Path
    center_x: float = 10.0
    center_y: float = 15.0
    center_z: float = 5.0
    box_size: float = 20.0
    exhaustiveness: int = 8
    n_poses: int = 3
    output_dir: Path = field(default_factory=lambda: Path("phase3_output"))
    vina_executable: str = "vina"
    lead_csv: Path = field(
        default_factory=lambda: Path("phase2_output/07_lead_scored.csv")
    )
    lactose_smiles: str = LACTOSE_SMILES
    rmsd_threshold_angstrom: float = 2.0


def minmax_normalize(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    """Min-max normalise a series to [0, 1]; constant series -> 0.5."""
    values = pd.to_numeric(series, errors="coerce")
    mn, mx = values.min(), values.max()
    if pd.isna(mn) or pd.isna(mx) or mx <= mn:
        return pd.Series(0.5, index=series.index, dtype=float)
    norm = (values - mn) / (mx - mn)
    if not higher_is_better:
        norm = 1.0 - norm
    return norm


def mol_to_pdbqt_string(mol: Chem.Mol, conf_id: int = 0, name: str = "LIG") -> str:
    """Generate a Vina/AD4-compatible ligand PDBQT string.

    We emit ATOM records ending with the atom Type token placed in columns 77-78.
    """

    def get_ad_type(atom: Chem.Atom) -> Optional[str]:
        sym = atom.GetSymbol()
        if sym == "H":
            return None
        if sym == "C":
            return "C"
        if sym == "N":
            return "NA"
        if sym == "O":
            return "OA"
        if sym in {"S", "P"}:
            return sym
        if sym in {"F", "Cl", "Br", "I"}:
            return sym
        return "C"

    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformer for PDBQT export")

    conf = mol.GetConformer(conf_id)

    lines: List[str] = [f"REMARK  Name = {name}\
"]
    serial = 0

    res_name = "LIG"
    chain_id = "A"
    res_seq = 1

    for atom in mol.GetAtoms():
        ad_type = get_ad_type(atom)
        if ad_type is None:
            continue

        serial += 1
        pos = conf.GetAtomPosition(atom.GetIdx())

        atom_name = atom.GetSymbol().upper().ljust(4)[:4]

        # Build the fixed-width ATOM line and force the atom type into columns 77-78
        atom_line = (
            f"ATOM  {serial:5d} {atom_name:>4s} {res_name:<3s} {chain_id}{res_seq:4d}    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00"
        )
        # pad to 76 chars so the next two chars occupy columns 77-78
        atom_line = atom_line.ljust(76) + f"{ad_type:>2s}\
"
        lines.append(atom_line)

    lines.append("TORSDOF 0\
")
    return "".join(lines)

def embed_ligand_3d(mol: Chem.Mol, random_seed: int = 42) -> Tuple[Optional[Chem.Mol], str]:
    """Embed and MMFF-optimise a 3D conformer."""
    try:
        work = Chem.Mol(mol)
        work = Chem.AddHs(work)
        params = AllChem.ETKDGv3()
        params.randomSeed = random_seed
        params.useSmallRingTorsions = True
        code = AllChem.EmbedMolecule(work, params)
        if code != 0:
            code = AllChem.EmbedMolecule(work, randomSeed=random_seed)
        if code != 0:
            return None, "prep_failed"

        try:
            AllChem.MMFFOptimizeMolecule(work, maxIters=500)
        except Exception:
            pass

        return work, "success"
    except Exception:
        return None, "prep_failed"


def parse_vina_affinity(stdout: str, stderr: str = "") -> Optional[float]:
    """Extract best (mode 1) binding affinity in kcal/mol from Vina output."""
    text = stdout + "\n" + stderr
    match = _VINA_AFFINITY_RE.search(text)
    if match:
        return float(match.group(1))
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                return float(parts[1])
            except ValueError:
                continue
    return None




def parse_pdbqt_coords(pdbqt_path: Path) -> np.ndarray:
    """Parse heavy-atom coordinates from a PDBQT file."""
    coords: List[List[float]] = []
    for line in pdbqt_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue
    if not coords:
        raise ValueError(f"No coordinates parsed from {pdbqt_path}")
    return np.array(coords, dtype=float)


def extract_ligand_coords_from_pdb(
    pdb_path: Path,
    resnames: Optional[Set[str]] = None,
) -> Tuple[np.ndarray, Chem.Mol]:
    """Extract the largest matching HETATM ligand from a crystal PDB."""
    resnames = resnames or set(LIGAND_RESNAMES)
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {pdb_path}")

    lines_by_res: Dict[str, List[str]] = {}
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("HETATM") or line.startswith("ATOM")):
            continue
        if len(line) < 22:
            continue
        resname = line[17:20].strip()
        if resname in resnames:
            lines_by_res.setdefault(resname, []).append(line)

    if not lines_by_res:
        raise ValueError(
            f"No ligand HETATM records found in {pdb_path} for residue names: {sorted(resnames)}"
        )

    best_res = max(lines_by_res.keys(), key=lambda r: len(lines_by_res[r]))
    block = "
".join(lines_by_res[best_res]) + "
"
    mol = Chem.MolFromPDBBlock(block, sanitize=True, removeHs=True)
    if mol is None:
        raise ValueError(f"RDKit could not parse ligand {best_res} from {pdb_path}")

    conf = mol.GetConformer()
    coords = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=float,
    )
    return coords, mol


def compute_rmsd(coords_a: np.ndarray, coords_b: np.ndarray, *, align: bool = False) -> float:
    """RMSD between two N x 3 coordinate sets (same atom count)."""
    if coords_a.shape != coords_b.shape:
        raise ValueError(f"Coordinate shape mismatch: {coords_a.shape} vs {coords_b.shape}")
    if coords_a.shape[0] < 1:
        raise ValueError("Need at least 1 atom for RMSD")

    if not align:
        diff = coords_a - coords_b
        return float(np.sqrt((diff * diff).sum() / coords_a.shape[0]))

    if coords_a.shape[0] < 3:
        raise ValueError("Kabsch alignment requires at least 3 atoms")

    a = coords_a - coords_a.mean(axis=0)
    b = coords_b - coords_b.mean(axis=0)
    cov = a.T @ b
    u, _, vt = np.linalg.svd(cov)
    rot = vt.T @ u.T
    if np.linalg.det(rot) < 0:
        vt[-1, :] *= -1
        rot = vt.T @ u.T
    b_aligned = b @ rot
    diff = a - b_aligned
    return float(np.sqrt((diff * diff).sum() / coords_a.shape[0]))


def align_and_rmsd(ref_mol: Chem.Mol, mob_mol: Chem.Mol) -> float:
    """MCS-based heavy-atom RMSD after alignment."""
    ref = Chem.Mol(ref_mol)
    mob = Chem.Mol(mob_mol)
    if ref.GetNumConformers() == 0 or mob.GetNumConformers() == 0:
        raise ValueError("Both molecules require 3D coordinates for RMSD")

    match = mob.GetSubstructMatch(ref)
    if not match:
        match = ref.GetSubstructMatch(mob)
        if not match:
            raise ValueError("No substructure match for RMSD calculation")
        ref_idx = list(range(len(match)))
        mob_idx = list(match)
    else:
        mob_idx = list(match)
        ref_idx = list(range(len(match)))

    rmsd = rdMolAlign.AlignMol(mob, ref, atomMap=list(zip(mob_idx, ref_idx)))
    return float(rmsd)


class VinaDocking:
    """AutoDock Vina docking pipeline for Phase 2 lead compounds."""

    def __init__(self, config: DockingConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or configure_logging(
            log_file=Path(config.output_dir) / "phase3_docking.log"
        )
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "ligands").mkdir(exist_ok=True)
        (self.config.output_dir / "poses").mkdir(exist_ok=True)

    def prepare_ligands(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        smiles_col = "product_smiles" if "product_smiles" in out.columns else "smiles"
        if smiles_col not in out.columns:
            raise KeyError("Input DataFrame must contain 'product_smiles' or 'smiles'")

        statuses: List[str] = []
        pdbqt_paths: List[Optional[str]] = []

        lig_dir = self.config.output_dir / "ligands"
        lig_dir.mkdir(parents=True, exist_ok=True)

        for _, row in out.iterrows():
            compound_id = str(row.get("compound_id", "unknown"))
            smiles = str(row[smiles_col])

            # Step 1: Attempt standard parsing
            mol = Chem.MolFromSmiles(smiles)

            # Step 2: Fallback with two-pass sanitization for triazole valence errors
            if mol is None:
                self.logger.error("=" * 80)
                self.logger.error(f"FAILED PARSE: {compound_id}")
                self.logger.error(f"SMILES: {smiles}")

                try:
                    raw = Chem.MolFromSmiles(smiles, sanitize=False)

                    if raw is not None:
                        self.logger.error("ATOM TABLE (before sanitization):")

                        for atom in raw.GetAtoms():
                            self.logger.error(
                                f"idx={atom.GetIdx():2d} "
                                f"symbol={atom.GetSymbol():2s} "
                                f"charge={atom.GetFormalCharge():2d} "
                                f"degree={atom.GetDegree():2d} "
                                f"explicitHs={atom.GetNumExplicitHs():2d} "
                                f"aromatic={atom.GetIsAromatic()}"
                            )

                        # Pass 1: Sanitize without strict property/valence validation
                        mask = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                        Chem.SanitizeMol(raw, sanitizeOps=mask)

                        # Graph repair: strip hypervalent explicit Hs from neutral ring nitrogens
                        repairs_made = 0
                        for atom in raw.GetAtoms():
                            if atom.GetAtomicNum() == 7:  # Nitrogen
                                if atom.GetExplicitValence() > 3 and atom.GetFormalCharge() == 0:
                                    if atom.GetNumExplicitHs() > 0:
                                        self.logger.error(
                                            f"Repairing N idx={atom.GetIdx()}: "
                                            f"explicitHs {atom.GetNumExplicitHs()} -> 0"
                                        )
                                        atom.SetNumExplicitHs(0)
                                        repairs_made += 1

                        self.logger.error(f"Repairs made: {repairs_made}")

                        # Pass 2: Re-verify properties strictly
                        raw.UpdatePropertyCache(strict=True)
                        Chem.SanitizeMol(raw)

                        self.logger.error("ATOM TABLE (after sanitization):")
                        for atom in raw.GetAtoms():
                            self.logger.error(
                                f"idx={atom.GetIdx():2d} "
                                f"symbol={atom.GetSymbol():2s} "
                                f"charge={atom.GetFormalCharge():2d} "
                                f"degree={atom.GetDegree():2d} "
                                f"explicitHs={atom.GetNumExplicitHs():2d} "
                                f"aromatic={atom.GetIsAromatic()}"
                            )

                        mol = raw
                        self.logger.error("SUCCESS: Fallback sanitization worked")
                    else:
                        self.logger.error("FAILED: MolFromSmiles(sanitize=False) returned None")
                except Exception as exc:
                    self.logger.error(f"Fallback sanitization failed: {exc}")
                    import traceback
                    self.logger.error(traceback.format_exc())

            if mol is None:
                statuses.append("prep_failed")
                pdbqt_paths.append(None)
                continue

            std = standardize_mol(mol)
            if std is not None:
                mol = std

            embedded, status = embed_ligand_3d(mol)
            if status != "success" or embedded is None:
                statuses.append("prep_failed")
                pdbqt_paths.append(None)
                continue

            try:
                # Always use manual PDBQT writer (Meeko disabled)
                pdbqt_text = mol_to_pdbqt_string(embedded, conf_id=0, name=compound_id)
                lig_path = lig_dir / f"{compound_id}.pdbqt"
                lig_path.write_text(pdbqt_text, encoding="utf-8")

                statuses.append("success")
                pdbqt_paths.append(str(lig_path))
            except Exception:
                statuses.append("prep_failed")
                pdbqt_paths.append(None)

        out["generation_status"] = statuses
        out["ligand_pdbqt"] = pdbqt_paths
        return out

    def _build_vina_command(self, ligand_pdbqt: Path, output_pdbqt: Path) -> List[str]:
        return [
            self.config.vina_executable,
            "--receptor",
            str(self.config.receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--center_x",
            str(self.config.center_x),
            "--center_y",
            str(self.config.center_y),
            "--center_z",
            str(self.config.center_z),
            "--size_x",
            str(self.config.box_size),
            "--size_y",
            str(self.config.box_size),
            "--size_z",
            str(self.config.box_size),
            "--exhaustiveness",
            str(self.config.exhaustiveness),
            "--num_modes",
            str(self.config.n_poses),
            "--out",
            str(output_pdbqt),
        ]

    def run_docking(self, prepared_df: pd.DataFrame) -> pd.DataFrame:
        if not self.config.receptor_pdbqt.exists():
            raise FileNotFoundError(
                f"Receptor PDBQT not found: {self.config.receptor_pdbqt}. Prepare 3ZSJ receptor before docking."
            )

        out = prepared_df.copy()
        vina_scores: List[Optional[float]] = []
        dock_statuses: List[str] = []
        pose_paths: List[Optional[str]] = []

        for _, row in out.iterrows():
            compound_id = str(row.get("compound_id", "unknown"))
            if row.get("generation_status") != "success" or not row.get("ligand_pdbqt"):
                vina_scores.append(None)
                dock_statuses.append("dock_skipped")
                pose_paths.append(None)
                continue

            lig_path = Path(row["ligand_pdbqt"])
            pose_path = self.config.output_dir / "poses" / f"{compound_id}_out.pdbqt"

            cmd = self._build_vina_command(lig_path, pose_path)

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=600,
                )
            except (FileNotFoundError, subprocess.TimeoutExpired):
                vina_scores.append(None)
                dock_statuses.append("dock_failed")
                pose_paths.append(None)
                continue

            affinity = parse_vina_affinity(result.stdout, result.stderr)
            if affinity is None or result.returncode != 0:
                vina_scores.append(None)
                dock_statuses.append("dock_failed")
                pose_paths.append(None)
                continue

            vina_scores.append(affinity)
            dock_statuses.append("dock_success")
            pose_paths.append(str(pose_path) if pose_path.exists() else None)

        out["vina_score"] = vina_scores
        out["docking_status"] = dock_statuses
        out["pose_pdbqt"] = pose_paths
        return out

    def merge_with_lead_scores(self, docking_df: pd.DataFrame, scored_df: pd.DataFrame) -> pd.DataFrame:
        merged = docking_df.merge(
            scored_df[["compound_id", "lead_score"]],
            on="compound_id",
            how="left",
            suffixes=("", "_phase2"),
        )

        docked = merged["vina_score"].notna()
        norm_vina = pd.Series(np.nan, index=merged.index, dtype=float)
        norm_lead = pd.Series(np.nan, index=merged.index, dtype=float)

        if docked.any():
            norm_vina.loc[docked] = minmax_normalize(
                merged.loc[docked, "vina_score"],
                higher_is_better=False,
            )

        has_lead = merged["lead_score"].notna()
        if has_lead.any():
            norm_lead.loc[has_lead] = minmax_normalize(
                merged.loc[has_lead, "lead_score"],
                higher_is_better=True,
            )

        merged["norm_vina_score"] = norm_vina
        merged["norm_lead_score"] = norm_lead
        merged["combined_score"] = np.where(
            docked & has_lead,
            0.4 * norm_vina + 0.6 * norm_lead,
            np.nan,
        )

        out_path = self.config.output_dir / "08_docking_results.csv"
        merged.sort_values(
            ["combined_score", "vina_score"],
            ascending=[False, True],
            na_position="last",
        ).to_csv(out_path, index=False)
        return merged

    def validate_receptor(self) -> bool:
        if not self.config.receptor_pdbqt.exists():
            raise FileNotFoundError(
                f"Receptor PDBQT required for validation: {self.config.receptor_pdbqt}"
            )
        if not self.config.receptor_pdb.exists():
            raise FileNotFoundError(
                f"Receptor PDB required for crystal ligand: {self.config.receptor_pdb}"
            )

        crystal_coords, crystal_mol = extract_ligand_coords_from_pdb(self.config.receptor_pdb)

        lactose = Chem.MolFromSmiles(self.config.lactose_smiles)
        if lactose is None:
            raise ValueError("Invalid lactose SMILES in DockingConfig")

        std = standardize_mol(lactose) or lactose
        embedded, status = embed_ligand_3d(std)
        if status != "success" or embedded is None:
            raise AssertionError("Lactose 3D preparation failed for receptor validation")

        lig_dir = self.config.output_dir / "validation"
        lig_dir.mkdir(parents=True, exist_ok=True)

        lactose_pdbqt = lig_dir / "lactose_redock.pdbqt"
        if mol_to_pdbqt_string_meeko is None:
            raise RuntimeError("Meeko PDBQT export not available: cannot import mol_to_pdbqt_string_meeko")

        lactose_pdbqt.write_text(
    mol_to_pdbqt_string(embedded, conf_id=0, name="lactose"),
    encoding="utf-8",
)



        pose_out = lig_dir / "lactose_redock_out.pdbqt"

        cmd = self._build_vina_command(lactose_pdbqt, pose_out)
        result = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=600)

        if result.returncode != 0 or not pose_out.exists():
            raise AssertionError(
                "Lactose redocking failed; cannot compute validation RMSD. "
                f"Vina return code: {result.returncode}
STDOUT:
{result.stdout}
STDERR:
{result.stderr}"
            )

        docked_coords = parse_pdbqt_coords(pose_out)

        if crystal_coords.shape[0] == docked_coords.shape[0]:
            rmsd = compute_rmsd(crystal_coords, docked_coords, align=False)
        else:
            docked_mol = Chem.MolFromPDBBlock(
                pose_out.read_text(encoding="utf-8", errors="replace"),
                sanitize=True,
                removeHs=True,
            )
            if docked_mol is None or docked_mol.GetNumConformers() == 0:
                raise AssertionError("Could not parse docked lactose pose for RMSD")
            rmsd = align_and_rmsd(crystal_mol, docked_mol)

        if rmsd >= self.config.rmsd_threshold_angstrom:
            raise AssertionError(
                f"Receptor validation failed: lactose redock RMSD {rmsd:.3f} A >= {self.config.rmsd_threshold_angstrom:.1f} A threshold"
            )

        return True


def run_phase3_pipeline(config: Optional[DockingConfig] = None, validate: bool = True) -> pd.DataFrame:
    if config is None:
        config = DockingConfig(
            receptor_pdb=Path("data/docking/3ZSJ.pdb"),
            receptor_pdbqt=Path("data/docking/3ZSJ.pdbqt"),
        )

    logger = configure_logging(log_file=config.output_dir / "phase3_docking.log")
    logger.info("=" * 70)
    logger.info("PHASE 3: GALECTIN-3 DOCKING (3ZSJ)")
    logger.info("=" * 70)

    lead_path = Path(config.lead_csv)
    if not lead_path.exists():
        raise FileNotFoundError(f"Lead file not found: {lead_path}. Run Phase 2 first.")

    leads = pd.read_csv(lead_path)
    scored = leads.copy()

    docker = VinaDocking(config, logger=logger)

    if validate:
        docker.validate_receptor()

    prepared = docker.prepare_ligands(leads)
    docked = docker.run_docking(prepared)
    merged = docker.merge_with_lead_scores(docked, scored)
    return merged


if __name__ == "__main__":
    cfg = DockingConfig(
        receptor_pdb=Path("data/docking/3ZSJ.pdb"),
        receptor_pdbqt=Path("data/docking/3ZSJ.pdbqt"),
        vina_executable=r"C:\\Users\\adamh\\Docking - Copy\\vina.exe",
    )
    run_phase3_pipeline(cfg)

