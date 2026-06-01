"""
Phase 3: Molecular docking against Galectin-3 (PDB: 3ZSJ) via AutoDock Vina
============================================================================

Workflow:
  1. Load Phase 2 lead list (07_lead_scored.csv)
  2. Standardise SMILES, embed 3D conformers, write ligand PDBQT (no Open Babel)
  3. Dock each ligand with Vina (subprocess)
  4. Merge docking affinities with Phase 2 lead scores → 08_docking_results.csv
  5. Optional receptor validation: redock lactose vs crystal pose (RMSD < 2.0 Å)

Requires:
  - AutoDock Vina on PATH (or path set in DockingConfig.vina_executable)
  - Prepared receptor PDBQT (3ZSJ) at paths given in DockingConfig
  - RDKit >= 2022.09.1
"""

from __future__ import annotations

import logging
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

try:
    from .glycolibrary_generator import configure_logging, standardize_mol
except ImportError:  # pragma: no cover
    from glycolibrary_generator import configure_logging, standardize_mol  # type: ignore


# β-D-galactopyranosyl-(1→4)-D-glucopyranose (lactose), CRD ligand in 3ZSJ
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
    lead_csv: Path = field(default_factory=lambda: Path("phase2_output/07_lead_scored.csv"))
    lactose_smiles: str = LACTOSE_SMILES
    rmsd_threshold_angstrom: float = 2.0


def minmax_normalize(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    """Min–max normalise a series to [0, 1]; constant series → 0.5."""
    values = pd.to_numeric(series, errors="coerce")
    mn, mx = values.min(), values.max()
    if pd.isna(mn) or pd.isna(mx) or mx <= mn:
        return pd.Series(0.5, index=series.index, dtype=float)
    norm = (values - mn) / (mx - mn)
    if not higher_is_better:
        norm = 1.0 - norm
    return norm


def assign_autodock_atom_type(atom: Chem.Atom) -> Optional[str]:
    """Map RDKit atom to AutoDock4/PDBQT atom type (heavy atoms only)."""
    symbol = atom.GetSymbol()
    if symbol == "H":
        return None
    if symbol in {"F", "Cl", "Br", "I"}:
        return symbol
    if symbol == "S":
        return "S"
    if symbol == "P":
        return "P"
    if symbol == "O":
        return "OA" if atom.GetTotalNumHs() > 0 else "O"
    if symbol == "N":
        return "NA" if atom.GetTotalNumHs() > 0 else "N"
    if symbol == "C":
        if atom.GetIsAromatic():
            return "A"
        return "C"
    return "C"


def is_rotatable_bond(bond: Chem.Bond, mol: Chem.Mol) -> bool:
    """Heuristic rotatable bond definition aligned with med-chem docking prep."""
    if bond.GetBondType() != Chem.BondType.SINGLE:
        return False
    if bond.IsInRing():
        return False
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    if a1.GetDegree() == 1 or a2.GetDegree() == 1:
        return False
    # Exclude amide C–N
    amide = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide):
        for match in mol.GetSubstructMatches(amide):
            c_idx, n_idx = match[0], match[1]
            if {a1.GetIdx(), a2.GetIdx()} == {c_idx, n_idx}:
                return False
    return True


def get_rotatable_bonds(mol: Chem.Mol) -> List[Tuple[int, int]]:
    """Return list of (begin_idx, end_idx) for rotatable bonds."""
    bonds: List[Tuple[int, int]] = []
    for bond in mol.GetBonds():
        if is_rotatable_bond(bond, mol):
            bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    return bonds


def _fragment_atom_indices(mol: Chem.Mol, break_a: int, break_b: int) -> Tuple[Set[int], Set[int]]:
    """Split molecule into two atom sets by removing bond a–b."""
    rw = Chem.RWMol(mol)
    if rw.GetBondBetweenAtoms(break_a, break_b) is not None:
        rw.RemoveBond(break_a, break_b)
    frags = Chem.GetMolFrags(rw, asMols=False)
    if len(frags) != 2:
        return set(frags[0]) if frags else set(), set()
    a_set, b_set = set(frags[0]), set(frags[1])
    if break_a in a_set:
        return a_set, b_set
    return b_set, a_set


def _pdbqt_atom_line(
    mol: Chem.Mol,
    conf_id: int,
    atom_idx: int,
    serial: int,
) -> str:
    atom = mol.GetAtomWithIdx(atom_idx)
    ad_type = assign_autodock_atom_type(atom)
    if ad_type is None:
        return ""
    pos = mol.GetConformer(conf_id).GetAtomPosition(atom_idx)
    return (
        f"ATOM  {serial:5d}  {ad_type:2s}  LIG A   1    "
        f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00    "
        f"{ad_type:2s}\n"
    )


def _write_pdbqt_branch(
    mol: Chem.Mol,
    conf_id: int,
    atom_indices: Set[int],
    lines: List[str],
    serial_counter: List[int],
    rotatable_bonds: List[Tuple[int, int]],
    depth: int = 0,
) -> None:
    """Recursively emit ROOT / BRANCH / ENDBRANCH blocks for a fragment."""
    remaining = set(atom_indices)
    local_rot = [
        (a, b)
        for a, b in rotatable_bonds
        if a in remaining and b in remaining
    ]

    if not local_rot:
        lines.append("ROOT\n")
        for idx in sorted(remaining):
            serial_counter[0] += 1
            line = _pdbqt_atom_line(mol, conf_id, idx, serial_counter[0])
            if line:
                lines.append(line)
        lines.append("ENDROOT\n")
        return

    # Pick bond separating smallest branch (AutoDock-style)
    best_bond = None
    best_branch_size = len(remaining) + 1
    for a, b in local_rot:
        side_a, side_b = _fragment_atom_indices(mol, a, b)
        branch_a = side_a & remaining
        branch_b = side_b & remaining
        if not branch_a or not branch_b:
            continue
        if len(branch_a) <= len(branch_b):
            root_atoms, branch_atoms, root_atom, branch_atom = (
                remaining - branch_a,
                branch_a,
                b,
                a,
            )
        else:
            root_atoms, branch_atoms, root_atom, branch_atom = (
                remaining - branch_b,
                branch_b,
                a,
                b,
            )
        if 0 < len(branch_atoms) < best_branch_size:
            best_branch_size = len(branch_atoms)
            best_bond = (root_atoms, branch_atoms, root_atom, branch_atom, a, b)

    if best_bond is None:
        lines.append("ROOT\n")
        for idx in sorted(remaining):
            serial_counter[0] += 1
            line = _pdbqt_atom_line(mol, conf_id, idx, serial_counter[0])
            if line:
                lines.append(line)
        lines.append("ENDROOT\n")
        return

    root_atoms, branch_atoms, root_atom, branch_atom, _, _ = best_bond

    lines.append("ROOT\n")
    for idx in sorted(root_atoms):
        serial_counter[0] += 1
        line = _pdbqt_atom_line(mol, conf_id, idx, serial_counter[0])
        if line:
            lines.append(line)
    lines.append("ENDROOT\n")

    lines.append(f"BRANCH {root_atom + 1:4d} {branch_atom + 1:4d}\n")
    _write_pdbqt_branch(
        mol, conf_id, branch_atoms, lines, serial_counter, rotatable_bonds, depth + 1
    )
    lines.append("ENDBRANCH\n")


def mol_to_pdbqt_string(mol: Chem.Mol, conf_id: int = 0, name: str = "LIG") -> str:
    """
    Format an RDKit conformer as an AutoDock/Vina ligand PDBQT string.

    Handles AD4 atom types, coordinates, and rotatable-bond tree (ROOT/BRANCH/TORSDOF).
    """
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformer for PDBQT export")

    rotatable = get_rotatable_bonds(mol)
    lines: List[str] = [f"REMARK  Name = {name}\n"]
    serial_counter = [0]
    all_heavy = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if assign_autodock_atom_type(a) is not None
    }
    _write_pdbqt_branch(mol, conf_id, all_heavy, lines, serial_counter, rotatable)
    lines.append(f"TORSDOF {len(rotatable)}\n")
    return "".join(lines)


def embed_ligand_3d(mol: Chem.Mol, random_seed: int = 42) -> Tuple[Optional[Chem.Mol], str]:
    """
    Embed and MMFF-optimise a 3D conformer.

    Returns (mol_with_conf, status) where status is 'success' or 'prep_failed'.
    """
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
    """
    Extract the largest matching HETATM ligand from a crystal PDB.

    Returns (N×3 coordinates, RDKit molecule from PDB block).
    """
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
            f"No ligand HETATM records found in {pdb_path} "
            f"for residue names: {sorted(resnames)}"
        )

    best_res = max(lines_by_res.keys(), key=lambda r: len(lines_by_res[r]))
    block = "\n".join(lines_by_res[best_res]) + "\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=True, removeHs=True)
    if mol is None:
        raise ValueError(f"RDKit could not parse ligand {best_res} from {pdb_path}")

    conf = mol.GetConformer()
    coords = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=float,
    )
    return coords, mol


def compute_rmsd(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    *,
    align: bool = False,
) -> float:
    """
    RMSD between two N×3 coordinate sets (same atom count).

    Parameters
    ----------
    align
        If True, apply Kabsch superposition before RMSD (for MCS-matched poses).
        If False, compute direct RMSD in the input frame (standard redock check
        when receptor is fixed and poses share the same coordinate system).
    """
    if coords_a.shape != coords_b.shape:
        raise ValueError(
            f"Coordinate shape mismatch: {coords_a.shape} vs {coords_b.shape}"
        )
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
    """MCS-based heavy-atom RMSD after alignment (for redocking validation)."""
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

    rmsd = rdMolAlign.AlignMol(
        mob,
        ref,
        atomMap=list(zip(mob_idx, ref_idx)),
    )
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
        """
        Generate 3D conformers from SMILES; tag failed rows as ``prep_failed``.

        Reads ``product_smiles`` (or ``smiles``) per row; writes ligand PDBQT files
        under ``phase3_output/ligands/``.
        """
        out = df.copy()
        smiles_col = "product_smiles" if "product_smiles" in out.columns else "smiles"
        if smiles_col not in out.columns:
            raise KeyError("Input DataFrame must contain 'product_smiles' or 'smiles'")

        statuses: List[str] = []
        pdbqt_paths: List[Optional[str]] = []

        for _, row in out.iterrows():
            compound_id = str(row.get("compound_id", "unknown"))
            smiles = row[smiles_col]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                statuses.append("prep_failed")
                pdbqt_paths.append(None)
                self.logger.warning(f"[prep] Invalid SMILES for {compound_id}")
                continue

            std = standardize_mol(mol)
            if std is not None:
                mol = std

            embedded, status = embed_ligand_3d(mol)
            if status != "success" or embedded is None:
                statuses.append("prep_failed")
                pdbqt_paths.append(None)
                self.logger.warning(f"[prep] 3D embedding failed for {compound_id}")
                continue

            try:
                pdbqt_text = mol_to_pdbqt_string(embedded, conf_id=0, name=compound_id)
                lig_path = self.config.output_dir / "ligands" / f"{compound_id}.pdbqt"
                lig_path.write_text(pdbqt_text, encoding="utf-8")
                statuses.append("success")
                pdbqt_paths.append(str(lig_path))
                self.logger.info(f"[prep] {compound_id} → {lig_path.name}")
            except Exception as exc:
                statuses.append("prep_failed")
                pdbqt_paths.append(None)
                self.logger.warning(f"[prep] PDBQT failed for {compound_id}: {exc}")

        out["generation_status"] = statuses
        out["ligand_pdbqt"] = pdbqt_paths
        return out

    def _build_vina_command(
        self,
        ligand_pdbqt: Path,
        output_pdbqt: Path,
        log_path: Path,
    ) -> List[str]:
        half = self.config.box_size / 2.0
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
            "--log",
            str(log_path),
        ]

    def run_docking(self, prepared_df: pd.DataFrame) -> pd.DataFrame:
        """
        Run Vina via subprocess; parse best-pose affinity (kcal/mol) per compound.

        Rows with ``generation_status != 'success'`` are marked ``dock_skipped``.
        Failed Vina runs are marked ``dock_failed``.
        """
        if not self.config.receptor_pdbqt.exists():
            raise FileNotFoundError(
                f"Receptor PDBQT not found: {self.config.receptor_pdbqt}. "
                "Prepare 3ZSJ receptor before docking."
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
            log_path = self.config.output_dir / "poses" / f"{compound_id}.log"

            cmd = self._build_vina_command(lig_path, pose_path, log_path)
            self.logger.info(f"[dock] Running Vina for {compound_id}")

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=600,
                )
            except FileNotFoundError:
                self.logger.error(
                    f"Vina executable not found: {self.config.vina_executable}"
                )
                vina_scores.append(None)
                dock_statuses.append("dock_failed")
                pose_paths.append(None)
                continue
            except subprocess.TimeoutExpired:
                self.logger.error(f"Vina timed out for {compound_id}")
                vina_scores.append(None)
                dock_statuses.append("dock_failed")
                pose_paths.append(None)
                continue

            log_text = ""
            if log_path.exists():
                log_text = log_path.read_text(encoding="utf-8", errors="replace")

            affinity = parse_vina_affinity(result.stdout, result.stderr + "\n" + log_text)

            if affinity is None or result.returncode != 0:
                self.logger.warning(
                    f"[dock] Vina failed for {compound_id} "
                    f"(returncode={result.returncode})"
                )
                vina_scores.append(None)
                dock_statuses.append("dock_failed")
                pose_paths.append(None)
                continue

            self.logger.info(f"[dock] {compound_id}: {affinity:.2f} kcal/mol")
            vina_scores.append(affinity)
            dock_statuses.append("dock_success")
            pose_paths.append(str(pose_path) if pose_path.exists() else None)

        out["vina_score"] = vina_scores
        out["docking_status"] = dock_statuses
        out["pose_pdbqt"] = pose_paths
        return out

    def merge_with_lead_scores(
        self,
        docking_df: pd.DataFrame,
        scored_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Join on ``compound_id``; add ``combined_score`` =
        0.4 × norm(vina_score) + 0.6 × norm(lead_score).

        Vina affinities are inverted so more negative (favourable) scores map to
        higher normalised values. Writes ``08_docking_results.csv`` to output_dir.
        """
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
            # More negative Vina affinity = better → higher_is_better after invert
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
        self.logger.info(f"[OK] Wrote merged docking results → {out_path}")
        return merged

    def validate_receptor(self) -> bool:
        """
        Redock lactose into 3ZSJ; assert RMSD < 2.0 Å vs crystal pose.

        Raises
        ------
        AssertionError
            If RMSD >= configured threshold or validation cannot complete.
        FileNotFoundError
            If receptor files are missing.
        """
        if not self.config.receptor_pdbqt.exists():
            raise FileNotFoundError(
                f"Receptor PDBQT required for validation: {self.config.receptor_pdbqt}"
            )
        if not self.config.receptor_pdb.exists():
            raise FileNotFoundError(
                f"Receptor PDB required for crystal ligand: {self.config.receptor_pdb}"
            )

        crystal_coords, crystal_mol = extract_ligand_coords_from_pdb(
            self.config.receptor_pdb
        )

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
        lactose_pdbqt.write_text(
            mol_to_pdbqt_string(embedded, conf_id=0, name="LACTOSE"),
            encoding="utf-8",
        )
        pose_out = lig_dir / "lactose_redock_out.pdbqt"
        log_out = lig_dir / "lactose_redock.log"

        cmd = self._build_vina_command(lactose_pdbqt, pose_out, log_out)
        self.logger.info("[validate] Redocking lactose for receptor validation")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                timeout=600,
            )
        except FileNotFoundError as exc:
            raise AssertionError(
                f"Vina not available for receptor validation: {exc}"
            ) from exc

        log_text = ""
        if log_out.exists():
            log_text = log_out.read_text(encoding="utf-8", errors="replace")

        affinity = parse_vina_affinity(result.stdout, result.stderr + "\n" + log_text)
        if affinity is None or not pose_out.exists():
            raise AssertionError(
                "Lactose redocking failed; cannot compute validation RMSD. "
                f"Vina return code: {result.returncode}"
            )

        self.logger.info(f"[validate] Lactose redock affinity: {affinity:.2f} kcal/mol")

        docked_coords = parse_pdbqt_coords(pose_out)
        n_crystal = crystal_coords.shape[0]
        n_docked = docked_coords.shape[0]

        if n_crystal == n_docked:
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

        self.logger.info(f"[validate] Lactose redock RMSD: {rmsd:.3f} Å")

        if rmsd >= self.config.rmsd_threshold_angstrom:
            raise AssertionError(
                f"Receptor validation failed: lactose redock RMSD {rmsd:.3f} Å "
                f">= {self.config.rmsd_threshold_angstrom:.1f} Å threshold"
            )

        self.logger.info(
            f"[OK] Receptor validation passed (RMSD {rmsd:.3f} Å < "
            f"{self.config.rmsd_threshold_angstrom:.1f} Å)"
        )
        return True


def run_phase3_pipeline(
    config: Optional[DockingConfig] = None,
    validate: bool = True,
) -> pd.DataFrame:
    """
    End-to-end Phase 3: prepare ligands → dock → merge scores → export CSV.

    Parameters
    ----------
    config
        Docking configuration; uses defaults if None.
    validate
        If True, run ``validate_receptor()`` before docking leads.

    Returns
    -------
    pd.DataFrame
        Merged results (also written to ``08_docking_results.csv``).
    """
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
        raise FileNotFoundError(
            f"Lead file not found: {lead_path}. Run Phase 2 first."
        )

    leads = pd.read_csv(lead_path)
    scored = leads.copy()

    docker = VinaDocking(config, logger=logger)

    if validate:
        logger.info("\n--- Receptor validation (lactose redock) ---")
        docker.validate_receptor()

    logger.info("\n--- Ligand preparation ---")
    prepared = docker.prepare_ligands(leads)

    logger.info("\n--- Vina docking ---")
    docked = docker.run_docking(prepared)

    logger.info("\n--- Merge with lead scores ---")
    merged = docker.merge_with_lead_scores(docked, scored)

    n_ok = (merged["docking_status"] == "dock_success").sum()
    logger.info(f"\nDocking complete: {n_ok}/{len(merged)} compounds docked successfully")
    logger.info("=" * 70)
    return merged


if __name__ == "__main__":
    cfg = DockingConfig(
        receptor_pdb=Path("data/docking/3ZSJ.pdb"),
        receptor_pdbqt=Path("data/docking/3ZSJ.pdbqt"),
    )
    run_phase3_pipeline(cfg)
