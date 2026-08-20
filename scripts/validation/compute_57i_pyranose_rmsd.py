# -*- coding: utf-8 -*-
"""
Reproduce the 57I pyranose-ring docking-validation RMSD reported in the manuscript (§3.4).

The 7RGX crystal structure is superposed onto 3ZSJ on common Ca atoms, the 57I
pyranose ring (C1-C5, O5) is extracted from the aligned crystal ligand, and its
RMSD against the corresponding atoms of the redocked pose is reported both
directly and after Kabsch (SVD) superposition.

Reported values (manuscript §3.4):
    Ca RMSD, 7RGX -> 3ZSJ            0.275 A
    pyranose RMSD, direct            2.17 A
    pyranose RMSD, Kabsch-aligned    0.48 A   <- value quoted in the manuscript

The redocked pose (data/docking/7RGX_57I_docked.pdbqt) is committed, so this
check needs neither AutoDock Vina nor a Phase 3 rerun. The two crystal
structures are public RCSB entries and are not tracked; fetch them first:

    curl -o data/docking/3ZSJ.pdb https://files.rcsb.org/download/3ZSJ.pdb
    curl -o data/docking/7RGX.pdb https://files.rcsb.org/download/7RGX.pdb

Usage:
    python scripts/validation/compute_57i_pyranose_rmsd.py

Exits 0 if the computed Kabsch RMSD matches the manuscript value, 1 otherwise.
"""
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer

PROJ = Path(__file__).resolve().parents[2]
PDB_3ZSJ = PROJ / "data" / "docking" / "3ZSJ.pdb"
PDB_7RGX = PROJ / "data" / "docking" / "7RGX.pdb"
PDBQT_DOCKED = PROJ / "data" / "docking" / "7RGX_57I_docked.pdbqt"

# Pyranose ring atoms, in the order used for the crystal/docked correspondence.
RING_ATOMS = ["C1", "C2", "C3", "C4", "C5", "O5"]

# Manuscript §3.4 value and the tolerance within which we consider it reproduced.
MANUSCRIPT_RMSD = 0.48
TOLERANCE = 0.05


def kabsch_rmsd(P, Q):
    """RMSD between two equal-length point sets after optimal superposition.

    P is the reference (n x 3), Q the mobile set. Both are centred, the optimal
    rotation is recovered by SVD of the covariance matrix, and a reflection
    guard keeps the result a proper rotation.
    """
    assert P.shape == Q.shape

    P_center = P - P.mean(axis=0)
    Q_center = Q - Q.mean(axis=0)

    H = np.dot(Q_center.T, P_center)
    U, S, VT = np.linalg.svd(H)
    R = np.dot(VT.T, U.T)

    if np.linalg.det(R) < 0:
        VT[-1, :] *= -1
        R = np.dot(VT.T, U.T)

    Q_aligned = np.dot(Q_center, R.T)
    return float(np.sqrt(np.mean(np.sum((P_center - Q_aligned) ** 2, axis=1))))


def check_inputs():
    missing = [p for p in (PDB_3ZSJ, PDB_7RGX, PDBQT_DOCKED) if not p.exists()]
    if missing:
        print("Missing required input files:")
        for p in missing:
            print(f"  {p.relative_to(PROJ)}")
        print("\nSee the module docstring for how to obtain them.")
        return False
    return True


def get_ca_atoms(structure):
    ca = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " " and "CA" in residue:
                    ca[residue.id[1]] = residue["CA"]
    return ca


def main():
    if not check_inputs():
        return 1

    print("=" * 72)
    print("57I PYRANOSE-RING DOCKING VALIDATION")
    print("=" * 72)

    # 1. Superpose 7RGX onto 3ZSJ so both ligands share a frame.
    parser = PDBParser(QUIET=True)
    s3zsj = parser.get_structure("3ZSJ", str(PDB_3ZSJ))
    s7rgx = parser.get_structure("7RGX", str(PDB_7RGX))

    ca3 = get_ca_atoms(s3zsj)
    ca7 = get_ca_atoms(s7rgx)
    common = sorted(set(ca3) & set(ca7))

    sup = Superimposer()
    sup.set_atoms([ca3[r] for r in common], [ca7[r] for r in common])
    sup.apply(list(s7rgx.get_atoms()))
    print(f"\n1. Superposed 7RGX onto 3ZSJ over {len(common)} Ca atoms")
    print(f"   Ca RMSD: {sup.rms:.3f} A")

    # 2. Pyranose ring of the aligned crystal ligand.
    crystal = {}
    for model in s7rgx:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == "57I":
                    for atom in residue:
                        name = atom.get_name().strip()
                        if name in RING_ATOMS:
                            crystal[name] = atom.get_vector().get_array()

    if len(crystal) != len(RING_ATOMS):
        print(f"\nERROR: expected {len(RING_ATOMS)} ring atoms in the crystal "
              f"ligand, found {sorted(crystal)}")
        return 1
    print(f"\n2. Crystal pyranose atoms: {RING_ATOMS}")

    # 3. Pyranose ring of the docked pose. The PDBQT ROOT block is the rigid
    #    core written first by the ligand preparation step, so its leading six
    #    atoms correspond to RING_ATOMS in the same order.
    docked = []
    in_model1 = False
    for line in PDBQT_DOCKED.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL 1"):
            in_model1 = True
            continue
        if line.startswith("MODEL 2") or line.startswith("ENDMDL"):
            break
        if not in_model1 or not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            docked.append(np.array([float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54])]))
        except (ValueError, IndexError):
            pass

    if len(docked) < len(RING_ATOMS):
        print(f"\nERROR: docked pose has only {len(docked)} atoms in MODEL 1")
        return 1
    docked_ring = np.array(docked[:len(RING_ATOMS)])
    print(f"3. Docked pyranose atoms: first {len(RING_ATOMS)} of ROOT block")

    # 4. RMSD, before and after optimal superposition.
    crystal_ring = np.array([crystal[name] for name in RING_ATOMS])
    rmsd_direct = float(np.sqrt(np.mean(np.sum((crystal_ring - docked_ring) ** 2, axis=1))))
    rmsd_kabsch = kabsch_rmsd(crystal_ring, docked_ring)

    print("\n" + "=" * 72)
    print("RESULT")
    print("=" * 72)
    print(f"  Pyranose RMSD, direct           {rmsd_direct:.2f} A")
    print(f"  Pyranose RMSD, Kabsch-aligned   {rmsd_kabsch:.2f} A")
    print(f"  Manuscript value (sec. 3.4)     {MANUSCRIPT_RMSD:.2f} A")

    delta = abs(rmsd_kabsch - MANUSCRIPT_RMSD)
    if delta <= TOLERANCE:
        print(f"\n  PASS: reproduces the manuscript value (delta {delta:.3f} A)")
        return 0

    print(f"\n  FAIL: differs from the manuscript value by {delta:.3f} A")
    return 1


if __name__ == "__main__":
    sys.exit(main())
