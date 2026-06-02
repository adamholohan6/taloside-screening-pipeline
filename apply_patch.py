import io
from pathlib import Path

p = Path("src/taloside_pipeline/phase3_docking_FIXED.py")
if not p.exists():
    raise SystemExit(f"File not found: {p}")

text = p.read_text(encoding="utf-8")

start_sig = "def mol_to_pdbqt_string(mol: Chem.Mol, conf_id: int = 0, name: str = \"LIG\") -> str:"
if start_sig not in text:
    raise SystemExit("Could not find mol_to_pdbqt_string signature in file.")

start_idx = text.index(start_sig)
next_def = text.find("\ndef embed_ligand_3d", start_idx)
if next_def == -1:
    raise SystemExit("Could not locate end of mol_to_pdbqt_string function.")

new_func = r'''def mol_to_pdbqt_string(mol: Chem.Mol, conf_id: int = 0, name: str = "LIG") -> str:
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

    lines: List[str] = [f"REMARK  Name = {name}\\n"]
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
        atom_line = atom_line.ljust(76) + f"{ad_type:>2s}\\n"
        lines.append(atom_line)

    lines.append("TORSDOF 0\\n")
    return "".join(lines)
'''

new_text = text[:start_idx] + new_func + text[next_def:]

# Replace Meeko check in validate_receptor (tolerant replacement)
new_text = new_text.replace("mol_to_pdbqt_string_meeko(embedded, conf_id=0, name=\"lactose\")",
                            "mol_to_pdbqt_string(embedded, conf_id=0, name=\"lactose\")")

p.write_text(new_text, encoding="utf-8")
print("Patch applied to", p)
