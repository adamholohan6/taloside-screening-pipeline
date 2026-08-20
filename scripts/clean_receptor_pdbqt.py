from __future__ import annotations

from pathlib import Path


FORBIDDEN_RESNAMES = {"BGC", "GAL", "HOH"}


def _resname(line: str) -> str:
    return line[17:20].strip() if len(line) >= 20 else ""


def clean_pdbqt(source: Path, target: Path) -> None:
    lines = []
    for line in source.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and _resname(line) in FORBIDDEN_RESNAMES:
            continue
        lines.append(line)
    target.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    source = root / "data" / "docking" / "3ZSJ.pdbqt"
    target = root / "data" / "docking" / "3ZSJ_clean.pdbqt"
    clean_pdbqt(source, target)
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
