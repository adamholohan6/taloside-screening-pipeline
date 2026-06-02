from pathlib import Path

p = Path("src/taloside_pipeline/phase3_docking_FIXED.py")
if not p.exists():
    raise SystemExit(f"File not found: {p}")

text = p.read_text(encoding="utf-8")

# Replace literal backslash-n sequences with real newlines
fixed = text.replace("\\n", "\n")

# Also ensure we didn't accidentally remove any intended escape sequences
p.write_text(fixed, encoding="utf-8")
print("Fixed file written to", p)
