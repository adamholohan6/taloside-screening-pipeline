# Taloside Screening Pipeline – Phase 2

**Virtual library generation, drug‑likeness filtering, and PAINS screening for taloside derivatives.**

## 📌 Overview
- Generates virtual taloside derivatives via click chemistry (16 compounds)
- Calculates MW, LogP, TPSA, H‑bond donors/acceptors, rotatable bonds
- Applies Lipinski’s Rule of 5 (carbohydrate‑adjusted thresholds)
- Screens for PAINS (Pan‑Assay Interference Compounds)

## ✨ Results
- **16 compounds** generated
- **100% pass Lipinski** (MW ≤600, LogP ≤4, HBD ≤6, HBA ≤12)
- **No PAINS flagged** – suitable for experimental screening

## 🚀 Quick Start
```bash
cd src/taloside_pipeline
python phase2_integration.py