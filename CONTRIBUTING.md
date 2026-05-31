# Contributing to Taloside Screening Pipeline

Thank you for your interest in contributing! This document provides guidelines.

## How to Contribute

### Reporting Issues
1. Check existing issues to avoid duplicates
2. Provide a clear description of the problem
3. Include error messages or logs if applicable
4. Specify Python and RDKit versions

### Suggesting Enhancements
1. Describe the use case
2. Explain expected behavior
3. Include examples if relevant

### Adding New Compounds

To add new taloside variants:

1. Generate SMILES string (e.g., from ChemDraw)
2. Edit `generate_library.py` → `load_compounds_from_dict()`
3. Add entry: `"Talo-XX": "YOUR_SMILES_STRING"`
4. Run tests: `python generate_library.py`
5. Verify output CSV

### Submitting Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Follow PEP 8 style guidelines
4. Add docstrings to all functions
5. Test thoroughly
6. Commit clearly: `git commit -m "Add feature description"`
7. Push: `git push origin feature/your-feature`
8. Open a pull request

## Code Standards

- **Python Version:** 3.6+
- **Style:** PEP 8
- **Linting:** All code should pass `flake8` checks
- **Documentation:** Docstrings for all functions
- **Testing:** Include usage examples

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR-USERNAME/taloside-screening-pipeline.git
cd taloside-screening-pipeline

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Test the script
python generate_library.py
```

## Questions?

Feel free to open a GitHub discussion or issue!
