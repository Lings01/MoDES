# Fresh Install Review (v0.5.0-beta)

## Environment
- Python: 3.10+ required (X | None syntax from ruff auto-fix)
- OS: Linux (CI: ubuntu-latest)
- CI: GitHub Actions, Python 3.10 and 3.11

## Steps
1. `git clone https://github.com/Lings01/MoDES.git`
2. `cd MoDES`
3. `pip install -e .`
4. `pip install -r requirements-dev.txt`
5. `python -m pytest -q` → 91 tests pass
6. `modes --help` → CLI works
7. `python examples/minimal_bulk/run_minimal.py` → runs successfully

## Known Issues
- Python 3.7/3.8/3.9 not supported (requires >=3.10)
- pyproject.toml breaks old pytest <7.0 (use Python >=3.10)
- Local tests: remove pyproject.toml if using Python <3.10

## CI Verification
- pytest: Python 3.10 ✅, Python 3.11 ✅
- lint: ruff check ✅
- build: wheel + twine + fresh install ✅
