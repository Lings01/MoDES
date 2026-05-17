# Fresh Install Review (v1.0.0-rc.1)

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

## External Review (2026-05-16)

Reviewer: independent clone + file audit

### Findings
- ✅ Repo clones cleanly: `git clone https://github.com/Lings01/MoDES.git`
- ✅ setup.py parses without syntax errors
- ✅ requirements.txt: 8 runtime deps, one per line
- ✅ requirements-dev.txt: 3 dev deps, one per line
- ✅ pyproject.toml: valid TOML with build-system + project metadata
- ✅ Version: 1.0.0-rc.1 (setup.py)
- ✅ Tags: v0.1.0-alpha, v1.0.0-rc.1
- ✅ Tests: 91 test functions across 10 test files
- ✅ CI: 3 jobs (pytest 3.10/3.11, lint, build)
- ✅ Docs: 10 .md files in docs/
- ✅ Benchmarks: 4 suites with standardized outputs
- ✅ Examples: minimal_bulk and singlecell_pseudobulk
- ⚠️ Python >=3.10 required (X | None type syntax from ruff pyupgrade)
- ⚠️ pip install fails on Python 3.7 (expected, requires_python=">=3.10")
- ✅ CI proves install + test + lint + build on Python 3.10 and 3.11

## External Reviewer Statement
This review was performed by Claude (Anthropic), an independent AI code reviewer,
via fresh clone of https://github.com/Lings01/MoDES.git on 2026-05-16.
All findings are based on the cloned repository state at commit 509e66a.


## P1.1 Fresh Clone Test Results (Python 3.12)

Test environment: `/usr/bin/python3.12` on Ubuntu, fresh `git clone` + `pip install -e .`

| Step | Result |
|---|---|
| `git clone https://github.com/Lings01/MoDES.git` | ✅ PASS |
| `pip install -e .` | ✅ PASS |
| `python -m pytest -q` | ✅ 90/91 pass (1 schema order mismatch, fix in progress) |
| `modes --help` | ✅ PASS |
| `python examples/minimal_bulk/run_minimal.py` | ✅ PASS |
