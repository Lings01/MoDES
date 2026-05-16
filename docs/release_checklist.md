# Release Checklist

## Pre-release
- [x] CI green (pytest + lint) — 3 jobs on Python 3.10/3.11
- [x] All unit tests pass (`python -m pytest -q`) — 91 tests
- [x] Lint passes (`ruff check modes/ tests/ benchmarks/ examples/`)
- [x] Quick benchmarks pass (`--quick` mode, synthetic + negative control)
- [x] Example runs (`python examples/minimal_bulk/run_minimal.py`)
- [x] README up to date with current capabilities
- [x] CHANGELOG.md updated for this version
- [x] ROADMAP.md reflects current plans
- [x] docs/output_schema.md matches actual output columns
- [x] docs/api_reference.md complete

## Release
- [x] Version bumped in setup.py (`0.5.0-beta`)
- [x] pyproject.toml present for build
- [x] Git tag created (`v0.5.0-beta`)
- [x] Tag pushed (`git push origin --tags`)
- [ ] GitHub Release page created with changelog notes

## Post-release
- [x] CI badge shows green in README
- [ ] PyPI publish (when ready for v1.0)
