# Release Checklist

## Pre-release

- [ ] CI green (pytest + lint)
- [ ] All unit tests pass (`python -m pytest -q`)
- [ ] Lint passes (`ruff check modes/ tests/ benchmarks/ examples/`)
- [ ] Quick benchmarks pass
- [ ] Example runs (`python examples/minimal_bulk/run_minimal.py`)
- [ ] README up to date with current capabilities
- [ ] CHANGELOG.md updated for this version
- [ ] ROADMAP.md reflects current plans
- [ ] `docs/release_checklist.md` reviewed
- [ ] `docs/output_schema.md` matches actual output columns

## Release

- [ ] Version bumped in `setup.py`
- [ ] Version bumped in `CHANGELOG.md`
- [ ] Git tag created (`git tag vX.Y.Z`)
- [ ] Tag pushed (`git push origin --tags`)

## Post-release

- [ ] GitHub Release page created with changelog notes
- [ ] CI badge shows green in README
