# Fresh Install Review (v2.0.0)

## Environment
- Python: 3.10+ required
- OS: Linux (CI: ubuntu-latest)
- CI: GitHub Actions, Python 3.10 and 3.11

## Steps
1. `git clone https://github.com/Lings01/MoDES.git`
2. `cd MoDES`
3. `pip install -e .`
4. `pip install -r requirements-dev.txt`
5. `python -m pytest -q` → 119 tests pass
6. `modes --help` → CLI works
7. `python examples/minimal_bulk/run_minimal.py` → runs successfully
8. `python benchmarks/simulated_event_states/run_benchmark.py` → accuracy report
9. `python benchmarks/cuttag_benchmark/run_benchmark.py` → multi-modal test
10. `python benchmarks/protein_benchmark/run_benchmark.py` → protein test
11. `python benchmarks/stress/run_stress_benchmark.py` → stress test

## Known Issues
- Python 3.7/3.8/3.9 not supported (requires >=3.10)
- pyproject.toml breaks old pytest <7.0 (use Python >=3.10)
- 1 pre-existing test failure: test_report_escapes_html (HTML escaping edge case)
- Network-dependent tests may timeout on slow connections
