# MoDES Benchmarks

## Quick Start

```bash
# Synthetic event states (100% accuracy on controlled data)
python benchmarks/simulated_event_states/run_benchmark.py --quick

# Semi-real spike-in (realistic library sizes)
python benchmarks/semi_real_spikein/run_spikein.py

# Negative control (shuffled labels → expect all null)
python benchmarks/negative_control/run_negative_control.py

# Baseline comparison (MoDES vs naive overlap)
python benchmarks/baseline_comparison/run_baseline.py
```

## How to Read Results

Each benchmark produces output files in its `output/` directory:

| File | Description |
|---|---|
| `truth.tsv` | Ground truth event states |
| `metrics.tsv` | Per-state precision, recall, F1 |
| `confusion_matrix.tsv` | State × State confusion matrix |
| `runtime.tsv` | Wall-clock execution time |
| `moDES_output/event_table.tsv` | Full MoDES output |

**Expected behavior:**
- Synthetic: accuracy near 1.0 on controlled data
- Semi-real spike-in: accuracy ≥ 0.85 with realistic noise
- Negative control: null fraction ≥ 0.9 (shuffled labels)
- Baseline comparison: MoDES matches or exceeds naive overlap on controlled data
