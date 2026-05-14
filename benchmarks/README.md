# Benchmarks

## simulated_event_states

Synthetic benchmark evaluating MoDES state recovery against known ground truth.

```bash
python benchmarks/simulated_event_states/run_benchmark.py
```

Output:
- `truth.tsv` — ground truth event states
- `moDES_output/` — MoDES output files
- `metrics.tsv` — per-state precision, recall, F1

### Metrics

| Metric | Description |
|---|---|
| Accuracy | Fraction of events correctly classified |
| Per-state precision | TP / (TP + FP) per state |
| Per-state recall | TP / (TP + FN) per state |
| Per-state F1 | Harmonic mean of precision and recall |

### Planned

- Semi-real spike-in benchmark (real count matrices with controlled effects)
- Negative control (shuffled labels)
- Baseline comparison (naive overlap vs MoDES)
- Runtime / memory profiling by data scale
