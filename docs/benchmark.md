# Benchmark

MoDES includes a synthetic benchmark suite for evaluating state recovery accuracy.

## Simulated Event States

```bash
python benchmarks/simulated_event_states/run_benchmark.py
```

Generates 50 events (10 per state: concordant, chromatin_primed, rna_only,
discordant_opposite, null) across 10 chromosomes with known ground truth,
runs MoDES, and computes per-state metrics.

### Output

| File | Description |
|---|---|
| `truth.tsv` | Ground truth event states |
| `moDES_output/event_table.tsv` | MoDES predicted states |
| `metrics.tsv` | Per-state precision, recall, F1 |

### Metrics

| Metric | Formula | Range |
|---|---|---|
| Accuracy | (correct) / (total) | [0, 1] |
| Precision | TP / (TP + FP) | [0, 1] |
| Recall | TP / (TP + FN) | [0, 1] |
| F1 | 2 × P × R / (P + R) | [0, 1] |

## Planned Benchmarks

- **Semi-real spike-in**: Real count matrices with controlled effect injection
- **Negative control**: Shuffled condition labels (expect all null)
- **Baseline comparison**: Naive overlap vs MoDES state classification
- **Runtime profiling**: By data scale (samples × events)
