"""Multi-condition contrast matrix + pseudotime lag inference (v2.0)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass


@dataclass
class ContrastResult:
    """Per-contrast event result."""
    contrast_id: str
    ref: str
    target: str
    coef: float
    se: float
    p_value: float
    fdr: float
    direction: int


def build_contrast_matrix(
    condition_values: np.ndarray,
    contrasts: list[tuple[str, str]],
) -> tuple[np.ndarray, list[str]]:
    """
    Build a design matrix for multi-condition comparisons.

    Parameters
    ----------
    condition_values : ndarray
        Categorical condition labels for each sample.
    contrasts : list of (reference, target) tuples.

    Returns
    -------
    X_contrast : ndarray, shape (n_samples, n_contrasts + 1) with intercept
    col_names : list of str
    """
    n = len(condition_values)
    categories = sorted(set(condition_values))
    cat_to_idx = {c: i for i, c in enumerate(categories)}

    contrast_vecs = []
    col_names = ["intercept"]
    for ref, tgt in contrasts:
        if ref not in cat_to_idx or tgt not in cat_to_idx:
            raise ValueError(f"Contrast ({ref}, {tgt}) uses unknown categories")
        vec = np.zeros(n)
        vec[condition_values == tgt] = 1.0
        vec[condition_values == ref] = 0.0
        # Sparse: only include samples from ref or tgt
        mask = (condition_values == ref) | (condition_values == tgt)
        contrast_vecs.append((vec, mask, f"{tgt}_vs_{ref}"))
        col_names.append(f"{tgt}_vs_{ref}")

    X = np.column_stack([np.ones(n)] + [cv[0] for cv in contrast_vecs])
    return X, col_names


def estimate_pseudotime_lag(
    atac_values: np.ndarray,
    rna_values: np.ndarray,
    pseudotime: np.ndarray,
    max_lag: int = 10,
) -> dict:
    """
    Estimate ATAC→RNA time delay via cross-correlation.

    Parameters
    ----------
    atac_values, rna_values : ndarray
        Per-cell/pseudotime-point values.
    pseudotime : ndarray
        Pseudotime values sorted ascending.
    max_lag : int
        Maximum lag to test.

    Returns
    -------
    dict with: best_lag, correlation, lag_profile
    """
    order = np.argsort(pseudotime)
    atac_sorted = atac_values[order]
    rna_sorted = rna_values[order]

    correlations = []
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            corr = np.corrcoef(atac_sorted[-lag:], rna_sorted[:lag])[0, 1]
        elif lag > 0:
            corr = np.corrcoef(atac_sorted[:-lag], rna_sorted[lag:])[0, 1]
        else:
            corr = np.corrcoef(atac_sorted, rna_sorted)[0, 1]
        correlations.append((lag, corr if not np.isnan(corr) else 0.0))

    correlations.sort(key=lambda x: abs(x[1]), reverse=True)
    best_lag = correlations[0][0]
    return {
        "best_lag": int(best_lag),
        "lag_correlation": float(correlations[0][1]),
        "lag_profile": {int(l): float(c) for l, c in correlations},
    }
