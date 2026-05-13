"""Shared test fixtures with synthetic ground truth data."""

import pytest
import numpy as np
import pandas as pd


def _generate_counts(means, size=100):
    """Generate NB counts given means."""
    if means <= 0:
        means = 0.1
    dispersion = 0.5
    n = means / dispersion
    p = 1 / (1 + dispersion)
    return np.random.negative_binomial(n, p, size=size)


def _build_tss_map(n_genes=30):
    """Build synthetic TSS map with positions matching peak coordinates."""
    tss_map = {}
    for i in range(n_genes):
        gene_name = f"gene_{i}"
        tss_pos = 1000 + i * 500 + 250  # center of peak
        tss_map[gene_name] = (gene_name, "chr1", tss_pos)
    return tss_map


def _build_synthetic_data(
    n_control=10,
    n_treatment=10,
    n_genes=30,
    n_peaks=60,
    seed=42,
) -> tuple:
    """
    Build synthetic bulk RNA+ATAC data with known ground truth events.

    Ground truth:
      - 5 concordant (ATAC up, RNA up in treatment)
      - 5 chromatin_primed (ATAC up, RNA null)
      - 5 rna_only (ATAC null, RNA up)
      - 5 null (both null)
    """
    rng = np.random.default_rng(seed)
    n_total = n_control + n_treatment

    condition = np.array(["control"] * n_control + ["treatment"] * n_treatment)
    # Ensure donors appear in both conditions to avoid confounding
    n_unique_donors = max(2, n_total // 4)
    donor = np.array([f"donor_{i % n_unique_donors}" for i in range(n_total)])
    # Batch: groups of 2 to avoid confounding with donor
    batch = np.array([
        "batch_A" if (i // 2) % 2 == 0 else "batch_B" for i in range(n_total)
    ])

    base_expression = rng.lognormal(mean=3.0, sigma=0.5, size=n_genes)
    base_accessibility = rng.lognormal(mean=2.0, sigma=0.5, size=n_peaks)

    gene_names = [f"gene_{i}" for i in range(n_genes)]
    peak_names = [f"chr1:{1000 + i * 500}-{1500 + i * 500}" for i in range(n_peaks)]

    rna = np.zeros((n_total, n_genes))
    atac = np.zeros((n_total, n_peaks))

    gt_events = []
    gene_idx = 0
    peak_idx = 0

    # 5 concordant events
    for i in range(5):
        gt_events.append({
            "event_id": f"concordant_{i}",
            "gene": gene_names[gene_idx],
            "peak_id": peak_names[peak_idx],
            "true_state": "concordant",
            "atac_effect": 1.5,
            "rna_effect": 1.2,
            "rna_after_atac": 0.2,
        })
        gene_idx += 1
        peak_idx += 1

    # 5 chromatin_primed events
    for i in range(5):
        gt_events.append({
            "event_id": f"primed_{i}",
            "gene": gene_names[gene_idx],
            "peak_id": peak_names[peak_idx],
            "true_state": "chromatin_primed",
            "atac_effect": 1.5,
            "rna_effect": 0.0,
            "rna_after_atac": 0.0,
        })
        gene_idx += 1
        peak_idx += 1

    # 5 rna_only events
    for i in range(5):
        gt_events.append({
            "event_id": f"rna_only_{i}",
            "gene": gene_names[gene_idx],
            "peak_id": peak_names[peak_idx],
            "true_state": "rna_only",
            "atac_effect": 0.0,
            "rna_effect": 1.2,
            "rna_after_atac": 1.0,
        })
        gene_idx += 1
        peak_idx += 1

    # 5 null events
    for i in range(5):
        gt_events.append({
            "event_id": f"null_{i}",
            "gene": gene_names[gene_idx],
            "peak_id": peak_names[peak_idx],
            "true_state": "null",
            "atac_effect": 0.0,
            "rna_effect": 0.0,
            "rna_after_atac": 0.0,
        })
        gene_idx += 1
        peak_idx += 1

    # Generate counts
    for s in range(n_total):
        is_treated = condition[s] == "treatment"
        for i in range(n_genes):
            mu = base_expression[i]
            if is_treated:
                for ev in gt_events:
                    if ev["gene"] == gene_names[i]:
                        mu *= np.exp(ev["rna_effect"])
                        break
            rna[s, i] = _generate_counts(mu, size=1)[0]
        for i in range(n_peaks):
            mu = base_accessibility[i]
            if is_treated:
                for ev in gt_events:
                    if ev["peak_id"] == peak_names[i]:
                        mu *= np.exp(ev["atac_effect"])
                        break
            atac[s, i] = _generate_counts(mu, size=1)[0]

    rna += rng.poisson(0.1, rna.shape)
    atac += rng.poisson(0.1, atac.shape)

    obs = pd.DataFrame({
        "condition": condition,
        "donor": donor,
        "batch": batch,
    }, index=[f"sample_{i}" for i in range(n_total)])

    rna_df = pd.DataFrame(rna, index=obs.index, columns=gene_names)
    atac_df = pd.DataFrame(atac, index=obs.index, columns=peak_names)
    gt_df = pd.DataFrame(gt_events)
    tss_map = _build_tss_map(n_genes)

    return rna_df, atac_df, obs, gt_df, tss_map


@pytest.fixture
def synthetic_tss_map():
    return _build_tss_map(30)


@pytest.fixture
def synthetic_bulk_data():
    """Return MoDEData and ground truth for synthetic bulk experiment."""
    from modes.data import MoDEData
    rna_df, atac_df, obs, gt, tss_map = _build_synthetic_data(seed=42)
    data = MoDEData(rna=rna_df, atac=atac_df, obs=obs)
    return data, gt, tss_map


@pytest.fixture
def synthetic_bulk_data_small():
    """Small version for faster tests."""
    from modes.data import MoDEData
    rna_df, atac_df, obs, gt, tss_map = _build_synthetic_data(
        n_control=5, n_treatment=5, n_genes=20, n_peaks=25, seed=123
    )
    data = MoDEData(rna=rna_df, atac=atac_df, obs=obs)
    return data, gt, tss_map
