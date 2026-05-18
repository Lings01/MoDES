"""Statistical utilities for MoDES."""

import numpy as np


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.

    Returns q-values (FDR-adjusted p-values).
    """
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    # Sort p-values and track original order
    order = np.argsort(p)
    sorted_p = p[order]
    # Compute BH critical values
    ranks = np.arange(1, n + 1)
    bh_values = sorted_p * n / ranks
    # Ensure monotonicity
    q_sorted = np.minimum.accumulate(bh_values[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    # Restore original order
    qvalues = np.empty(n)
    qvalues[order] = q_sorted
    return qvalues


def empirical_bayes_moderate(
    coefs: np.ndarray, ses: np.ndarray, df: int = 3
) -> np.ndarray:
    """
    Limma-style empirical Bayes variance moderation.

    Pools variance estimates across all features to produce
    more stable moderated standard errors.

    Parameters
    ----------
    coefs : ndarray
        Coefficient estimates.
    ses : ndarray
        Standard errors of coefficients.
    df : int
        Prior degrees of freedom (default 3).

    Returns
    -------
    moderated_ses : ndarray
        Moderated standard errors.
    """
    coefs = np.asarray(coefs, dtype=float)
    ses = np.asarray(ses, dtype=float)

    # Remove NaN/Inf entries for prior estimation
    valid = np.isfinite(coefs) & np.isfinite(ses) & (ses > 0)
    if valid.sum() < 5:
        return ses.copy()

    s2 = ses[valid] ** 2
    # Method of moments to estimate prior parameters
    # Assume s2 ~ InverseGamma(df0/2, df0*s02/2)
    # log(s2) is approximately normal with mean = log(s02) + psi(df0/2) - log(df0/2)
    # and variance = trigamma(df0/2)
    log_s2 = np.log(s2)

    # Estimate s02 as exp(mean(log_s2)) -- geometric mean
    s02 = np.exp(np.mean(log_s2))

    # Estimate df0 by matching variance of log(s2) to trigamma(df0/2)
    var_log_s2 = np.var(log_s2)
    trigamma_val = max(var_log_s2, 1e-10)

    df0 = _trigamma_inverse(trigamma_val) * 2

    # Clamp df0 to reasonable range
    df0 = max(0.5, min(df0, 50.0))

    # Compute posterior variances
    n = df  # residual df from GLM
    s2_post = (df0 * s02 + n * s2) / (df0 + n)

    moderated_ses = ses.copy()
    moderated_ses[valid] = np.sqrt(s2_post)
    return moderated_ses


def _trigamma_inverse(x: float) -> float:
    """Approximate inverse of the trigamma function."""
    if x >= 10:
        return 1.0 / (x - 0.5 / (x - 0.25 / x))
    elif x >= 1:
        return 1.0 / x + 0.5
    elif x > 0.5:
        return 1.0 / x + 0.3
    else:
        return 1.0 / max(x, 1e-10)


def compute_quality_score(
    counts: np.ndarray,
    batch_labels: np.ndarray = None,
) -> dict:
    """
    Compute quality components for a feature. Returns dict with:
      - detection_score: fraction of non-zero samples [0, 1]
      - depth_score: log-scaled mean count [0, 1]
      - batch_score: Kruskal-Wallis p-value [0, 1] (1 = no batch effect)
      - quality_score: weighted geometric mean [0, 1]
    """
    counts = np.asarray(counts, dtype=float)
    n = len(counts)
    if n == 0:
        return {"detection_score": 0.0, "depth_score": 0.0,
                "batch_score": 1.0, "quality_score": 0.0}

    nonzero_frac = np.mean(counts > 0)
    mean_count = np.mean(counts)
    # Reference log-scale: log(1 + ref_count) as anchor so depth_score reflects
    # meaningful dynamic range instead of always saturating at ~1.
    ref_count = max(np.percentile(counts, 95), 100.0)
    depth_score = min(np.log(mean_count + 1) / max(np.log(ref_count + 1), 0.01), 1.0) if mean_count > 0 else 0.0

    batch_score = 1.0
    if batch_labels is not None and len(np.unique(batch_labels)) > 1:
        try:
            from scipy.stats import kruskal
            groups = [counts[batch_labels == b] for b in np.unique(batch_labels)]
            _, pval = kruskal(*groups)
            batch_score = float(np.clip(pval, 0.0, 1.0))
        except Exception:
            pass

    scores = [nonzero_frac, depth_score, batch_score]
    weights = [0.4, 0.3, 0.3]
    log_score = sum(w * np.log(max(s, 1e-10)) for w, s in zip(weights, scores))
    quality = np.exp(log_score)

    return {
        "detection_score": float(nonzero_frac),
        "depth_score": float(depth_score),
        "batch_score": float(batch_score),
        "quality_score": float(np.clip(quality, 0.0, 1.0)),
    }


def compute_feature_quality_scores(
    count_matrix: np.ndarray,
    metadata: np.ndarray = None,
    batch_labels: np.ndarray = None,
) -> np.ndarray:
    """
    Compute quality scores for all features in a count matrix.

    Parameters
    ----------
    count_matrix : ndarray of shape (n_samples, n_features)
    metadata : ndarray, optional
    batch_labels : ndarray of shape (n_samples,), optional

    Returns
    -------
    qualities : ndarray of shape (n_features,)
    """
    n_features = count_matrix.shape[1]
    qualities = np.zeros(n_features)
    for i in range(n_features):
        qc = compute_quality_score(
            count_matrix[:, i],
            batch_labels=batch_labels,
        )
        qualities[i] = qc["quality_score"]
    return qualities


def parse_genomic_interval(feature_id: str) -> tuple[str, int, int] | None:
    """Parse a genomic interval string into (chr, start, end).

    Supports formats:
      - chr1:100-200
      - chr1:100-200|H3K27ac
      - chr1_100_200

    Returns None if not parseable.
    """
    import re
    # Strip suffix after | (e.g., CUT&Tag target annotation)
    clean = str(feature_id).split("|")[0]
    m = re.match(r"(chr[\w]+)[:\-_](\d+)[:\-_](\d+)", clean)
    if m:
        return (m.group(1), int(m.group(2)), int(m.group(3)))
    return None


def interval_overlap(
    region_a: str, region_b: str, min_reciprocal: float = 0.5
) -> dict | None:
    """Compute overlap between two genomic interval feature IDs.

    Returns dict with overlap_bp, reciprocal_overlap_a, reciprocal_overlap_b,
    min_reciprocal_overlap, and match (bool) if both are parseable intervals.
    Returns None if either is not a genomic interval.
    """
    a = parse_genomic_interval(region_a)
    b = parse_genomic_interval(region_b)
    if a is None or b is None:
        return None

    chrom_a, start_a, end_a = a
    chrom_b, start_b, end_b = b

    if chrom_a != chrom_b:
        return {
            "overlap_bp": 0, "reciprocal_overlap_a": 0.0,
            "reciprocal_overlap_b": 0.0, "min_reciprocal_overlap": 0.0,
            "match": False, "method": "genomic_interval",
        }

    width_a = end_a - start_a
    width_b = end_b - start_b
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    overlap_bp = max(0, overlap_end - overlap_start)

    recip_a = overlap_bp / max(width_a, 1)
    recip_b = overlap_bp / max(width_b, 1)
    min_recip = min(recip_a, recip_b)

    return {
        "overlap_bp": overlap_bp,
        "reciprocal_overlap_a": round(recip_a, 4),
        "reciprocal_overlap_b": round(recip_b, 4),
        "min_reciprocal_overlap": round(min_recip, 4),
        "match": min_recip >= min_reciprocal,
        "method": "genomic_interval",
    }
