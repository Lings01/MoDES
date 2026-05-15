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
) -> float:
    """
    Compute quality score q_e in [0, 1] for a feature.

    Factors:
      - mean non-zero proportion across samples
      - mean count level (log-scale, normalized)
      - batch association (Cramer's V, if batch provided)

    Returns a score in [0, 1] where 1 is highest quality.
    """
    counts = np.asarray(counts, dtype=float)
    n = len(counts)

    if n == 0:
        return 0.0

    # Fraction of samples with non-zero counts
    nonzero_frac = np.mean(counts > 0)

    # Mean expression level (log-scale, normalized to [0,1] via logistic)
    mean_count = np.mean(counts)
    if mean_count > 0:
        expr_score = 1.0 - np.exp(-mean_count / max(np.mean(counts[counts > 0]) if np.any(counts > 0) else 1, 1))
    else:
        expr_score = 0.0

    # Batch association score (1 = no batch effect)
    batch_score = 1.0
    if batch_labels is not None and len(np.unique(batch_labels)) > 1:
        try:
            from scipy.stats import kruskal
            groups = [counts[batch_labels == b] for b in np.unique(batch_labels)]
            _, pval = kruskal(*groups)
            batch_score = pval  # high p-value = no batch effect = high quality
        except Exception:
            pass

    # Weighted geometric mean
    scores = [nonzero_frac, expr_score, batch_score]
    weights = [0.4, 0.3, 0.3]
    log_score = sum(w * np.log(max(s, 1e-10)) for w, s in zip(weights, scores))
    quality = np.exp(log_score)

    return float(np.clip(quality, 0.0, 1.0))


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
        qualities[i] = compute_quality_score(
            count_matrix[:, i], metadata, batch_labels
        )
    return qualities
