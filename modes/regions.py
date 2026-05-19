"""Genomic region matching for multi-modal feature alignment.

Provides interval-overlap-based matching between query regions (e.g., ATAC peaks)
and target regions (e.g., CUT&Tag or ChIP peaks). Used by EvidenceBuilder and
_build_modality_evidence for cross-assay feature matching.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from modes.utils import parse_genomic_interval


def match_regions_by_overlap(
    query_regions: list[str],
    target_regions: list[str],
    min_reciprocal_overlap: float = 0.5,
    min_overlap_bp: int = 50,
) -> pd.DataFrame:
    """Match query regions to target regions by reciprocal interval overlap.

    Parameters
    ----------
    query_regions : list of str
        Query region IDs (e.g., ATAC peaks in "chr1:100-200" format).
    target_regions : list of str
        Target region IDs (e.g., CUT&Tag peaks, may have |suffix).
    min_reciprocal_overlap : float
        Minimum min(overlap/width_query, overlap/width_target) for a match.
    min_overlap_bp : int
        Minimum absolute overlap in base pairs.

    Returns
    -------
    DataFrame with columns:
        query_region, target_region, overlap_bp, reciprocal_overlap_query,
        reciprocal_overlap_target, min_reciprocal_overlap, region_match_score,
        match
    """
    # Parse all regions
    query_parsed = {}
    for q in query_regions:
        parsed = parse_genomic_interval(q)
        if parsed:
            query_parsed[q] = parsed

    target_parsed = {}
    for t in target_regions:
        parsed = parse_genomic_interval(t)
        if parsed:
            target_parsed[t] = parsed

    rows = []
    for q_name, (q_chr, q_start, q_end) in query_parsed.items():
        q_width = q_end - q_start
        for t_name, (t_chr, t_start, t_end) in target_parsed.items():
            if q_chr != t_chr:
                continue
            overlap_start = max(q_start, t_start)
            overlap_end = min(q_end, t_end)
            overlap_bp = max(0, overlap_end - overlap_start)
            if overlap_bp < min_overlap_bp:
                continue

            recip_q = overlap_bp / max(q_width, 1)
            recip_t = overlap_bp / max(t_end - t_start, 1)
            min_recip = min(recip_q, recip_t)
            match = min_recip >= min_reciprocal_overlap

            rows.append({
                "query_region": q_name,
                "target_region": t_name,
                "overlap_bp": overlap_bp,
                "reciprocal_overlap_query": round(recip_q, 4),
                "reciprocal_overlap_target": round(recip_t, 4),
                "min_reciprocal_overlap": round(min_recip, 4),
                "region_match_score": round(min_recip, 4) if match else 0.0,
                "match": match,
            })

    return pd.DataFrame(rows)


def find_best_region_match(
    query: str,
    target_dict: dict,
    min_reciprocal_overlap: float = 0.5,
    min_overlap_bp: int = 50,
) -> tuple | None:
    """Find the best-matching target region for a single query.

    Returns (target_key, target_value, overlap_info_dict) or None.
    Target_dict maps region_id -> value (e.g., ModalityEffect).
    """
    best_key = None
    best_value = None
    best_info = None
    best_score = 0.0

    q_parsed = parse_genomic_interval(query)
    if q_parsed is None:
        return None

    q_chr, q_start, q_end = q_parsed
    q_width = q_end - q_start

    for t_name in target_dict:
        t_parsed = parse_genomic_interval(t_name)
        if t_parsed is None:
            continue
        t_chr, t_start, t_end = t_parsed
        if t_chr != q_chr:
            continue

        overlap_start = max(q_start, t_start)
        overlap_end = min(q_end, t_end)
        overlap_bp = max(0, overlap_end - overlap_start)
        if overlap_bp < min_overlap_bp:
            continue

        recip_q = overlap_bp / max(q_width, 1)
        recip_t = overlap_bp / max(t_end - t_start, 1)
        min_recip = min(recip_q, recip_t)

        if min_recip >= min_reciprocal_overlap and min_recip > best_score:
            best_key = t_name
            best_value = target_dict[t_name]
            best_info = {
                "overlap_bp": overlap_bp,
                "reciprocal_overlap_query": round(recip_q, 4),
                "reciprocal_overlap_target": round(recip_t, 4),
                "min_reciprocal_overlap": round(min_recip, 4),
                "region_match_score": round(min_recip, 4),
            }
            best_score = min_recip

    if best_key is not None:
        return (best_key, best_value, best_info)
    return None
