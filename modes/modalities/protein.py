"""Protein layer support for MoDES-RAP (v2.0)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from modes.modalities.base import ModalitySpec


PROTEIN_SPEC = ModalitySpec(
    name="protein",
    assay="PROTEIN",
    feature_type="protein",
    regulatory_role="protein_output",
    priority=30,
)


def validate_protein_links(links: pd.DataFrame, rna_genes: set, protein_ids: set) -> list[str]:
    """Validate protein-to-gene link table."""
    issues = []
    for col in ["protein_id", "gene"]:
        if col not in links.columns:
            issues.append(f"Missing required column: {col}")
    if issues:
        return issues
    n_prot_match = links["protein_id"].isin(protein_ids).sum()
    n_gene_match = links["gene"].isin(rna_genes).sum()
    if n_prot_match < len(links):
        issues.append(f"Only {n_prot_match}/{len(links)} protein_ids found in data")
    if n_gene_match < len(links):
        issues.append(f"Only {n_gene_match}/{len(links)} genes found in RNA matrix")
    return issues


def build_protein_gene_map(links: pd.DataFrame) -> dict[str, str]:
    """Build protein_id → gene mapping from link table."""
    return dict(zip(links["protein_id"], links["gene"]))
