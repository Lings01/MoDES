"""State grammar: pluggable rule registry for multi-modal event states.

Each modality contributes state rules. The StateClassifier composes them
into a final biological state + per-modality state components.
"""

from __future__ import annotations

# --- State rule grammar ---

# RNA+ATAC core states (MoDES-RA)
RA_STATES = {
    "concordant": {
        "required": {"atac": "sig", "rna": "sig"},
        "direction": "same",
        "description": "Local chromatin opening drives transcription",
    },
    "chromatin_primed": {
        "required": {"atac": "sig", "rna": "not_sig"},
        "description": "Chromatin open, transcription not started",
    },
    "rna_only": {
        "required": {"rna": "sig", "atac": "not_sig"},
        "description": "RNA change not explained by local chromatin",
    },
    "discordant_opposite": {
        "required": {"atac": "sig", "rna": "sig"},
        "direction": "opposite",
        "description": "Opposite direction across layers",
    },
    "null": {
        "required": {"atac": "not_sig", "rna": "not_sig"},
        "description": "No significant change",
    },
}

# CUT&Tag activating mark states (e.g., H3K27ac, H3K4me3)
EPI_ACTIVATING_STATES = {
    "epigenomic_concordant": {
        "required": {"epi_activating": "sig", "rna": "sig"},
        "direction": "same",
        "description": "Active mark and RNA both up",
    },
    "active_enhancer_primed": {
        "required": {"epi_activating": "sig", "rna": "not_sig"},
        "description": "Active histone mark present, RNA not yet responding",
    },
    "mark_only": {
        "required": {"epi_activating": "sig", "atac": "not_sig", "rna": "not_sig"},
        "description": "Mark change without accessibility or RNA change",
    },
}

# CUT&Tag repressive mark states (e.g., H3K27me3, H3K9me3)
EPI_REPRESSIVE_STATES = {
    "repressive_concordant": {
        "required": {"epi_repressive": "sig", "rna": "sig"},
        "direction": "opposite",
        "description": "Repressive mark up, RNA down",
    },
    "derepression": {
        "required": {"epi_repressive": "sig", "rna": "sig"},
        "direction": "same",
        "description": "Repressive mark down, RNA up",
    },
    "repressive_primed": {
        "required": {"epi_repressive": "sig", "rna": "not_sig"},
        "description": "Repressive mark established, RNA unchanged",
    },
}

# Protein layer states (MoDES-RAP)
PROTEIN_STATES = {
    "full_activation": {
        "required": {"atac": "sig", "rna": "sig", "protein": "sig"},
        "direction": "same",
        "description": "Complete regulatory chain: chromatin→RNA→protein",
    },
    "protein_buffered": {
        "required": {"rna": "sig", "protein": "not_sig"},
        "description": "RNA changes but protein does not",
    },
    "protein_memory": {
        "required": {"rna": "not_sig", "protein": "sig"},
        "description": "Protein persists after RNA returns to baseline",
    },
    "protein_opposite": {
        "required": {"rna": "sig", "protein": "sig"},
        "direction": "opposite",
        "description": "RNA and protein move in opposite directions",
    },
}

# Spatial states (MoDES-Spatial)
SPATIAL_STATES = {
    "spatial_region_specific": {
        "required": {"spatial_region": "sig"},
        "description": "Event signal localized to specific anatomical region",
    },
    "spatial_niche_driven": {
        "required": {"spatial_neighbor": "sig"},
        "description": "Event explained by neighbor cell composition",
    },
    "cell_intrinsic": {
        "required": {"spatial_neighbor": "not_sig"},
        "description": "Event independent of spatial context",
    },
    "spatial_edge_artifact": {
        "required": {"spatial_edge": "sig"},
        "description": "High signal near tissue edge — possible artifact",
    },
}


def get_all_states(
    include_atac: bool = True,
    include_epigenomic: bool = False,
    include_protein: bool = False,
    include_spatial: bool = False,
) -> dict:
    """Return the set of applicable states based on available modalities."""
    states = dict(RA_STATES)
    if include_epigenomic:
        states.update(EPI_ACTIVATING_STATES)
        states.update(EPI_REPRESSIVE_STATES)
    if include_protein:
        states.update(PROTEIN_STATES)
    if include_spatial:
        states.update(SPATIAL_STATES)
    return states
