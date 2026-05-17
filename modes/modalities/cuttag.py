"""CUT&Tag / CUT&RUN / ChIP-seq target registry and helpers."""

from __future__ import annotations

from modes.modalities.base import ModalitySpec

# Target registry: epigenomic mark → biological semantics
CUTTAG_REGISTRY: dict[str, dict] = {
    "H3K27ac": {
        "role": "activating_enhancer",
        "expected_rna_direction": 1,
        "peak_type": "narrow_or_broad",
    },
    "H3K4me1": {
        "role": "enhancer_priming",
        "expected_rna_direction": 1,
        "peak_type": "broad",
    },
    "H3K4me3": {
        "role": "active_promoter",
        "expected_rna_direction": 1,
        "peak_type": "narrow",
    },
    "H3K27me3": {
        "role": "polycomb_repression",
        "expected_rna_direction": -1,
        "peak_type": "broad",
    },
    "H3K9me3": {
        "role": "heterochromatin_repression",
        "expected_rna_direction": -1,
        "peak_type": "broad",
    },
    "H3K36me3": {
        "role": "transcription_elongation",
        "expected_rna_direction": 1,
        "peak_type": "broad",
    },
    "CTCF": {
        "role": "insulator_or_architectural",
        "expected_rna_direction": None,
        "peak_type": "narrow",
    },
    "RAD21": {
        "role": "cohesin_component",
        "expected_rna_direction": None,
        "peak_type": "narrow",
    },
    # Generic TF (overridden by user-provided role)
    "TF": {
        "role": "transcription_factor",
        "expected_rna_direction": None,
        "peak_type": "narrow",
    },
}


def get_cuttag_target_info(target: str) -> dict:
    """Look up registry info for a CUT&Tag target."""
    if target in CUTTAG_REGISTRY:
        return dict(CUTTAG_REGISTRY[target])
    # Fallback for unknown targets
    return {
        "role": "unknown",
        "expected_rna_direction": None,
        "peak_type": "unknown",
    }


def make_cuttag_spec(
    name: str,
    target: str,
    assay: str = "CUTTAG",
) -> ModalitySpec:
    """Create a ModalitySpec from a CUT&Tag target name."""
    info = get_cuttag_target_info(target)
    return ModalitySpec(
        name=name,
        assay=assay.upper(),
        feature_type="region",
        target=target,
        regulatory_role=info["role"],
        expected_rna_direction=info.get("expected_rna_direction"),
        peak_type=info.get("peak_type"),
        priority=30,
    )


def validate_cuttag_features(features_df) -> list[str]:
    """Validate a CUT&Tag features DataFrame. Returns list of issues."""
    issues = []
    required = ["feature_id", "chr", "start", "end", "assay", "target"]
    for col in required:
        if col not in features_df.columns:
            issues.append(f"Missing required column: {col}")
    if "target" in features_df.columns:
        unknown = set(features_df["target"].unique()) - set(CUTTAG_REGISTRY.keys()) - {"TF"}
        if unknown:
            issues.append(
                f"Unrecognized targets (will use fallback): {sorted(unknown)}"
            )
    return issues
