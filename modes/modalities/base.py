"""Base modality specification for multi-modal MoDES."""

from dataclasses import dataclass, field


@dataclass
class ModalitySpec:
    """Specification for one molecular modality/layer."""

    name: str
    assay: str  # "RNA", "ATAC", "CUTTAG", "CUTRUN", "CHIPSEQ", "PROTEIN"
    feature_type: str  # "gene", "region", "protein", "mark"
    target: str | None = None  # "H3K27ac", "CTCF", "CD4"
    regulatory_role: str = "unknown"
    expected_rna_direction: int | None = None  # +1: activating, -1: repressive
    peak_type: str | None = None  # "narrow", "broad"
    normalization: str = "library_size"
    control: str | None = None  # IgG, input, spike-in
    priority: int = 0  # ordering in evidence vector

    def is_epigenomic(self) -> bool:
        return self.assay in ("CUTTAG", "CUTRUN", "CHIPSEQ")

    def is_activating(self) -> bool:
        return self.expected_rna_direction == 1

    def is_repressive(self) -> bool:
        return self.expected_rna_direction == -1


# Built-in RNA and ATAC specs
RNA_SPEC = ModalitySpec(
    name="rna", assay="RNA", feature_type="gene",
    regulatory_role="transcript_output", priority=10,
)

ATAC_SPEC = ModalitySpec(
    name="atac", assay="ATAC", feature_type="region",
    regulatory_role="chromatin_accessibility", priority=20,
    expected_rna_direction=1,
)
