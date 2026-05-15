"""Internal type definitions for MoDES."""

from dataclasses import dataclass, field


@dataclass
class EventCandidate:
    """A single candidate regulatory event linking a peak to a gene."""

    event_id: str
    gene: str
    peak_id: str
    peak_chr: str
    peak_start: int
    peak_end: int
    tf_name: str | None = None
    context: str = ""
    link_source: str = ""
    distance_to_tss: int = 0


@dataclass
class ModalityEffect:
    """Effect size estimate for one modality."""

    coef: float
    se: float
    z_score: float
    p_value: float
    fdr: float = 1.0
    direction: int = 0
    convergence: bool = True
    model_summary: dict = field(default_factory=dict)


@dataclass
class EventEvidence:
    """Evidence vector D_e for a single event."""

    event_id: str
    context: str
    effect_atac: ModalityEffect
    effect_rna: ModalityEffect
    effect_rna_given_atac: ModalityEffect | None = None
    quality_score: float = 0.5

    @property
    def z_atac(self) -> float:
        return self.effect_atac.z_score

    @property
    def z_rna(self) -> float:
        return self.effect_rna.z_score

    @property
    def z_rna_given_atac(self) -> float:
        if self.effect_rna_given_atac is None:
            return 0.0
        return self.effect_rna_given_atac.z_score

    @property
    def evidence_vector(self) -> list[float]:
        return [self.z_atac, self.z_rna, self.z_rna_given_atac, self.quality_score]


@dataclass
class EventState:
    """Classified regulatory state for an event."""

    event_id: str
    context: str
    state: str
    state_confidence: float = 1.0

    @property
    def local_fdr(self) -> float:
        return 1.0 - self.state_confidence


@dataclass
class EventResult:
    """Combined result for a single event -- the main output row."""

    event_id: str
    tf_name: str | None
    peak_id: str
    gene: str
    context: str
    atac_coef: float
    atac_se: float
    atac_pval: float
    atac_fdr: float
    atac_direction: int
    rna_coef: float
    rna_se: float
    rna_pval: float
    rna_fdr: float
    rna_direction: int
    rna_after_atac_coef: float
    rna_after_atac_se: float
    rna_after_atac_pval: float
    rna_after_atac_fdr: float
    state: str
    state_confidence: float
    quality_score: float
    artifact_risk: str = "low"
    artifact_reason: str = ""
    event_pval: float = 1.0
    event_fdr: float = 1.0
