"""StateRule grammar: declarative evidence rules for multi-modal state classification.

Each state is defined by RequiredEvidence (must be significant in a given direction),
NeutralEvidence (allowed but not required), and ForbiddenEvidence (must NOT be present).

The StateClassifier scores all applicable rules simultaneously and selects the best-matching
state by assignment_score, rather than using priority-based if/else chains.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass(frozen=True)
class RequiredEvidence:
    """Evidence that MUST be significant in the specified direction for this state."""
    modality: str
    direction: int  # +1 (up), -1 (down)
    role: str | None = None
    target: str | None = None


@dataclass(frozen=True)
class NeutralEvidence:
    """Evidence that is allowed but not required; absence is not penalized."""
    modality: str
    role: str | None = None


@dataclass(frozen=True)
class ForbiddenEvidence:
    """Evidence that MUST NOT be significant in the specified direction."""
    modality: str
    direction: int
    role: str | None = None


@dataclass(frozen=True)
class StateRule:
    """Declarative rule defining the evidence pattern for one biological state."""
    name: str
    required: Sequence[RequiredEvidence]
    neutral: Sequence[NeutralEvidence] = ()
    forbidden: Sequence[ForbiddenEvidence] = ()
    description: str = ""
    interpretation_strength: str = "association"  # association | hypothesis


# ── RNA+ATAC core states ──────────────────────────────────────────────

CONCORDANT = StateRule(
    name="concordant",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
    ],
    description="ATAC and RNA change concordantly in the same direction.",
)

CHROMATIN_PRIMED = StateRule(
    name="chromatin_primed",
    required=[RequiredEvidence("atac", +1)],
    neutral=[NeutralEvidence("rna")],
    description="ATAC changes while RNA does not show significant change.",
)

RNA_ONLY = StateRule(
    name="rna_only",
    required=[RequiredEvidence("rna", +1)],
    neutral=[NeutralEvidence("atac")],
    description="RNA changes without corresponding local chromatin change.",
)

DISCORDANT_OPPOSITE = StateRule(
    name="discordant_opposite",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", -1),
    ],
    description="ATAC and RNA change in opposite directions.",
)

NULL = StateRule(
    name="null",
    required=[],
    neutral=[
        NeutralEvidence("atac"),
        NeutralEvidence("rna"),
    ],
    description="No significant change detected for this event under the tested contrast.",
)

# ── CUT&Tag activating mark states ────────────────────────────────────

EPIGENOMIC_CONCORDANT = StateRule(
    name="epigenomic_concordant",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
        RequiredEvidence("cuttag_activating", +1, role="activating_mark"),
    ],
    description="ATAC, RNA, and activating histone mark change concordantly.",
)

ACTIVE_ENHANCER_PRIMED = StateRule(
    name="active_enhancer_primed",
    required=[RequiredEvidence("cuttag_activating", +1, role="activating_mark")],
    neutral=[NeutralEvidence("rna")],
    description="Activating histone mark present; transcription not yet responding.",
)

MARK_ONLY = StateRule(
    name="mark_only",
    required=[RequiredEvidence("cuttag_activating", +1, role="activating_mark")],
    neutral=[
        NeutralEvidence("atac"),
        NeutralEvidence("rna"),
    ],
    description="Histone mark change without detectable chromatin or RNA change.",
)

# ── CUT&Tag repressive mark states ────────────────────────────────────

REPRESSIVE_CONCORDANT = StateRule(
    name="repressive_concordant",
    required=[
        RequiredEvidence("cuttag_repressive", +1, role="repressive_mark"),
        RequiredEvidence("rna", -1),
    ],
    description="Repressive histone mark gain and RNA decrease.",
)

DEREPRESSION = StateRule(
    name="derepression",
    required=[
        RequiredEvidence("cuttag_repressive", -1, role="repressive_mark"),
        RequiredEvidence("rna", +1),
    ],
    description="Loss of repressive mark accompanied by RNA increase.",
)

REPRESSIVE_PRIMED = StateRule(
    name="repressive_primed",
    required=[RequiredEvidence("cuttag_repressive", +1, role="repressive_mark")],
    neutral=[NeutralEvidence("rna")],
    description="Repressive histone mark gain; RNA not yet responding.",
)

# ── Protein states ────────────────────────────────────────────────────

FULL_ACTIVATION = StateRule(
    name="full_activation",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
        RequiredEvidence("protein", +1),
    ],
    description="ATAC, RNA, and protein layers show concordant differential signal.",
)

PROTEIN_BUFFERED = StateRule(
    name="protein_buffered",
    required=[RequiredEvidence("rna", +1)],
    neutral=[NeutralEvidence("protein")],
    description="RNA changes but protein level does not, suggesting post-transcriptional buffering.",
)

PROTEIN_MEMORY = StateRule(
    name="protein_memory",
    required=[RequiredEvidence("protein", +1)],
    neutral=[NeutralEvidence("rna")],
    description="Protein level changes persist while RNA returns to baseline.",
)

PROTEIN_OPPOSITE = StateRule(
    name="protein_opposite",
    required=[
        RequiredEvidence("rna", +1),
        RequiredEvidence("protein", -1),
    ],
    description="RNA and protein changes in opposite directions.",
)

# ── Spatial states ────────────────────────────────────────────────────

SPATIAL_REGION_SPECIFIC = StateRule(
    name="spatial_region_specific",
    required=[RequiredEvidence("spatial", +1, role="spatial_autocorrelation")],
    neutral=[NeutralEvidence("rna"), NeutralEvidence("atac")],
    description="Event shows spatial autocorrelation suggesting region-specific regulation.",
)

SPATIAL_NICHE_DRIVEN = StateRule(
    name="spatial_niche_driven",
    required=[
        RequiredEvidence("spatial", +1, role="spatial_autocorrelation"),
        RequiredEvidence("spatial", +1, role="neighbor_effect"),
    ],
    description="Event is strongly driven by local niche rather than cell-intrinsic signals.",
)

CELL_INTRINSIC = StateRule(
    name="cell_intrinsic",
    required=[],
    neutral=[
        NeutralEvidence("rna"),
        NeutralEvidence("atac"),
        NeutralEvidence("spatial"),
    ],
    description="No spatial dependence detected; event appears cell-intrinsic.",
)

SPATIAL_EDGE_ARTIFACT = StateRule(
    name="spatial_edge_artifact",
    required=[RequiredEvidence("spatial", +1, role="edge_artifact")],
    description="Event is concentrated at tissue edges; likely technical artifact.",
)

# ── Ambiguous / unresolved ────────────────────────────────────────────

MIXED_EVIDENCE = StateRule(
    name="mixed_evidence",
    required=[],
    description="Multiple states scored similarly; evidence pattern is ambiguous.",
)

UNRESOLVED = StateRule(
    name="unresolved",
    required=[],
    description="Insufficient or invalid evidence to assign any state.",
)

# ── Rule registries ───────────────────────────────────────────────────

# Ordered by specificity: most specific (most required evidence) scored first
RA_RULES: list[StateRule] = [CONCORDANT, DISCORDANT_OPPOSITE, CHROMATIN_PRIMED, RNA_ONLY, NULL]

EPIGENOMIC_RULES: list[StateRule] = [
    EPIGENOMIC_CONCORDANT, ACTIVE_ENHANCER_PRIMED, MARK_ONLY,
    REPRESSIVE_CONCORDANT, DEREPRESSION, REPRESSIVE_PRIMED,
]

PROTEIN_RULES: list[StateRule] = [
    FULL_ACTIVATION, PROTEIN_BUFFERED, PROTEIN_MEMORY, PROTEIN_OPPOSITE,
]

SPATIAL_RULES: list[StateRule] = [
    SPATIAL_REGION_SPECIFIC, SPATIAL_NICHE_DRIVEN, CELL_INTRINSIC, SPATIAL_EDGE_ARTIFACT,
]


def get_all_rules(
    include_epigenomic: bool = False,
    include_protein: bool = False,
    include_spatial: bool = False,
) -> list[StateRule]:
    """Return all applicable StateRules given the available modalities."""
    rules = list(RA_RULES)
    if include_epigenomic:
        rules = list(EPIGENOMIC_RULES) + rules
    if include_protein:
        rules = list(PROTEIN_RULES) + rules
    if include_spatial:
        rules = list(SPATIAL_RULES) + rules
    return rules


ALL_STATE_NAMES: set[str] = {
    r.name for rules in [RA_RULES, EPIGENOMIC_RULES, PROTEIN_RULES, SPATIAL_RULES]
    for r in rules
} | {"mixed_evidence", "unresolved"}


# ── Directed p-value helper ───────────────────────────────────────────

def directed_pvalue(pval: float, coef: float, expected_direction: int) -> float:
    """Compute a directional evidence score from a two-sided p-value.

    This is NOT a one-sided test; it is a directional evidence score that
    down-weights p-values when the effect direction matches expectation,
    and returns 1.0 when the direction is opposite.

    Parameters
    ----------
    pval : float
        Two-sided p-value from the GLM.
    coef : float
        Estimated coefficient (log fold change).
    expected_direction : int
        +1 (expected up), -1 (expected down), or 0 (any direction).

    Returns
    -------
    float
        Directional evidence score in [pval/2, 1.0].
    """
    if expected_direction == 0:
        return pval
    if np.sign(coef) == expected_direction:
        return min(pval / 2.0, 1.0)
    return 1.0


import numpy as np
