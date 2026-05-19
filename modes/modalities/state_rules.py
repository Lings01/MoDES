"""StateRule grammar: declarative evidence rules for multi-modal state annotation.

v2.0: Full evidence semantics with RequiredEvidence (must be sig+direction),
RequiredAbsentEvidence (must be measured but NOT sig), ForbiddenEvidence
(must NOT be sig in given direction), OptionalEvidence (bonus if present),
and MissingPolicy (what happens when a modality is not measured).

Each StateRule carries a state_family (e.g., "concordant") and a specific
state name (e.g., "concordant_activation").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass(frozen=True)
class RequiredEvidence:
    """Modality must be significant in the specified direction."""
    modality: str
    direction: int  # +1 (up), -1 (down)
    role: str | None = None
    target: str | None = None


@dataclass(frozen=True)
class RequiredAbsentEvidence:
    """Modality must be measured (available) but NOT statistically significant."""
    modality: str
    role: str | None = None
    target: str | None = None
    require_available: bool = True  # if True, modality must exist in data


@dataclass(frozen=True)
class ForbiddenEvidence:
    """Modality must NOT be significant, or must NOT be significant in a given direction."""
    modality: str
    direction: int | None = None  # None = any direction is forbidden
    role: str | None = None
    target: str | None = None


@dataclass(frozen=True)
class OptionalEvidence:
    """Bonus if significant in the specified direction; no penalty if absent."""
    modality: str
    direction: int | None = None
    role: str | None = None
    target: str | None = None
    bonus: float = 0.1  # multiplicative bonus factor


@dataclass(frozen=True)
class MissingPolicy:
    """Policy for when a modality is not measured in the data."""
    modality: str
    allowed: bool = False  # if False, rule cannot be triggered
    penalty: float = 0.5  # multiplicative penalty if allowed


@dataclass(frozen=True)
class StateRule:
    """Declarative rule defining multi-modal evidence pattern for one state."""
    name: str
    required: Sequence[RequiredEvidence] = ()
    required_absent: Sequence[RequiredAbsentEvidence] = ()
    forbidden: Sequence[ForbiddenEvidence] = ()
    optional: Sequence[OptionalEvidence] = ()
    missing_policy: Sequence[MissingPolicy] = ()
    state_family: str = ""
    description: str = ""
    interpretation_strength: str = "association"


# ── RNA+ATAC core: full-direction states ─────────────────────────────

CONCORDANT_ACTIVATION = StateRule(
    name="concordant_activation",
    state_family="concordant",
    required=[RequiredEvidence("atac", +1), RequiredEvidence("rna", +1)],
    description="ATAC and RNA both increase under the tested contrast.",
)

CONCORDANT_REPRESSION = StateRule(
    name="concordant_repression",
    state_family="concordant",
    required=[RequiredEvidence("atac", -1), RequiredEvidence("rna", -1)],
    description="ATAC and RNA both decrease under the tested contrast.",
)

DISCORDANT_OPENING_REPRESSION = StateRule(
    name="discordant_opening_repression",
    state_family="discordant",
    required=[RequiredEvidence("atac", +1), RequiredEvidence("rna", -1)],
    description="ATAC increases while RNA decreases.",
)

DISCORDANT_CLOSING_ACTIVATION = StateRule(
    name="discordant_closing_activation",
    state_family="discordant",
    required=[RequiredEvidence("atac", -1), RequiredEvidence("rna", +1)],
    description="ATAC decreases while RNA increases.",
)

CHROMATIN_OPEN_PRIMED = StateRule(
    name="chromatin_open_primed",
    state_family="chromatin_primed",
    required=[RequiredEvidence("atac", +1)],
    required_absent=[RequiredAbsentEvidence("rna")],
    description="ATAC increases without corresponding RNA change.",
)

CHROMATIN_CLOSED_PRIMED = StateRule(
    name="chromatin_closed_primed",
    state_family="chromatin_primed",
    required=[RequiredEvidence("atac", -1)],
    required_absent=[RequiredAbsentEvidence("rna")],
    description="ATAC decreases without corresponding RNA change.",
)

RNA_UP_ONLY = StateRule(
    name="rna_up_only",
    state_family="rna_only",
    required=[RequiredEvidence("rna", +1)],
    required_absent=[RequiredAbsentEvidence("atac")],
    description="RNA increases without local chromatin change.",
)

RNA_DOWN_ONLY = StateRule(
    name="rna_down_only",
    state_family="rna_only",
    required=[RequiredEvidence("rna", -1)],
    required_absent=[RequiredAbsentEvidence("atac")],
    description="RNA decreases without local chromatin change.",
)

NULL = StateRule(
    name="null",
    state_family="null",
    required_absent=[RequiredAbsentEvidence("atac"), RequiredAbsentEvidence("rna")],
    description="No significant change detected.",
)

# ── CUT&Tag activating marks ────────────────────────────────────────

EPIGENOMIC_CONCORDANT_ACTIVATION = StateRule(
    name="epigenomic_concordant_activation",
    state_family="epigenomic_concordant",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
        RequiredEvidence("cuttag_activating", +1, role="activating_mark"),
    ],
    description="ATAC, RNA, and activating mark increase concordantly.",
)

ACTIVE_ENHANCER_PRIMED = StateRule(
    name="active_enhancer_primed",
    state_family="epigenomic_primed",
    required=[RequiredEvidence("cuttag_activating", +1, role="activating_mark")],
    required_absent=[RequiredAbsentEvidence("rna")],
    optional=[OptionalEvidence("atac", +1)],
    description="Activating mark present; RNA not responding. ATAC optional.",
)

MARK_ONLY = StateRule(
    name="mark_only",
    state_family="epigenomic_only",
    required=[RequiredEvidence("cuttag_activating", +1, role="activating_mark")],
    required_absent=[
        RequiredAbsentEvidence("rna"),
        RequiredAbsentEvidence("atac"),
    ],
    description="Histone mark change without detectable chromatin or RNA change.",
)

# ── CUT&Tag repressive marks ────────────────────────────────────────

REPRESSIVE_CONCORDANT = StateRule(
    name="repressive_concordant",
    state_family="epigenomic_repressive",
    required=[
        RequiredEvidence("cuttag_repressive", +1, role="repressive_mark"),
        RequiredEvidence("rna", -1),
    ],
    description="Repressive mark gain with RNA decrease.",
)

DEREPRESSION = StateRule(
    name="derepression",
    state_family="epigenomic_derepression",
    required=[
        RequiredEvidence("cuttag_repressive", -1, role="repressive_mark"),
        RequiredEvidence("rna", +1),
    ],
    description="Loss of repressive mark with RNA increase.",
)

REPRESSIVE_PRIMED = StateRule(
    name="repressive_primed",
    state_family="epigenomic_primed",
    required=[RequiredEvidence("cuttag_repressive", +1, role="repressive_mark")],
    required_absent=[RequiredAbsentEvidence("rna")],
    description="Repressive mark gain; RNA not responding.",
)

# ── Protein states ──────────────────────────────────────────────────

FULL_ACTIVATION = StateRule(
    name="full_activation_up",
    state_family="full_activation",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
        RequiredEvidence("protein", +1),
    ],
    description="ATAC, RNA, and protein show concordant increase.",
)

PROTEIN_BUFFERED_UP = StateRule(
    name="protein_buffered_up",
    state_family="protein_buffered",
    required=[RequiredEvidence("rna", +1)],
    required_absent=[RequiredAbsentEvidence("protein", require_available=True)],
    description="RNA increases but protein does not change (post-transcriptional discordance).",
)

PROTEIN_BUFFERED_DOWN = StateRule(
    name="protein_buffered_down",
    state_family="protein_buffered",
    required=[RequiredEvidence("rna", -1)],
    required_absent=[RequiredAbsentEvidence("protein", require_available=True)],
    description="RNA decreases but protein does not change.",
)

PROTEIN_MEMORY_UP = StateRule(
    name="protein_memory_up",
    state_family="protein_memory",
    required=[RequiredEvidence("protein", +1)],
    required_absent=[RequiredAbsentEvidence("rna", require_available=True)],
    description="Protein increases while RNA is at baseline.",
)

PROTEIN_MEMORY_DOWN = StateRule(
    name="protein_memory_down",
    state_family="protein_memory",
    required=[RequiredEvidence("protein", -1)],
    required_absent=[RequiredAbsentEvidence("rna", require_available=True)],
    description="Protein decreases while RNA is at baseline.",
)

PROTEIN_OPPOSITE = StateRule(
    name="protein_opposite_up_down",
    state_family="protein_opposite",
    required=[RequiredEvidence("rna", +1), RequiredEvidence("protein", -1)],
    description="RNA and protein change in opposite directions.",
)

# ── Spatial states ──────────────────────────────────────────────────

SPATIAL_REGION_SPECIFIC = StateRule(
    name="spatial_region_specific",
    state_family="spatial",
    required=[RequiredEvidence("spatial", +1, role="spatial_autocorrelation")],
    description="Event shows spatial autocorrelation.",
)

SPATIAL_NICHE_DRIVEN = StateRule(
    name="spatial_niche_driven",
    state_family="spatial",
    required=[
        RequiredEvidence("spatial", +1, role="spatial_autocorrelation"),
        RequiredEvidence("spatial", +1, role="neighbor_effect"),
    ],
    description="Event driven by local niche signals.",
)

CELL_INTRINSIC = StateRule(
    name="cell_intrinsic",
    state_family="spatial",
    required_absent=[RequiredAbsentEvidence("spatial")],
    description="No spatial dependence; event appears cell-intrinsic.",
)

SPATIAL_EDGE_ARTIFACT = StateRule(
    name="spatial_edge_artifact",
    state_family="spatial_artifact",
    required=[RequiredEvidence("spatial", +1, role="edge_artifact")],
    description="Event concentrated at tissue edges; likely technical artifact.",
)

# ── Ambiguous / unresolved ──────────────────────────────────────────

MIXED_EVIDENCE = StateRule(
    name="mixed_evidence", state_family="ambiguous",
    description="Multiple states scored similarly.",
)

UNRESOLVED = StateRule(
    name="unresolved", state_family="unresolved",
    description="Insufficient evidence to assign any state.",
)

# ── Rule registries ─────────────────────────────────────────────────

RA_RULES: list[StateRule] = [
    CONCORDANT_ACTIVATION, CONCORDANT_REPRESSION,
    DISCORDANT_OPENING_REPRESSION, DISCORDANT_CLOSING_ACTIVATION,
    CHROMATIN_OPEN_PRIMED, CHROMATIN_CLOSED_PRIMED,
    RNA_UP_ONLY, RNA_DOWN_ONLY,
    NULL,
]

EPIGENOMIC_RULES: list[StateRule] = [
    EPIGENOMIC_CONCORDANT_ACTIVATION,
    ACTIVE_ENHANCER_PRIMED, MARK_ONLY,
    REPRESSIVE_CONCORDANT, DEREPRESSION, REPRESSIVE_PRIMED,
]

PROTEIN_RULES: list[StateRule] = [
    FULL_ACTIVATION,
    PROTEIN_BUFFERED_UP, PROTEIN_BUFFERED_DOWN,
    PROTEIN_MEMORY_UP, PROTEIN_MEMORY_DOWN,
    PROTEIN_OPPOSITE,
]

SPATIAL_RULES: list[StateRule] = [
    SPATIAL_REGION_SPECIFIC, SPATIAL_NICHE_DRIVEN,
    CELL_INTRINSIC, SPATIAL_EDGE_ARTIFACT,
]


def get_all_rules(
    include_epigenomic: bool = False,
    include_protein: bool = False,
    include_spatial: bool = False,
) -> list[StateRule]:
    """Return all applicable StateRules given available modalities."""
    rules = list(RA_RULES)
    if include_epigenomic:
        rules = list(EPIGENOMIC_RULES) + rules
    if include_protein:
        rules = list(PROTEIN_RULES) + rules
    if include_spatial:
        rules = list(SPATIAL_RULES) + rules
    return rules


ALL_STATE_NAMES: set[str] = {
    r.name for rules_list in [RA_RULES, EPIGENOMIC_RULES, PROTEIN_RULES, SPATIAL_RULES]
    for r in rules_list
} | {"mixed_evidence", "unresolved"}

ALL_STATE_FAMILIES: set[str] = {
    r.state_family for rules_list in [RA_RULES, EPIGENOMIC_RULES, PROTEIN_RULES, SPATIAL_RULES]
    for r in rules_list if r.state_family
} | {"ambiguous", "unresolved"}


# ── Directed evidence score helper ──────────────────────────────────

import numpy as np


def directed_score(pval: float, coef: float, expected_direction: int) -> float:
    """Directional evidence score (NOT a one-sided p-value).

    Returns a score in [0, inf) where higher = stronger evidence:
      - When direction matches: -log10(pval/2)
      - When direction opposite: 0.0

    This is a ranking-oriented evidence score, NOT a calibrated p-value.
    """
    if expected_direction == 0:
        return max(-np.log10(max(pval, 1e-15)), 0.0)
    if np.sign(coef) == expected_direction:
        return max(-np.log10(max(pval / 2.0, 1e-15)), 0.0)
    return 0.0
