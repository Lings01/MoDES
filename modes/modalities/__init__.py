"""MoDES modality abstraction layer — v2.0 multi-modal support."""

from modes.modalities.base import ModalitySpec
from modes.modalities.cuttag import (
    CUTTAG_REGISTRY,
    get_cuttag_target_info,
    validate_cuttag_features,
    make_cuttag_spec,
)

__all__ = [
    "ModalitySpec",
    "CUTTAG_REGISTRY",
    "get_cuttag_target_info",
    "validate_cuttag_features",
    "make_cuttag_spec",
]
