"""MoDES: Multi-Omics Discordance/Event State inference."""

from modes.core import MoDES, MoDESResult
from modes.data import MoDEData
from modes.events import EventCandidateBuilder

__all__ = [
    "MoDES",
    "MoDESResult",
    "MoDEData",
    "EventCandidateBuilder",
]
__version__ = "2.0.0"
