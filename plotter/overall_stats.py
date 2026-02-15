from __future__ import annotations
from dataclasses import dataclass


@dataclass
class OverallStats:
    total_targets: int
    matched: int
    unmatched: int
    total_target_len: int
    total_decoy_len: int
    total_target_k: int
    total_target_r: int
    total_decoy_k: int
    total_decoy_r: int