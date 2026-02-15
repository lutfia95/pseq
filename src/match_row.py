from __future__ import annotations
from dataclasses import dataclass
@dataclass
class MatchRow:
    target_name: str # just holder to save the real name
    decoy_name: str # same! 
    target_len: int # sequence length of the target, this is only for plotting later ;) 
    decoy_len: int # same! 
    target_k: int # number of k in target
    target_r: int # number of r in target
    decoy_k: int # same! 
    decoy_r: int # same! 
    target_kr: int # same! 
    decoy_kr: int # same! 
    len_diff: int # length differences, only for plotting! ;) 
    kr_diff: int # same 
    score: int # hand-made assigned score