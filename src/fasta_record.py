
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class FastaRecord:
    header: str #biopy? as no need for sequence analysis, we will stay with lazy reading ;) 
    seq: str

    @property
    def name(self) -> str:
        return self.header.split()[0]

    def length(self) -> int:
        return len(self.seq)

    def count_aa(self, aa: str) -> int:
        aa = aa.upper()
        return sum(1 for c in self.seq if c == aa)