from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Iterable
from .fasta_record import FastaRecord
from .match_row import MatchRow

class DecoyMatcher:
    def __init__(
        self,
        decoys: List[FastaRecord],
        len_weight: int = 1,
        kr_weight: int = 10,
        max_len_window: int = 2000,
        max_kr_window: int = 2000,
    ):
        self.decoys = decoys
        self.len_weight = int(len_weight)
        self.kr_weight = int(kr_weight)
        self.max_len_window = int(max_len_window)
        self.max_kr_window = int(max_kr_window)

        self._by_len_kr: Dict[Tuple[int, int], List[int]] = {}
        self._used: List[bool] = [False] * len(decoys)

        for i, rec in enumerate(decoys):
            L = len(rec.seq)
            kr = self._kr_count(rec.seq)
            self._by_len_kr.setdefault((L, kr), []).append(i)

    @staticmethod # function connected to the class but not need for instance data, same error M2? def _...(self, ) solved the issue! 
    def _kr_count(seq: str) -> int:
        return sum(1 for c in seq if c == "K" or c == "R")

    @staticmethod
    def _k_count(seq: str) -> int:
        return sum(1 for c in seq if c == "K")

    @staticmethod
    def _r_count(seq: str) -> int:
        return sum(1 for c in seq if c == "R")

    def _pop_first_unused(self, idxs: List[int]) -> Optional[int]:
        while idxs:
            i = idxs.pop()
            if not self._used[i]:
                self._used[i] = True
                return i
        return None

    def _iter_search_shell(
        self, L: int, kr: int, dL: int, dKR: int
    ) -> Iterable[Tuple[int, int]]:
        if dL == 0 and dKR == 0:
            yield (L, kr)
            return

        lens = [L] if dL == 0 else [L - dL, L + dL] if dL > 0 else []
        krs = [kr] if dKR == 0 else [kr - dKR, kr + dKR] if dKR > 0 else []

        for LL in lens:
            for KK in krs:
                yield (LL, KK)

        if dL > 0 and dKR > 0:
            yield (L - dL, kr)
            yield (L + dL, kr)
            yield (L, kr - dKR)
            yield (L, kr + dKR)
            
    # ToDo: fix if human proteins, however too many for loops!  O(L .K .(1 + B)) = O(L.K.B)
    def find_best_match(self, target: FastaRecord) -> Optional[Tuple[int, int, int, int]]:
        tL = len(target.seq)
        tkr = self._kr_count(target.seq)

        best_i = None
        best_score = None
        best_dL = None
        best_dKR = None

        max_dL = self.max_len_window
        max_dKR = self.max_kr_window
        # assume L = max_len_window = max_dL && K = max_kr_window = max_dKR
        # worst-case run time per target: O(max(L, K))
        for radius in range(0, max(max_dL, max_dKR) + 1):
            dL_min = max(0, radius - max_dKR)
            dL_max = min(max_dL, radius)
            for dL in range(dL_min, dL_max + 1):
                dKR = radius - dL
                if dKR < 0 or dKR > max_dKR:
                    continue

                for key in self._iter_search_shell(tL, tkr, dL, dKR): # O(1)
                    idxs = self._by_len_kr.get(key)
                    if not idxs:
                        continue

                    i = None
                    for j in reversed(idxs):
                        if not self._used[j]:
                            i = j
                            break
                    if i is None:
                        continue

                    score = self.len_weight * abs(key[0] - tL) + self.kr_weight * abs(
                        key[1] - tkr
                    )

                    if best_score is None or score < best_score:
                        best_score = score
                        best_i = i
                        best_dL = abs(key[0] - tL)
                        best_dKR = abs(key[1] - tkr)

            if best_i is not None and radius > (best_score // max(1, min(self.len_weight, self.kr_weight))):
                break

        if best_i is None:
            return None

        self._used[best_i] = True
        return best_i, best_score or 0, best_dL or 0, best_dKR or 0

    def match_all(self, targets: List[FastaRecord]) -> Tuple[List[Tuple[FastaRecord, FastaRecord]], List[MatchRow], List[FastaRecord]]:
        pairs: List[Tuple[FastaRecord, FastaRecord]] = []
        rows: List[MatchRow] = []
        unmatched: List[FastaRecord] = []

        for t in targets:
            res = self.find_best_match(t)
            if res is None:
                unmatched.append(t)
                continue

            i, score, _, _ = res
            d = self.decoys[i]

            tL = len(t.seq)
            dL = len(d.seq)

            tK = self._k_count(t.seq)
            tR = self._r_count(t.seq)
            dK = self._k_count(d.seq)
            dR = self._r_count(d.seq)

            tKR = tK + tR
            dKR = dK + dR

            row = MatchRow(
                target_name=t.name,
                decoy_name=d.name,
                target_len=tL,
                decoy_len=dL,
                target_k=tK,
                target_r=tR,
                decoy_k=dK,
                decoy_r=dR,
                target_kr=tKR,
                decoy_kr=dKR,
                len_diff=abs(tL - dL),
                kr_diff=abs(tKR - dKR),
                score=score,
            )

            pairs.append((t, d))
            rows.append(row)

        return pairs, rows, unmatched
