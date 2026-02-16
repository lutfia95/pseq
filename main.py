from __future__ import annotations

import argparse
from typing import Dict, List, Optional, Tuple, Iterable

from src.arch import print_arch
from src.fasta_record import FastaRecord
from src.match_row import MatchRow
from src.decoy_matcher import DecoyMatcher



def read_fasta(path: str) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    header: Optional[str] = None # somehow Optional is giving errors on mac M2?! 
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(chunks).replace(" ", "").replace("\t", "").upper()
                    records.append(FastaRecord(header=header[1:].strip(), seq=seq))
                header = line
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seq = "".join(chunks).replace(" ", "").replace("\t", "").upper()
            records.append(FastaRecord(header=header[1:].strip(), seq=seq))

    return records


def write_fasta_pairs(path: str, pairs: List[Tuple[FastaRecord, FastaRecord]]) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for t, d in pairs:
            out.write(f">{t.header}\n")
            out.write(f"{t.seq}\n")
            out.write(f">{d.header}\n")
            out.write(f"{d.seq}\n")


def write_summary(path: str, rows: List[MatchRow], unmatched: List[FastaRecord]) -> None:
    total_targets = len(rows) + len(unmatched)
    matched = len(rows)

    tot_t_len = sum(r.target_len for r in rows)
    tot_d_len = sum(r.decoy_len for r in rows)

    tot_t_k = sum(r.target_k for r in rows)
    tot_t_r = sum(r.target_r for r in rows)
    tot_d_k = sum(r.decoy_k for r in rows)
    tot_d_r = sum(r.decoy_r for r in rows)

    with open(path, "w", encoding="utf-8") as out:
        out.write("# overall\n")
        out.write("total_targets\tmatched\tunmatched\t")
        out.write("total_target_len\ttotal_decoy_len\t")
        out.write("total_target_K\ttotal_target_R\t")
        out.write("total_decoy_K\ttotal_decoy_R\n")
        out.write(
            f"{total_targets}\t{matched}\t{len(unmatched)}\t"
            f"{tot_t_len}\t{tot_d_len}\t"
            f"{tot_t_k}\t{tot_t_r}\t"
            f"{tot_d_k}\t{tot_d_r}\n"
        )

        out.write("\n# per_target\n")
        out.write(
            "target_name\tdecoy_name\t"
            "target_len\tdecoy_len\tlen_diff\t"
            "target_K\ttarget_R\ttarget_KR\t"
            "decoy_K\tdecoy_R\tdecoy_KR\tkr_diff\t"
            "score\n"
        )
        for r in rows:
            out.write(
                f"{r.target_name}\t{r.decoy_name}\t"
                f"{r.target_len}\t{r.decoy_len}\t{r.len_diff}\t"
                f"{r.target_k}\t{r.target_r}\t{r.target_kr}\t"
                f"{r.decoy_k}\t{r.decoy_r}\t{r.decoy_kr}\t{r.kr_diff}\t"
                f"{r.score}\n"
            )

        if unmatched:
            out.write("\n# unmatched_targets\n")
            out.write("target_name\ttarget_len\ttarget_K\ttarget_R\ttarget_KR\n")
            for t in unmatched:
                seq = t.seq
                tL = len(seq)
                tK = sum(1 for c in seq if c == "K")
                tR = sum(1 for c in seq if c == "R")
                out.write(f"{t.name}\t{tL}\t{tK}\t{tR}\t{tK+tR}\n")


def cmd_generate_decoy(args):
    print(f'[WARNING] This is only your system check;) \n ')
    print_arch()

    print(f'[INFO] Reading targets fasta:\n {args.targets}')
    targets = read_fasta(args.targets)
    print(f'[INFO] Finished reading targets fasta:\n {args.targets}')

    print(f'[INFO] Reading decoys fasta:\n {args.decoys}')
    decoys = read_fasta(args.decoys)
    print(f'[INFO] Finished reading decoys fasta:\n {args.decoys}')

    matcher = DecoyMatcher(
        decoys=decoys,
        len_weight=args.len_weight,
        kr_weight=args.kr_weight,
        max_len_window=args.max_len_window,
        max_kr_window=args.max_kr_window,
    )

    print(f'[INFO] Matching...')
    pairs, rows, unmatched = matcher.match_all(targets)
    print(f'[INFO] Finished target matching!')

    print(f'[INFO] Writing {len(pairs)} to {args.out_fasta}')
    write_fasta_pairs(args.out_fasta, pairs)
    write_summary(args.out_summary, rows, unmatched)
    print(f'[INFO] Finished writing {args.out_fasta} and {args.out_summary}!')
    


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser(
        "generate_decoy",
        help="Match each target protein to a unique decoy with similar length and similar (K+R) content.",
    )
    g.add_argument("--targets", required=True, help="Target FASTA")
    g.add_argument("--decoys", required=True, help="Decoy FASTA (real decoys)")
    g.add_argument("--out-fasta", required=True, help="Output FASTA with interleaved target/decoy records")
    g.add_argument("--out-summary", required=True, help="Output TSV summary")
    g.add_argument("--len-weight", type=int, default=1, help="Weight for length difference in scoring")
    g.add_argument("--kr-weight", type=int, default=10, help="Weight for (K+R) difference in scoring")
    g.add_argument("--max-len-window", type=int, default=2000, help="Max length search window (aa)")
    g.add_argument("--max-kr-window", type=int, default=2000, help="Max (K+R) search window (count)")
    g.set_defaults(func=cmd_generate_decoy)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

# python main.py generate_decoy --targets ./data/targets.fasta --decoys ./data/ecoli.fasta --out-fasta ./output/mixed.fasta --out-summary ./output/summary.tsv
if __name__ == "__main__":
    main()