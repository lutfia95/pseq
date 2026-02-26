from __future__ import annotations
#python .\src\fl_sites.py --fasta .\FLAna\bacsu_sorc\bacsu_sorc\bsu_sorc.fasta --sites-dir .\FLAna\bacsu_sorc\bacsu_sorc\reports\ --outdir .\output\bacsu_sorc --log-y --style 2
import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

STYLE2_LEN_RE = re.compile(r"_([0-9]+)\b")
FASTA_BOUNDARY_RE = re.compile(r"&&\s*\d+\s*&&\s*(\d+)\s*$")
ACCESSION_RE = re.compile(r"\b(?:sp|tr)\|([A-Z0-9]+)\|", re.IGNORECASE)
POS_RE = re.compile(r"\((\d+)\)")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Classify site-level results as BACSU vs SORC-tail using concatenation boundary from FASTA headers."
    )
    p.add_argument(
        "--fasta",
        required=True,
        help="Path to the FASTA used for the search.",
    )
    p.add_argument(
        "--sites-dir",
        required=True,
        help="Directory containing pLink3 *sites.csv files (e.g., result_*.filtered_*_sites.csv).",
    )
    p.add_argument(
        "--outdir",
        default="sites_boundary_plots",
        help="Output directory for plots and summaries (default: sites_boundary_plots).",
    )
    p.add_argument(
        "--log-y",
        action="store_true",
        help="Use log scale for count plots.",
    )
    p.add_argument(
        "--style",
        type=int,
        choices=[1, 2],
        default=1,
        help=(
            "FASTA boundary style. "
            "1: headers end with '&& <num> && <boundary>'. "
            "2: headers encode lengths like 'YUKE_BACSU_97-A9FUC4_SORC5_107' and boundary=first length."
        ),
    )
    return p.parse_args()

def read_fasta_boundaries(fasta_path: Path, style: int) -> dict[str, int]:
    """
    style 1: boundary from header ending: '&& <n> && <boundary>'
    style 2: boundary from first *_<len> token in the FASTA identifier after the last '|'
             e.g. sp|C0SP85-A9FUC4|YUKE_BACSU_97-A9FUC4_SORC5_107  -> boundary=97
    """
    boundaries: dict[str, int] = {}

    with fasta_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            primary = header.split()[0]  # full FASTA ID token (no whitespace)

            boundary: int | None = None
            if style == 1:
                m = FASTA_BOUNDARY_RE.search(header)
                if m:
                    boundary = int(m.group(1))
            elif style == 2:
                # For style2, we expect something like:
                # sp|C0SP85-A9FUC4|YUKE_BACSU_97-A9FUC4_SORC5_107
                # We take the part after the last '|' and read the FIRST _<number> as boundary.
                after_pipe = primary.split("|")[-1]
                m = STYLE2_LEN_RE.search(after_pipe)
                if m:
                    boundary = int(m.group(1))

            if boundary is None:
                continue

            boundaries[primary] = boundary

            am = ACCESSION_RE.search(primary)
            if am:
                boundaries[am.group(1)] = boundary

    if not boundaries:
        raise ValueError(
            "No boundaries found in FASTA. "
            "For style=1 expected headers ending with '&& <num> && <boundary>'. "
            "For style=2 expected IDs like '...|NAME_BACSU_<len>-..._SORC..._<len>'."
        )
    return boundaries

def read_plink_sites_only(csv_path: Path) -> pd.DataFrame:
    
    keep: list[str] = []
    with csv_path.open("r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n")
        keep.append(header)
        for line in f:
            if not line.strip():
                continue
            if line.startswith(","):
                continue
            keep.append(line.rstrip("\n"))
    return pd.read_csv(StringIO("\n".join(keep)), sep=",", engine="python")


def infer_type_label(filename: str) -> str:
    s = filename.lower()
    if "cross-linked_sites" in s or "cross-linked" in s:
        return "cross-linked"
    if "mono-linked_sites" in s or "mono-linked" in s:
        return "mono-linked"
    if "loop-linked_sites" in s or "loop-linked" in s:
        return "loop-linked"
    if "regular_sites" in s or "linear_sites" in s or "regular" in s:
        return "linear"
    return "unknown-type"


@dataclass
class ParsedSite:
    proteins: list[str]          # protein keys (best-effort)
    positions: list[int]         # positions extracted from "(n)"
    per_side: list[str]          # "BACSU" / "SORC" / "UNKNOWN" per extracted position
    site_class: str              # "BACSU-only" / "SORC-only" / "MIXED" / "UNKNOWN"


PROT_SPLIT_RE = re.compile(r"-(?=(?:sp|tr)\|)", re.IGNORECASE)
TRAIL_PARENS_RE = re.compile(r"(?:\(\d+\))+$")  # remove one or many trailing (n)


def extract_protein_keys_and_positions(protein_field: object) -> tuple[list[str], list[int]]:
    if not isinstance(protein_field, str):
        return ([], [])
    s = protein_field.strip()
    if not s:
        return ([], [])

    parts = [p.strip() for p in PROT_SPLIT_RE.split(s) if p.strip()]

    proteins: list[str] = []
    for part in parts:
        primary = part.split()[0]
        primary = TRAIL_PARENS_RE.sub("", primary)
        proteins.append(primary)

        am = ACCESSION_RE.search(primary)
        if am:
            proteins.append(am.group(1))

    seen = set()
    proteins2: list[str] = []
    for p in proteins:
        if p not in seen:
            proteins2.append(p)
            seen.add(p)

    positions = [int(x) for x in POS_RE.findall(s)]
    return (proteins2, positions)

def classify_site(protein_field: object, boundaries: dict[str, int]) -> ParsedSite:
    proteins, positions = extract_protein_keys_and_positions(protein_field)

    if not proteins or not positions:
        return ParsedSite(proteins=proteins, positions=positions, per_side=[], site_class="UNKNOWN")

    # Map positions to sides. We must decide which boundary to use.
    # For cross-links, Protein field usually contains two proteins; for mono/loop, one.
    # We assign boundaries using the FIRST matching protein key we can find in boundaries.
    # If there are two proteins and one boundary missing, some sides bSORCme UNKNOWN.
    #
    # Strategy:
    # - If we can find 2+ distinct proteins with boundaries, use the first two boundary values in order.
    # - Otherwise use the first boundary we can find for all positions.
    bvals: list[int] = []
    for pk in proteins:
        if pk in boundaries:
            bvals.append(boundaries[pk])
    # unique boundaries in order
    buniq = []
    for b in bvals:
        if b not in buniq:
            buniq.append(b)

    # Assign boundaries per position:
    # - If we have at least 2 distinct boundaries and at least 2 proteins, we try to split positions evenly:
    #   first half -> boundary1, sSORCnd half -> boundary2.
    # This matches common formatting: (posA)-(posB) in the same field.
    # - Else use boundary1 for all.
    per_side: list[str] = []
    if not buniq:
        per_side = ["UNKNOWN"] * len(positions)
        site_class = "UNKNOWN"
        return ParsedSite(proteins=proteins, positions=positions, per_side=per_side, site_class=site_class)

    if len(buniq) >= 2 and len(positions) >= 2:
        b_for_pos = [buniq[0], buniq[1]] + [buniq[1]] * (len(positions) - 2)
    else:
        b_for_pos = [buniq[0]] * len(positions)

    for pos, b in zip(positions, b_for_pos):
        per_side.append("SORC" if pos > b else "BACSU")

    if all(x == "BACSU" for x in per_side):
        site_class = "BACSU-only"
    elif all(x == "SORC" for x in per_side):
        site_class = "SORC-only"
    elif any(x == "UNKNOWN" for x in per_side):
        site_class = "UNKNOWN"
    else:
        site_class = "MIXED"

    return ParsedSite(proteins=proteins, positions=positions, per_side=per_side, site_class=site_class)


def support_bin(n: object) -> str:
    if pd.isna(n):
        return "NA"
    try:
        v = int(n)
    except Exception:
        return "NA"
    if v == 1:
        return "1"
    if v == 2:
        return "2"
    if 3 <= v <= 5:
        return "3-5"
    if 6 <= v <= 10:
        return "6-10"
    return "11+"


def savefig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def main() -> None:
    args = parse_args()
    fasta_path = Path(args.fasta)
    sites_dir = Path(args.sites_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    boundaries = read_fasta_boundaries(fasta_path, style=args.style)

    sites_files = sorted(sites_dir.glob("*sites.csv"))
    if not sites_files:
        raise SystemExit(f"No *sites.csv files found in: {sites_dir}")

    all_rows: list[pd.DataFrame] = []

    for fp in sites_files:
        df = read_plink_sites_only(fp)
        if "Protein" not in df.columns:
            if "Proteins" in df.columns:
                df = df.rename(columns={"Proteins": "Protein"})
            else:
                continue

        df = df.copy()
        df["type_label"] = infer_type_label(fp.name)

        for c in ["Unique_Peptide_Number", "Spectrum_Number"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        parsed = df["Protein"].map(lambda s: classify_site(s, boundaries))
        df["site_class"] = parsed.map(lambda x: x.site_class)
        df["positions"] = parsed.map(lambda x: ";".join(str(p) for p in x.positions))
        df["sides"] = parsed.map(lambda x: ";".join(x.per_side))
        df["protein_keys"] = parsed.map(lambda x: ";".join(x.proteins))

        all_rows.append(df)

    sites_all = pd.concat(all_rows, ignore_index=True)

    annotated_path = outdir / "sites_annotated.tsv"
    sites_all.to_csv(annotated_path, sep="\t", index=False)

    summary = (
        sites_all.groupby(["type_label", "site_class"], observed=True)
        .size()
        .rename("n_sites")
        .reset_index()
    )
    summary_path = outdir / "summary_sites_by_type.csv"
    summary.to_csv(summary_path, index=False)

    pivot = summary.pivot(index="type_label", columns="site_class", values="n_sites").fillna(0).astype(int)

    for c in ["BACSU-only", "SORC-only", "MIXED", "UNKNOWN"]:
        if c not in pivot.columns:
            pivot[c] = 0
    pivot = pivot[["BACSU-only", "SORC-only", "MIXED", "UNKNOWN"]]
    
    want_types = ["cross-linked", "mono-linked", "loop-linked"]
    want_cols = ["BACSU-only", "SORC-only", "MIXED"]

    p2 = pivot.copy()
    p2 = p2.reindex([t for t in want_types if t in p2.index]).fillna(0).astype(int)

    for c in want_cols:
        if c not in p2.columns:
            p2[c] = 0
    p2 = p2[want_cols]

    x = np.arange(len(p2.index))
    w = 0.26

    plt.figure(figsize=(12.0, 5.8))
    b1 = plt.bar(x - w, p2["BACSU-only"].to_numpy(), width=w, label="BACSU-only")
    b2 = plt.bar(x,      p2["SORC-only"].to_numpy(), width=w, label="SORC-only")
    b3 = plt.bar(x + w,  p2["MIXED"].to_numpy(),     width=w, label="MIXED")

    plt.xticks(x, p2.index.astype(str), rotation=20, ha="right")
    plt.ylabel("Unique sites (count)")
    plt.title("Unique sites by type and class (BACSU vs SORC vs MIXED)")
    plt.legend()
    plt.box(False)

    if args.log_y:
        plt.yscale("log")
        ymax = max(1, int(p2.to_numpy().max()))
        plt.ylim(1, ymax * 3)
        # labels on log scale
        for bars in (b1, b2, b3):
            for rect in bars:
                h = rect.get_height()
                if h <= 0:
                    continue
                plt.text(rect.get_x() + rect.get_width()/2, h * 1.08, f"{int(h)}",
                        ha="center", va="bottom", fontsize=9, clip_on=False)
    else:
        ymax = max(1, int(p2.to_numpy().max()))
        plt.ylim(0, ymax * 1.25)
        for bars in (b1, b2, b3):
            for rect in bars:
                h = rect.get_height()
                if h <= 0:
                    continue
                plt.text(rect.get_x() + rect.get_width()/2, h + ymax * 0.02, f"{int(h)}",
                        ha="center", va="bottom", fontsize=9, clip_on=False)

    savefig(outdir / "sites_grouped_BACSU_SORC_MIXED.png")

    preferred = ["cross-linked", "mono-linked", "loop-linked", "linear", "unknown-type"]
    idx = [t for t in preferred if t in pivot.index] + [t for t in pivot.index if t not in preferred]
    pivot = pivot.reindex(idx)

    for c in ["BACSU-only", "SORC-only", "MIXED", "UNKNOWN"]:
        if c not in pivot.columns:
            pivot[c] = 0
    pivot = pivot[["BACSU-only", "SORC-only", "MIXED", "UNKNOWN"]]

    x = np.arange(len(pivot.index))
    bottom = np.zeros(len(pivot.index), dtype=int)

    plt.figure(figsize=(11.0, 6.0))
    for col in pivot.columns:
        vals = pivot[col].to_numpy()
        plt.bar(x, vals, bottom=bottom, label=col)
        bottom += vals

    plt.xticks(x, pivot.index.astype(str), rotation=20, ha="right")
    plt.ylabel("Unique sites (count)")
    plt.title("Site classes by type (stacked counts)")
    plt.legend()
    plt.box(False)
    if args.log_y:
        plt.yscale("log")
        ymax = max(1, int(pivot.to_numpy().max()))
        plt.ylim(1, ymax * 3)
    else:
        ymax = max(1, int(pivot.to_numpy().max()))
        plt.ylim(0, ymax * 1.25)

    totals = pivot.sum(axis=1).to_numpy()
    for i, t in enumerate(totals):
        if args.log_y:
            y = max(1, t) * 1.08
        else:
            y = t + ymax * 0.02
        plt.text(i, y, f"{int(t)}", ha="center", va="bottom", fontsize=9, clip_on=False)

    savefig(outdir / "sites_stacked_counts.png")

    frac = pivot.div(pivot.sum(axis=1), axis=0).replace([np.inf, -np.inf], np.nan).fillna(0)
    x = np.arange(len(frac.index))
    bottom = np.zeros(len(frac.index), dtype=float)

    plt.figure(figsize=(11.0, 6.0))
    for col in frac.columns:
        vals = frac[col].to_numpy()
        plt.bar(x, vals, bottom=bottom, label=col)
        bottom += vals

    plt.xticks(x, frac.index.astype(str), rotation=20, ha="right")
    plt.ylabel("Fraction of sites")
    plt.ylim(0, 1)
    plt.title("Site classes by type (stacked fractions)")
    plt.legend()
    plt.box(False)
    savefig(outdir / "sites_stacked_fractions.png")

    if "Spectrum_Number" in sites_all.columns:
        tmp = sites_all.copy()
        tmp["support_bin"] = tmp["Spectrum_Number"].map(support_bin)
        tmp["is_wrong"] = tmp["site_class"].isin(["SORC-only", "MIXED"]).astype(int)

        bins_order = ["1", "2", "3-5", "6-10", "11+", "NA"]

        g = (
            tmp.groupby(["type_label", "support_bin"], observed=True)["is_wrong"]
            .agg(["count", "mean"])
            .reset_index()
        )

        for t in pivot.index:
            gt = g[g["type_label"] == t].copy()
            if gt.empty:
                continue
            gt["support_bin"] = pd.Categorical(gt["support_bin"], categories=bins_order, ordered=True)
            gt = gt.sort_values("support_bin")

            plt.figure(figsize=(9.0, 5.0))
            bars = plt.bar(gt["support_bin"].astype(str), gt["mean"].to_numpy())
            plt.bar_label(bars, padding=3, fmt="%.3f")
            plt.ylabel("Wrong-site fraction (SORC-only or MIXED)")
            plt.ylim(0, min(0.5, max(0.05, float(np.nanmax(gt["mean"])) * 1.2 if np.isfinite(np.nanmax(gt["mean"])) else 0.1)))
            plt.title(f"Wrong-site fraction vs spectral support — {t}")
            plt.box(False)
            savefig(outdir / f"wrong_fraction_vs_support_{t}.png")

        g.to_csv(outdir / "wrong_fraction_vs_support_bins.csv", index=False)

    print(f"[OK] Wrote annotated sites: {annotated_path}")
    print(f"[OK] Wrote summary:        {summary_path}")
    print(f"[OK] Wrote plots to:       {outdir.resolve()}")


if __name__ == "__main__":
    main()