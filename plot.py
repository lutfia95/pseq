from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plotter.overall_stats import OverallStats

# conda install -y python-dateutil

def _ensure_dir(path: str) -> str:
    path = os.path.abspath(path)
    os.makedirs(path, exist_ok=True)
    return path


def read_summary_tsv(path: str) -> Tuple[OverallStats, pd.DataFrame]:
    lines: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.rstrip("\n") for ln in f]

    overall: Optional[OverallStats] = None

    try:
        i_overall = next(i for i, ln in enumerate(lines) if ln.strip() == "# overall")
        header = lines[i_overall + 1].strip()
        values = lines[i_overall + 2].strip()
        cols = header.split("\t")
        vals = values.split("\t")
        d = {c: int(v) for c, v in zip(cols, vals)}
        overall = OverallStats(
            total_targets=d["total_targets"],
            matched=d["matched"],
            unmatched=d["unmatched"],
            total_target_len=d["total_target_len"],
            total_decoy_len=d["total_decoy_len"],
            total_target_k=d["total_target_K"],
            total_target_r=d["total_target_R"],
            total_decoy_k=d["total_decoy_K"],
            total_decoy_r=d["total_decoy_R"],
        )
    except StopIteration:
        raise RuntimeError("Could not find '# overall' section in summary.tsv")
    except Exception as e:
        raise RuntimeError(f"Failed parsing overall section: {e}")

    try:
        i_pt = next(i for i, ln in enumerate(lines) if ln.strip() == "# per_target")
        pt_header = lines[i_pt + 1].strip()
        pt_cols = pt_header.split("\t")
        pt_rows: List[List[str]] = []
        for ln in lines[i_pt + 2 :]:
            if not ln.strip():
                break
            if ln.startswith("#"):
                break
            pt_rows.append(ln.split("\t"))

        df = pd.DataFrame(pt_rows, columns=pt_cols)
    except StopIteration:
        raise RuntimeError("Could not find '# per_target' section in summary.tsv")
    except Exception as e:
        raise RuntimeError(f"Failed parsing per_target section: {e}")

    num_cols = [
        "target_len",
        "decoy_len",
        "len_diff",
        "target_K",
        "target_R",
        "target_KR",
        "decoy_K",
        "decoy_R",
        "decoy_KR",
        "kr_diff",
        "score",
    ]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["target_len", "decoy_len", "target_KR", "decoy_KR", "len_diff", "kr_diff", "score"])
    return overall, df


def _savefig(outdir: str, fname: str) -> None:
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, fname), dpi=200)
    plt.close()


def plot_all(summary_path: str, outdir: str, top_n_outliers: int = 30) -> None:
    overall, df = read_summary_tsv(summary_path)
    outdir = _ensure_dir(outdir)

    x_len = df["target_len"].to_numpy()
    y_len = df["decoy_len"].to_numpy()
    x_kr = df["target_KR"].to_numpy()
    y_kr = df["decoy_KR"].to_numpy()
    len_diff = df["len_diff"].to_numpy()
    kr_diff = df["kr_diff"].to_numpy()

    max_len = int(max(x_len.max(), y_len.max()))
    max_kr = int(max(x_kr.max(), y_kr.max()))

    plt.figure()
    plt.scatter(x_len, y_len, s=8, alpha=0.35)
    plt.plot([0, max_len], [0, max_len])
    plt.xlabel("target_len")
    plt.ylabel("decoy_len")
    plt.title("Target length vs Decoy length")
    _savefig(outdir, "01_scatter_target_vs_decoy_len.png")

    plt.figure()
    bins = min(200, max(10, int(np.sqrt(len(len_diff)))))
    plt.hist(len_diff, bins=bins)
    plt.xlabel("len_diff (aa)")
    plt.ylabel("count")
    plt.title("Length difference distribution")
    _savefig(outdir, "02_hist_len_diff.png")

    plt.figure()
    s = np.sort(len_diff)
    y = np.arange(1, len(s) + 1) / len(s)
    plt.step(s, y, where="post")
    plt.xlabel("len_diff (aa)")
    plt.ylabel("CDF")
    plt.title("CDF of length difference")
    _savefig(outdir, "03_cdf_len_diff.png")

    plt.figure()
    plt.scatter(x_kr, y_kr, s=8, alpha=0.35)
    plt.plot([0, max_kr], [0, max_kr])
    plt.xlabel("target_KR")
    plt.ylabel("decoy_KR")
    plt.title("Target (K+R) vs Decoy (K+R)")
    _savefig(outdir, "04_scatter_target_vs_decoy_KR.png")

    plt.figure()
    bins = min(200, max(10, int(np.sqrt(len(kr_diff)))))
    plt.hist(kr_diff, bins=bins)
    plt.xlabel("kr_diff (count)")
    plt.ylabel("count")
    plt.title("(K+R) difference distribution")
    _savefig(outdir, "05_hist_kr_diff.png")

    plt.figure()
    s = np.sort(kr_diff)
    y = np.arange(1, len(s) + 1) / len(s)
    plt.step(s, y, where="post")
    plt.xlabel("kr_diff (count)")
    plt.ylabel("CDF")
    plt.title("CDF of (K+R) difference")
    _savefig(outdir, "06_cdf_kr_diff.png")

    plt.figure()
    plt.hexbin(x_len, kr_diff, gridsize=60, mincnt=1)
    plt.xlabel("target_len")
    plt.ylabel("kr_diff (count)")
    plt.title("KR mismatch vs target length (hexbin density)")
    cb = plt.colorbar()
    cb.set_label("count")
    _savefig(outdir, "07_hexbin_target_len_vs_kr_diff.png")

    plt.figure()
    labels = ["K", "R"]
    target_vals = [overall.total_target_k, overall.total_target_r]
    decoy_vals = [overall.total_decoy_k, overall.total_decoy_r]
    x = np.arange(len(labels))
    w = 0.38
    plt.bar(x - w / 2, target_vals, width=w, label="targets")
    plt.bar(x + w / 2, decoy_vals, width=w, label="decoys")
    plt.xticks(x, labels)
    plt.ylabel("total count")
    plt.title("Overall composition totals (K and R)")
    plt.legend()
    _savefig(outdir, "08_bar_overall_K_R_totals.png")

    plt.figure()
    s = df["score"].to_numpy()
    bins = min(200, max(10, int(np.sqrt(len(s)))))
    plt.hist(s, bins=bins)
    plt.xlabel("score")
    plt.ylabel("count")
    plt.title("Score distribution")
    _savefig(outdir, "09_hist_score.png")

    outliers = df.sort_values(["kr_diff", "len_diff", "score"], ascending=[False, False, False]).head(top_n_outliers)
    plt.figure(figsize=(10, 6))
    idx = np.arange(len(outliers))
    plt.bar(idx, outliers["score"].to_numpy())
    plt.xticks(idx, outliers["target_name"].tolist(), rotation=90, fontsize=7)
    plt.ylabel("score")
    plt.title(f"Top {len(outliers)} worst matches by (kr_diff, len_diff, score)")
    _savefig(outdir, "10_bar_top_outliers_by_score.png")

    outliers_path = os.path.join(outdir, "top_outliers.tsv")
    outliers.to_csv(outliers_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, help="Path to summary.tsv produced by matcher")
    ap.add_argument("--outdir", required=True, help="Output directory for plots")
    ap.add_argument("--top-n-outliers", type=int, default=30, help="How many worst matches to export/plot")
    args = ap.parse_args()
    plot_all(args.summary, args.outdir, top_n_outliers=args.top_n_outliers)


if __name__ == "__main__":
    main()