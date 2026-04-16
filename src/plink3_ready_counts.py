#python plink3_ready_counts.py --ready-dir data/plink3/ready --outdir output/plink3_ready_counts
import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count pLink3 site-level identifications from numbered result folders "
            "inside a ready directory."
        )
    )
    parser.add_argument(
        "--ready-dir",
        required=True,
        help="Path to the pLink3 ready directory containing numbered folders such as 0, 1, ..., 17.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory where the summary tables will be written.",
    )
    return parser.parse_args()


def get_numbered_result_dirs(ready_dir: Path) -> list[Path]:
    return sorted(
        [path for path in ready_dir.iterdir() if path.is_dir() and path.name.isdigit()],
        key=lambda path: int(path.name),
    )


def find_single_report(report_dir: Path, suffix: str) -> Path:
    matches = sorted(report_dir.glob(f"result_*.{suffix}"))
    if not matches:
        raise FileNotFoundError(f"No report matching '*.{suffix}' found in {report_dir}")
    if len(matches) > 1:
        raise FileExistsError(
            f"Expected one report matching '*.{suffix}' in {report_dir}, found {len(matches)}"
        )
    return matches[0]


def read_site_rows(csv_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            return rows

        reader = csv.DictReader([header_line, *(
            line for line in handle if line.strip() and not line.startswith(",")
        )])
        for row in reader:
            rows.append({key: (value or "").strip() for key, value in row.items()})
    return rows


def summarize_result_dir(result_dir: Path) -> dict[str, int | str]:
    report_dir = result_dir / "reports"
    cross_rows = read_site_rows(find_single_report(report_dir, "filtered_cross-linked_sites.csv"))
    loop_rows = read_site_rows(find_single_report(report_dir, "filtered_loop-linked_sites.csv"))
    mono_rows = read_site_rows(find_single_report(report_dir, "filtered_mono-linked_sites.csv"))

    inter_count = sum(1 for row in cross_rows if row.get("Protein_Type") == "Inter-Protein")
    intra_count = sum(1 for row in cross_rows if row.get("Protein_Type") == "Intra-Protein")

    return {
        "file": result_dir.name,
        "cross_linked": len(cross_rows),
        "inter_protein_cross_linked": inter_count,
        "intra_protein_cross_linked": intra_count,
        "loop_linked": len(loop_rows),
        "mono_linked": len(mono_rows),
    }


def write_tsv(path: Path, rows: list[dict[str, int | str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    ready_dir = Path(args.ready_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    result_dirs = get_numbered_result_dirs(ready_dir)
    if not result_dirs:
        raise FileNotFoundError(f"No numbered result directories found in {ready_dir}")

    per_result_rows = [summarize_result_dir(result_dir) for result_dir in result_dirs]

    total_row = {
        "file": "ALL",
        "cross_linked": sum(int(row["cross_linked"]) for row in per_result_rows),
        "inter_protein_cross_linked": sum(
            int(row["inter_protein_cross_linked"]) for row in per_result_rows
        ),
        "intra_protein_cross_linked": sum(
            int(row["intra_protein_cross_linked"]) for row in per_result_rows
        ),
        "loop_linked": sum(int(row["loop_linked"]) for row in per_result_rows),
        "mono_linked": sum(int(row["mono_linked"]) for row in per_result_rows),
    }

    fieldnames = [
        "file",
        "cross_linked",
        "inter_protein_cross_linked",
        "intra_protein_cross_linked",
        "loop_linked",
        "mono_linked",
    ]

    write_tsv(outdir / "plink3_ready_counts_all.tsv", [total_row], fieldnames)
    write_tsv(outdir / "plink3_ready_counts_by_file.tsv", per_result_rows, fieldnames)

    print(f"Wrote {outdir / 'plink3_ready_counts_all.tsv'}")
    print(f"Wrote {outdir / 'plink3_ready_counts_by_file.tsv'}")


if __name__ == "__main__":
    main()
