import csv
import re
import sys

KEY_RE = re.compile(r"\b([A-Z]{2})=")  # OS= OX= GN= AC= SS= PC=? more? 

def parse_fasta(fasta_path):
    db = {}
    cur_id = None
    cur_hdr = None
    seq_chunks = []

    def commit():
        nonlocal cur_id, cur_hdr, seq_chunks
        if not cur_id:
            return
        protein_name, meta = parse_header(cur_hdr)
        db[cur_id] = {
            "protein_name": protein_name,
            "OS": meta.get("OS", ""),
            "OX": meta.get("OX", ""),
            "GN": meta.get("GN", ""),
            "AC": meta.get("AC", ""),
            "SS": meta.get("SS", ""),
            "PC": meta.get("PC", ""),
            "seq": "".join(seq_chunks),
        }

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                commit()
                seq_chunks = []
                cur_hdr = line[1:].strip()
                cur_id = cur_hdr.split(None, 1)[0] if cur_hdr else None
            else:
                if cur_id is not None:
                    seq_chunks.append(line.strip())

    commit()
    return db

def parse_header(hdr):
    if not hdr:
        return "", {}

    parts = hdr.split(None, 1)
    rest = parts[1] if len(parts) > 1 else ""

    m = KEY_RE.search(rest)
    if m:
        protein_name = rest[:m.start()].strip()
        meta_str = rest[m.start():].strip()
    else:
        protein_name = rest.strip()
        meta_str = ""

    meta = {}
    if meta_str:
        keys = list(KEY_RE.finditer(meta_str))
        for i, km in enumerate(keys):
            k = km.group(1)
            v_start = km.end()
            v_end = keys[i + 1].start() if i + 1 < len(keys) else len(meta_str)
            meta[k] = meta_str[v_start:v_end].strip()

    return protein_name, meta

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.csv> <input.fasta> <output.csv>", file=sys.stderr)
        sys.exit(1)

    csv_in, fasta_in, csv_out = sys.argv[1], sys.argv[2], sys.argv[3]

    fasta_db = parse_fasta(fasta_in)

    # new_cols = ["OS", "OX", "GN", "AC", "SS", "PC", "protein_name", "seq"]
    new_cols = ["OS", "OX", "GN", "AC", "SS", "PC", "protein_name", "seq", "matched_id"] # if more than on id found! 
    def first_hit_id(protein_group_value, fasta_db):
        ids = [x.strip() for x in (protein_group_value or "").split(";") if x.strip()]
        for pid in ids:
            if pid in fasta_db:
                return pid
        return ""
    not_found = []

    with open(csv_in, "r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin)
        if not reader.fieldnames:
            print("Empty/invalid CSV.", file=sys.stderr)
            sys.exit(1)

        out_fieldnames = list(reader.fieldnames) + [c for c in new_cols if c not in reader.fieldnames]

        with open(csv_out, "w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fieldnames, quoting=csv.QUOTE_MINIMAL)
            writer.writeheader()

            for row in reader:
                pg_raw = (row.get("Protein.Group") or "").strip()
                hit_id = first_hit_id(pg_raw, fasta_db)
                hit = fasta_db.get(hit_id) if hit_id else None
                if hit is None:
                    not_found.append(pg_raw)
                    for c in new_cols:
                        row[c] = ""
                else:
                    row["matched_id"] = hit_id
                    row["OS"] = hit["OS"]
                    row["OX"] = hit["OX"]
                    row["GN"] = hit["GN"]
                    row["AC"] = hit["AC"]
                    row["SS"] = hit["SS"]
                    row["PC"] = hit["PC"]
                    row["protein_name"] = hit["protein_name"]
                    row["seq"] = hit["seq"]
                writer.writerow(row)

    for x in not_found:
        if x:
            print("not found: ", x)
    print(len([x for x in not_found if x]))

if __name__ == "__main__":
    main()