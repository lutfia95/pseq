import re
import sys

# > header_seq1 && header_seq2 && idx_appended_seq_in_file && AA_IDX_seq1 (where the append occurs)
def parse_fasta(path):
    records = []
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records

def tag_from_header(header):
    m = re.search(r"\|([A-Z0-9]+)_([A-Z0-9]+)\b", header)
    if m:
        return m.group(2)
    m = re.search(r"\b([A-Z0-9]+)_([A-Z0-9]+)\b", header)
    if m:
        return m.group(2)
    return ""

def is_bsub(tag):
    return tag == "BACSU"

def is_ecoli(tag):
    return tag.startswith("ECO")

def is_yeast(tag):
    return tag == "YEAST"

def is_valid_append_target(tag):
    return is_ecoli(tag) or is_yeast(tag)

def merge_pairs(records):
    if len(records) % 2 != 0:
        raise SystemExit(f"FASTA has odd number of records ({len(records)}). Need pairs.")

    out = []
    for i in range(0, len(records), 2):
        h1, s1 = records[i]
        h2, s2 = records[i + 1]
        t1 = tag_from_header(h1)
        t2 = tag_from_header(h2)

        if not is_bsub(t1):
            raise SystemExit(f"Record {i+1} is not BACSU: {h1}")
        if not is_valid_append_target(t2):
            raise SystemExit(f"Record {i+2} is not ECO* or YEAST: {h2}")

        eco_index = i + 2
        merge_idx = len(s1)
        new_header = f">{h1} && {h2} && {eco_index} && {merge_idx}"
        out.append((new_header, s1 + s2))
    return out

def wrap_seq(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python append_seq.py input.fasta [output.fasta]")
    inp = sys.argv[1]
    outp = sys.argv[2] if len(sys.argv) > 2 else None

    records = parse_fasta(inp)
    merged = merge_pairs(records)

    text = []
    for h, s in merged:
        text.append(h)
        text.append(wrap_seq(s))
    text = "\n".join(text) + "\n"

    if outp:
        with open(outp, "w", encoding="utf-8") as f:
            f.write(text)
    else:
        sys.stdout.write(text)

if __name__ == "__main__":
    main()
