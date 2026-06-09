"""Build clean sample_metadata.csv from the raw metadata.csv."""
import re, csv, os, ast, sys

DATA_DIR = r"D:\PRJNA301507\data"
raw = os.path.join(DATA_DIR, "metadata.csv")
out = os.path.join(DATA_DIR, "sample_metadata.csv")

rows = []
with open(raw, encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    header = next(reader)
    for line in reader:
        if len(line) < len(header):
            continue
        d = dict(zip(header, line))
        cell_source_raw = d.get("CellSource", "")
        cell_type_raw  = d.get("CellType", "")
        chain_raw      = d.get("Chain", "")
        comment        = d.get("Comment", "")
        run_id         = d.get("RunId", "").strip()

        # parse python-list-like strings safely
        def parse_list(s):
            s = s.strip()
            if s.startswith("[") and s.endswith("]"):
                try:
                    return [str(x).strip() for x in ast.literal_eval(s)]
                except Exception:
                    return [s.strip("[]").strip("'\"")]
            return [s.strip("[]").strip("'\"")]

        cell_sources = parse_list(cell_source_raw)
        cell_types   = parse_list(cell_type_raw)
        chains       = parse_list(chain_raw)

        if "TRA" in chains:
            continue  # keep only TRB

        # determine compartment
        compartment = "Other"
        for src in cell_sources:
            sl = src.lower()
            if "blood" in sl:
                compartment = "Blood"
            elif "tumor" in sl:
                compartment = "Tumor"
            elif "lymph" in sl:
                compartment = "Lymph node"

        # determine subset
        subset = "Other"
        ct_lower = " ".join(cell_types).lower()
        if "treg" in ct_lower or ("cd25+" in ct_lower and "cd127-" in ct_lower):
            subset = "Treg"
        elif "memory" in ct_lower or ("cd45ro+" in ct_lower and "cd25-" in ct_lower):
            subset = "CD4_Memory"
        elif "cd45ra+" in ct_lower:
            subset = "Other"

        # extract patient ID from Comment
        m_pat = re.search(r"Patient\s+(\d+)", comment, re.IGNORECASE)
        patient_id = f"Patient{m_pat.group(1)}" if m_pat else ""

        # check the file exists
        fname = f"{run_id}.csv"
        if not os.path.exists(os.path.join(DATA_DIR, fname)):
            continue

        rows.append({
            "file_name": fname,
            "sample_id": run_id,
            "patient_id": patient_id,
            "subset": subset,
            "compartment": compartment,
            "comment": comment,
        })

# deduplicate by sample_id (keep first)
seen = set()
dedup = []
for r in rows:
    if r["sample_id"] not in seen:
        seen.add(r["sample_id"])
        dedup.append(r)

with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["file_name","sample_id","patient_id","subset","compartment","comment"])
    w.writeheader()
    w.writerows(dedup)

print(f"Wrote {len(dedup)} samples to {out}")
subsets = {}
for r in dedup:
    subsets.setdefault(r["subset"], set()).add(r["patient_id"])
for s in sorted(subsets):
    print(f"  {s}: {len(subsets[s])} patients, {sum(1 for r in dedup if r['subset']==s)} samples")
