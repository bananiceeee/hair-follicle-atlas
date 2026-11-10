import csv
import json
import sys
from pathlib import Path

SETS = [
    ("HighCHD_DEGup_top50_symbol", "High-CHD-up DEG"),
    ("LowCHD_DEGdown_top50_symbol", "Low-CHD-up DEG"),
    ("CHDpos_WGCNAdarkgreen_top50_by_absMM_symbol", "CHD-positive WGCNA"),
    ("CHDneg_WGCNAmagenta_top50_by_absMM_symbol", "CHD-negative WGCNA"),
    ("ATP6V1B1_All_top50", "ATP6V1B1 metagene"),
    ("GJA1_All_top50", "GJA1 metagene"),
    ("GPRC5D_All_top50", "GPRC5D metagene"),
]


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


source = Path(sys.argv[1])
input_rows = read_csv(source / "input_members_top50.csv")
mapped_rows = read_csv(source / "mapped_members_top50.csv")
assert {row["program_id"] for row in input_rows} == {key for key, _ in SETS}
assert {row["program"] for row in mapped_rows} == {label for _, label in SETS}
assert all(row["top_n"] == "50" for row in mapped_rows)
result = []
for key, label in SETS:
    rows = [row for row in input_rows if row["program_id"] == key]
    assert len(rows) == len({row["query_gene"].upper() for row in rows}) == 50
    mapped = {row["gene"] for row in mapped_rows if row["program"] == label}
    assert {row["matched_gene"] for row in rows if row["match_type"] != "missing"} == mapped
    result.append({
        "label": label,
        "members": {row["query_gene"].upper(): row["match_type"] != "missing" for row in rows},
    })

target = Path(__file__).resolve().parents[1] / "docs" / "data" / "goat-top50.json"
target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print("Top50 input / goat-skin matched:", [(entry["label"], sum(entry["members"].values())) for entry in result])
