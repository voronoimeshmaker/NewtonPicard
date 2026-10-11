"""Skip a pending campaign once every configured method/family is recorded.
Nonconvergence is a scientific result; completeness is based on unique case IDs.
"""
from pathlib import Path
import csv
import sys
import yaml
root = Path(__file__).resolve().parent
for name in sys.argv[1:]:
    cfg = yaml.safe_load((root / name / "config.yaml").read_text())
    expected = {str(i) for i in range(1, cfg["sampling"]["sample_count"] + 1)}
    required = {(m, f) for m in cfg["methods"] for f in cfg["mesh"]["families"]}
    completed = set()
    for run in sorted((root / name / "Result").glob("run_*")):
        if not (run / "summary.csv").exists() or not (run / "cases.csv").exists():
            continue
        with (run / "cases.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        for method, family in required:
            ids = [q.get("case_id") for q in rows
                   if q.get("method") == method and q.get("mesh_family") == family]
            if len(ids) == len(expected) and set(ids) == expected:
                completed.add((method, family))
    if required - completed:
        print(name)
