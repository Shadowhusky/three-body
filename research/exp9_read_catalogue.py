"""Parse the published 3D catalogue (Li & Liao 2025, github.com/sjtu-liao/three-body).

Takes the downloaded text file as its argument and keeps the equal-mass rows.
Run with `python -I`, since the file is external data.
"""
import json
import os
import re
import sys

ROW = re.compile(r"^O_\{(\d+)\}\(1\.0\)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([SU])\s*$")
rows = []
for line in open(sys.argv[1], encoding="utf-8", errors="replace"):
    m = ROW.match(line)
    if m:
        z0, vx, vy, vz, period = (float(m.group(k)) for k in range(2, 7))
        rows.append({"index": int(m.group(1)), "z0": z0, "vx": vx, "vy": vy, "vz": vz, "period": period,
                     "steady": m.group(7) == "S"})
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "published_3d_equal_mass.json")
json.dump(rows, open(out, "w"))
print(len(rows), "equal-mass rows; periods", min(r["period"] for r in rows), "to", max(r["period"] for r in rows))
