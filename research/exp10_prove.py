"""Experiment 10: computer-assisted existence proofs.

Closure to 1e-12 in floating point is numerical evidence, not proof. Each orbit
is refined to about 50 digits and then validated with rigorous ball arithmetic
(proof.py): a unique exact solution is shown to lie within rho of the initial
conditions we publish. Results go to results/proofs.json.
"""
import glob
import json
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np

import family
import proof
import tb

RES = os.path.join(tb.HERE, "results")
SETTINGS = ((384, 45, 1e-70), (512, 60, 1e-95))  # (bits, series order, step tolerance); second is the retry


def job(args):
    label, z, turn, note = args
    t0 = time.time()
    out = {"label": label, "turn": f"{abs(turn[0])}/{turn[1]}", "note": note, "proved": False}
    for bits, order, tol in SETTINGS:
        try:
            proof.set_precision(bits)
            zz, size = proof.newton(list(z), turn, order, tol)
            report = proof.certify(zz, turn, order, tol)
        except (ArithmeticError, ZeroDivisionError, ValueError) as err:
            out["error"] = f"{type(err).__name__}: {err}"
            continue
        out.update(report)
        out["newton_residual"] = float(size)
        if report["proved"]:
            out.pop("error", None)
            break
    out["half_lap"] = float(z[4])
    out["closes_after_time"] = f"at most {2 * turn[1]} half-laps"
    out["seconds"] = round(time.time() - t0)
    print(label, "PROVED" if out["proved"] else "not proved", out.get("rho"), out.get("contraction"),
          out.get("error", ""), f"{out['seconds']}s", flush=True)
    return out


def family_jobs(fractions):
    rows = np.load(os.path.join(RES, "family_trace.npy"))
    checked = json.load(open(os.path.join(RES, "exp9_catalogue.json")))
    jobs = []
    for frac in fractions:
        f = Fraction(frac)
        z, _ = family.member(rows, -np.pi * f.numerator / f.denominator)
        published = checked["family_member_published"].get(str(f))
        jobs.append((f"family turn {f}", z, (-f.numerator, f.denominator),
                     "in the published catalogue" if published else "not in the published catalogue"))
    return jobs


def absent_jobs():
    """Inertially periodic orbits from the search that are absent from the catalogue,
    cheapest first (long periods and close approaches need many more validated steps).
    """
    checked = json.load(open(os.path.join(RES, "exp9_catalogue.json")))
    found = []
    for path in sorted(glob.glob(os.path.join(RES, "spatial_orbits_seed*.json"))):
        found += json.load(open(path))
    absent = {round(u["period"], 6) for u in checked["unmatched"] if not u["built"]}
    picks = sorted((o for o in found if round(o["period"], 6) in absent), key=lambda o: o["period"] / o["min_sep"] ** 0.5)
    return [(f"spatial orbit, period {o['period']:.4f}, spin {o['scaled_L']:.3f}",
             np.array([*o["params"], o["period"] / 2]), (0, 1), "not in the published catalogue") for o in picks]


if __name__ == "__main__":
    # usage: exp10_prove.py family 1/3 1/4 ...   |   exp10_prove.py absent <how many> <time budget in seconds>
    mode = sys.argv[1] if len(sys.argv) > 1 else "family"
    if mode == "family":
        jobs, budget = family_jobs(sys.argv[2:] or ["1/3", "1/4", "1/5", "2/5", "3/7", "4/9", "2/7", "1/6"]), 1e9
    else:
        jobs, budget = absent_jobs()[:int(sys.argv[2])], float(sys.argv[3])
    out_file = os.path.join(RES, "proofs.json")
    results = {r["label"]: r for r in (json.load(open(out_file)) if os.path.exists(out_file) else [])}
    jobs = [j for j in jobs if not results.get(j[0], {}).get("proved")]
    print(len(jobs), "proofs to attempt", flush=True)
    t0 = time.time()
    pool = Pool(8)
    pending = [pool.apply_async(job, (j,)) for j in jobs]
    while pending and time.time() - t0 < budget:
        for p in [p for p in pending if p.ready()]:
            r = p.get()
            results[r["label"]] = r
            pending.remove(p)
            json.dump(list(results.values()), open(out_file, "w"), indent=1)
        time.sleep(2)
    pool.terminate()
    print(sum(r["proved"] for r in results.values()), "proved in total;", len(pending), "not finished in the time budget")
