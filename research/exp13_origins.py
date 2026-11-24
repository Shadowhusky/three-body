"""Experiment 13: where do the 3D families come from?

Families from experiment 12 that could be followed down to zero angular
momentum are re-traced, and the end of each is compared with our list of planar
zero-angular-momentum periodic orbits (experiment 3) and with the collinear
Schubart orbit. At zero angular momentum the isosceles parameters (y, a, b, c)
reduce to the planar slice with (p1, p2) = (a, c), so an end can be identified
both by its initial condition and by its scale-invariant period.
"""
import json
import os
from multiprocessing import Pool

import numpy as np

import family
import tb

RES = os.path.join(tb.HERE, "results")
SCHUBART = 4.69808380 * (1 / 6) ** 1.5 * 3 ** 2.5  # Henon's period (m = 1/3, E = -1/6) in our scale-invariant units
TOL = 3e-6  # relative period agreement required; the ends sit at spin ~1e-3, where periods are within ~1e-7 of their limits
NAMES = {9.2377: "Figure-eight", 19.9325: "Butterfly I", 19.9313: "Butterfly II", 24.2056: "Moth I", 29.6021: "Yin-yang I",
         34.8022: "Butterfly III", 36.7714: "Dragonfly", 39.6436: "Goggles", 42.7821: "Moth II", 53.8569: "Moth III"}


def retrace(f):
    seed = f["seed"]
    rows = family.trace_arclength(np.r_[seed["q"], seed["t"]], seed["phi"], max_points=600, budget=120)
    ends = []
    for row in (rows[0], rows[-1]):
        L, tau, _ = family.invariants(row[:5])
        if L < 3e-3:
            ends.append({"params": [float(v) for v in row[:4]], "spin": float(L), "period": float(tau), "phi_over_pi": float(row[5] / np.pi)})
    return f["id"], ends


def identify(end, planar):
    """Planar orbit (or Schubart multiple) that this zero-spin end coincides with."""
    y, a, b, c = end["params"]
    for k in (1, 2, 3, 4):
        if abs(end["period"] - k * SCHUBART) < TOL * end["period"] and abs(c) < 0.05:
            return {"orbit": "Schubart" if k == 1 else f"Schubart x{k}", "period_gap": abs(end["period"] / (k * SCHUBART) - 1), "by": "period, collinear start"}
    best = None
    for o in planar:
        for k in (1, 2, 3, 4):  # the family's relative period may be a whole fraction of the planar period
            gap = abs(end["period"] * k / o["scaled_period"] - 1)
            if gap < TOL and (best is None or gap < best["period_gap"]):
                same_start = abs(abs(a) - o["params"][0]) < 0.01 and abs(abs(c) - o["params"][1]) < 0.01 and abs(y) < 0.01 and abs(b) < 0.01
                name = next((n for T, n in NAMES.items() if abs(T - o["scaled_period"]) < 2e-4), f"planar orbit with T* = {o['scaled_period']:.4f}")
                best = {"orbit": name, "fraction_of_planar_period": f"1/{k}" if k > 1 else "1", "period_gap": gap,
                        "by": "period and initial condition" if same_start else "period"}
    return best


if __name__ == "__main__":
    fams = json.load(open(os.path.join(RES, "exp12_families.json")))["families"]
    planar = json.load(open(os.path.join(RES, "planar_orbits.json")))
    todo = [f for f in fams if "collinear limit" in f["ends"]]
    with Pool(8) as pool:
        traced = dict(pool.map(retrace, todo, chunksize=1))
    report = []
    for f in todo:
        for end in traced[f["id"]]:
            found = identify(end, planar)
            report.append({"family": f["id"], "members": len(f["members_ours"]) + len(f["members_catalogue"]),
                           "catalogue_members": f["members_catalogue"], **end, "origin": found})
            print(f"family {f['id']:2d} ({report[-1]['members']} members {f['members_catalogue']}): spin {end['spin']:.4f}, period {end['period']:.4f}"
                  f" -> {found['orbit'] + ' (' + found['by'] + ', gap %.1e)' % found['period_gap'] if found else 'no planar orbit matched'}")
    matched = [r for r in report if r["origin"]]
    json.dump({"schubart_period": SCHUBART, "ends": report, "ends_total": len(report), "ends_identified": len(matched),
               "families": len({r["family"] for r in report}), "families_identified": len({r["family"] for r in matched})},
              open(os.path.join(RES, "exp13_origins.json"), "w"), indent=1)
    print(f"{len(matched)} of {len(report)} zero-spin ends identified, in {len({r['family'] for r in matched})} of {len({r['family'] for r in report})} families")
