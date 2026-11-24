"""Experiment 12: do the catalogued 3D periodic orbits organise into families?

Every periodic orbit of this symmetric class lies on a one-parameter family of
relative periodic orbits (family.py). Members with a rational rotation angle
are the ones that close in the inertial frame, and those are what a search for
periodic orbits finds. Here we take all our 3D orbits and all equal-mass
orbits of the published catalogue (Li & Liao 2025), continue families from
them, and record which orbits lie on which family.

An orbit is assigned to a family only after solving for the family member at
the orbit's own rotation angle and finding the same scale-invariant angular
momentum and relative period to 1e-6.
"""
import glob
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

import family
import tb

RES = os.path.join(tb.HERE, "results")
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 1500.0
T_CAP = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0  # longest relative half-period analysed
BATCH = 16


def load_orbits():
    orbits = []
    seen = set()
    ours = []
    for path in sorted(glob.glob(os.path.join(RES, "spatial_orbits_seed*.json"))):
        ours += json.load(open(path))
    for o in ours:
        key = (round(o["scaled_period"], 4), round(o["scaled_L"], 4))
        if key not in seen:
            seen.add(key)
            orbits.append({"source": "ours", "label": f"T={o['period']:.4f}", "q": o["params"], "period": o["period"],
                           "stable": o["growth_per_period"] < 1.02})
    published = os.path.join(RES, "published_3d_equal_mass.json")
    if os.path.exists(published):
        for r in json.load(open(published)):
            orbits.append({"source": "catalogue", "label": f"O{r['index']}", "q": [r["z0"] / 3, -r["vx"], r["vz"], r["vy"]],
                           "period": r["period"], "stable": r["steady"]})
    return orbits


def signature(o):
    """(scale-invariant L, relative period, |phi|) of an orbit, with its symmetric epoch."""
    try:
        t, phi, miss = family.symmetric_epoch(o["q"], o["period"])
        L, tau, _ = family.invariants(np.r_[o["q"], t])
    except (RuntimeError, np.linalg.LinAlgError):
        return None
    if miss > 1e-6:
        return None
    return {"t": t, "phi": phi, "L": float(L), "tau": float(tau)}


def trace(args):
    q, t, phi = args
    try:
        rows = family.trace_arclength(np.r_[q, t], phi, max_points=400, budget=90)
    except (RuntimeError, np.linalg.LinAlgError):
        return None
    inv = np.array([family.invariants(r[:5])[:2] for r in rows])
    return rows, inv


def on_curve(sig, rows, inv):
    """Does the orbit with this signature lie on the traced family? Solve the family
    at the orbit's own rotation angle from the nearest traced member and compare invariants."""
    # coarse test against the curve, densified so sparse stretches do not hide true members
    fine = np.linspace(0, len(rows) - 1, 20 * len(rows))
    L, tau, phi_c = (np.interp(fine, np.arange(len(rows)), col) for col in (inv[:, 0], inv[:, 1], rows[:, 5]))
    folded = np.abs((phi_c + np.pi / 2) % np.pi - np.pi / 2)
    d = np.hypot(np.hypot((L - sig["L"]) / max(sig["L"], 0.05), (tau - sig["tau"]) / sig["tau"]),
                 (folded - abs(sig["phi"])) / np.pi)
    if d.min() > 0.02:
        return False
    k = int(round(fine[int(np.argmin(d))]))
    targets = [s * sig["phi"] + n * np.pi for s in (1, -1) for n in range(-3, 4)]
    phi = min(targets, key=lambda p: abs(p - rows[k, 5]))
    try:
        z, miss = family.solve(rows[k, :5], phi, iters=30)
        L, tau, _ = family.invariants(z)
    except (RuntimeError, np.linalg.LinAlgError):
        return False
    return miss < 1e-10 and abs(L - sig["L"]) < 1e-6 and abs(tau - sig["tau"]) < 1e-6 * sig["tau"] + 1e-6


def ending(row, inv_row):
    if inv_row[0] < 3e-3:
        return "collinear limit"
    if abs(row[3]) < 2e-3:
        return "planar orbit"
    return "not followed further"


if __name__ == "__main__":
    t0 = time.time()
    orbits = load_orbits()
    cache = os.path.join(RES, "exp12_signatures.json")
    if os.path.exists(cache):
        sigs = json.load(open(cache))
    else:
        with Pool(tb.NCPU) as pool:
            sigs = pool.map(signature, orbits, chunksize=8)
        json.dump(sigs, open(cache, "w"))
    for o, s in zip(orbits, sigs):
        o["sig"] = s
    usable = [o for o in orbits if o["sig"] and o["sig"]["t"] <= T_CAP]
    print(f"{len(orbits)} orbits ({sum(o['source'] == 'ours' for o in orbits)} ours, "
          f"{sum(o['source'] == 'catalogue' for o in orbits)} catalogue); {len(usable)} have a relative half-period "
          f"up to {T_CAP} and are analysed", flush=True)

    remaining = sorted(usable, key=lambda o: o["sig"]["t"])
    families = []
    with Pool(tb.NCPU) as pool:
        while remaining and time.time() - t0 < BUDGET:
            seeds = remaining[:BATCH]
            traced = pool.map(trace, [(o["q"], o["sig"]["t"], o["sig"]["phi"]) for o in seeds], chunksize=1)
            for seed, result in zip(seeds, traced):
                if seed not in remaining:
                    continue  # an earlier family in this batch already claimed it
                remaining.remove(seed)
                if result is None or len(result[0]) < 5:
                    seed["family"] = None
                    continue
                rows, inv = result
                members = [seed] + [o for o in remaining if on_curve(o["sig"], rows, inv)]
                for o in members:
                    o["family"] = len(families)
                    if o in remaining:
                        remaining.remove(o)
                families.append({
                    "id": len(families), "points": len(rows),
                    "spin_range": [float(inv[:, 0].min()), float(inv[:, 0].max())],
                    "period_range": [float(inv[:, 1].min()), float(inv[:, 1].max())],
                    "phi_over_pi_range": [float(rows[:, 5].min() / np.pi), float(rows[:, 5].max() / np.pi)],
                    "ends": [ending(rows[0], inv[0]), ending(rows[-1], inv[-1])],
                    "members_ours": [o["label"] for o in members if o["source"] == "ours"],
                    "members_catalogue": [o["label"] for o in members if o["source"] == "catalogue"],
                    "stable_members": sum(o["stable"] for o in members),
                    "curve": [[float(inv[k, 0]), float(inv[k, 1]), float(rows[k, 5])] for k in range(len(rows))],
                    "seed": {"q": [float(v) for v in seed["q"]], "t": seed["sig"]["t"], "phi": seed["sig"]["phi"]},
                })
            print(f"{time.time() - t0:.0f}s: {len(families)} families, {len(remaining)} orbits unassigned", flush=True)

    assigned = [o for o in usable if o.get("family") is not None]
    sizes = sorted((len(f["members_ours"]) + len(f["members_catalogue"]) for f in families), reverse=True)
    summary = {
        "orbits": len(orbits), "relative_half_period_cap": T_CAP, "analysed": len(usable),
        "analysed_ours": sum(o["source"] == "ours" for o in usable),
        "analysed_catalogue": sum(o["source"] == "catalogue" for o in usable), "assigned": len(assigned),
        "unassigned": len(remaining), "families": len(families),
        "families_with_several_members": sum(s > 1 for s in sizes), "largest": sizes[:10],
        "catalogue_assigned": sum(o["source"] == "catalogue" for o in assigned),
        "catalogue_in_shared_families": sum(len(f["members_catalogue"]) for f in families
                                             if len(f["members_ours"]) + len(f["members_catalogue"]) > 1),
        "ours_assigned": sum(o["source"] == "ours" for o in assigned),
    }
    json.dump({"summary": summary, "families": families}, open(os.path.join(RES, "exp12_families.json"), "w"))
    print(json.dumps(summary, indent=1))
