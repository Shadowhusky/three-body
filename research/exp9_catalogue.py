"""Experiment 9: comparison with the published 3D catalogue (Li & Liao 2025).

Which of our 3D orbits are already in the catalogue, and which catalogue
entries lie on the family from experiment 8?

Orbits are matched by two scale-invariant quantities, T |E|^(3/2) and
L |E|^(1/2), which do not depend on either project's choice of units or of
initial configuration.
"""
import glob
import json
import os
from fractions import Fraction

import numpy as np

import family
import tb

RES = os.path.join(tb.HERE, "results")
published = json.load(open(os.path.join(RES, "published_3d_equal_mass.json")))


def their_state(r):
    return np.array([-1, 0, 0, 1, 0, 0, 0, 0, r["z0"], r["vx"], r["vy"], r["vz"], r["vx"], r["vy"], -r["vz"],
                     -2 * r["vx"], -2 * r["vy"], 0], dtype=float)


for r in published:
    s = their_state(r)
    E = tb.energy(s)[0]
    r["scaled_period"] = r["period"] * abs(E) ** 1.5
    r["scaled_L"] = float(np.linalg.norm(tb.angular_momentum(s)[0]) * abs(E) ** 0.5)
pub = np.array([[r["scaled_period"], r["scaled_L"]] for r in published])


def in_catalogue(o):
    """Index of the matching published orbit (possibly traversed a whole number of times), or None."""
    close = np.abs(pub[:, 1] - o["scaled_L"]) < 2e-5
    ratio = o["scaled_period"] / pub[:, 0]
    whole = np.abs(ratio - np.round(ratio)) < 1e-5 * np.maximum(1, ratio)
    inverse = np.abs(1 / ratio - np.round(1 / ratio)) < 1e-5 * np.maximum(1, 1 / ratio)
    hit = np.nonzero(close & (whole | inverse) & (np.round(ratio) + np.round(1 / ratio) >= 1))[0]
    return int(published[hit[0]]["index"]) if len(hit) else None


ours = []
for f in sorted(glob.glob(os.path.join(RES, "spatial_orbits_seed*.json"))):
    ours += json.load(open(f))
ours += [dict(o, built=True) for o in json.load(open(os.path.join(RES, "family_orbits.json"))) if not o["already_in_search"]]
seen, unique = set(), []
for o in ours:
    key = (round(o["scaled_period"], 4), round(o["scaled_L"], 4))
    if key not in seen:
        seen.add(key)
        unique.append(o)
ours = unique
for o in ours:
    o["published_index"] = in_catalogue(o)

matched = [o for o in ours if o["published_index"]]
unmatched = [o for o in ours if not o["published_index"]]
t_max = max(o["period"] for o in ours if not o.get("built"))
pub_in_range = [r for r in published if r["period"] <= t_max]
found_idx = {o["published_index"] for o in matched}
print(f"ours: {len(ours)} spatial orbits; {len(matched)} are in the catalogue, {len(unmatched)} are not")
print(f"catalogue: {len(published)} equal-mass orbits, {len(pub_in_range)} with period <= {t_max:.1f}; "
      f"we found {sum(r['index'] in found_idx for r in pub_in_range)} of those")
print("unmatched steady:", sum(o["growth_per_period"] < 1.02 for o in unmatched),
      " unmatched with out-of-plane speed |b| > 0.1:", sum(abs(o["params"][2]) > 0.1 for o in unmatched))
print("published |vz| quantiles (50/90/99/max):", np.round(np.quantile([abs(r["vz"]) for r in published], [0.5, 0.9, 0.99, 1]), 4))

# which published orbits lie on the family?
rows = np.load(os.path.join(RES, "family_trace.npy"))
rows = rows[np.argsort(rows[:, 0])]
curve = []
for r in rows:
    s0 = tb.isosceles_state(*r[1:5])
    E = tb.energy(s0)[0]
    curve.append((np.linalg.norm(tb.angular_momentum(s0)[0]) * abs(E) ** 0.5, 2 * r[5] * abs(E) ** 1.5))
curve = np.array(curve)
order = np.argsort(curve[:, 0])
on_family = []
for r in published:
    if not curve[:, 0].min() <= r["scaled_L"] <= curve[:, 0].max():
        continue
    phi = np.interp(r["scaled_L"], curve[order, 0], rows[order, 0])
    frac = Fraction(abs(phi) / np.pi).limit_denominator(80)
    z, miss = family.member(rows, -np.pi * frac.numerator / frac.denominator)
    if miss > 1e-10:
        continue
    s0 = tb.isosceles_state(*z[:4])
    E = tb.energy(s0)[0]
    if abs(np.linalg.norm(tb.angular_momentum(s0)[0]) * abs(E) ** 0.5 - r["scaled_L"]) > 2e-6:
        continue
    half_laps = r["scaled_period"] / (z[4] * abs(E) ** 1.5)
    if abs(half_laps * 2 - round(half_laps * 2)) > 1e-4:
        continue
    on_family.append({"index": r["index"], "period": r["period"], "turn": str(frac), "steady": r["steady"]})
print(f"published orbits on the family: {len(on_family)}")
print("  turns:", sorted({m["turn"] for m in on_family}, key=Fraction))
turns_published = {m["turn"] for m in on_family}
family_ours = json.load(open(os.path.join(RES, "exp8_family.json")))
status = {str(Fraction(m["turn"])): str(Fraction(m["turn"])) in turns_published
          for m in family_ours["made_to_order"] + family_ours["on_curve"]}
print("  our family members -> already published?", dict(sorted(status.items(), key=lambda kv: Fraction(kv[0]))))

json.dump({"ours": len(ours), "in_catalogue": len(matched), "not_in_catalogue": len(unmatched),
           "catalogue_size": len(published), "catalogue_in_our_period_range": len(pub_in_range),
           "of_those_we_found": sum(r["index"] in found_idx for r in pub_in_range),
           "period_range_max": t_max, "published_on_family": on_family, "family_member_published": status,
           "unmatched": [{"params": o["params"], "period": o["period"], "scaled_period": o["scaled_period"],
                          "scaled_L": o["scaled_L"], "growth": o["growth_per_period"], "built": bool(o.get("built")),
                          "turn": o.get("turn")} for o in unmatched],
           "matched": [{"period": o["period"], "published_index": o["published_index"]} for o in matched]},
          open(os.path.join(RES, "exp9_catalogue.json"), "w"), indent=1)
