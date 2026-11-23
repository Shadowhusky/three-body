"""Experiment 8: a family of relative periodic orbits.

Continue one stable 3D orbit from the search as a relative periodic orbit, i.e.
one that is periodic in a rotating frame (family.py), with the rotation angle
as the continuation parameter. Then:
  1. check which stable 3D orbits found by the search lie on this family, and
     whether their rotation angles are the rational multiples of pi required
     for them to close in the inertial frame;
  2. compute members at rational angles the search did not find, and check
     that they close in the inertial frame;
  3. take a member with an irrational angle (it never closes) and check that
     one relative period plus the rotation predicts it far ahead.

The family was first computed by Michalodimitrakis (1981, A&A 93, 212).
"""
import glob
import json
import os
from fractions import Fraction

import numpy as np

import family
import periodic
import tb

found = []
for f in sorted(glob.glob(os.path.join(tb.HERE, "results", "spatial_orbits_seed*.json"))):
    found += json.load(open(f))
seed = min(found, key=lambda o: o["scaled_L"] if o["growth_per_period"] < 1.02 else np.inf)
z0 = np.array([*seed["params"], seed["period"] / 2])
print("seed:", seed["params"], "period", seed["period"])

rows = np.vstack([family.trace(z0, +1)[::-1], [[0.0, *z0, 0.0]], family.trace(z0, -1)])
np.save(os.path.join(tb.HERE, "results", "family_trace.npy"), rows)
print(f"traced {len(rows)} members, turning angle/pi from {rows[0, 0] / np.pi:.4f} to {rows[-1, 0] / np.pi:.4f}")


def invariants(z):
    s0 = tb.isosceles_state(*z[:4])
    E = tb.energy(s0)[0]
    return np.linalg.norm(tb.angular_momentum(s0)[0]) * abs(E) ** 0.5, 2 * z[4] * abs(E) ** 1.5


curve = np.array([invariants(r[1:6]) for r in rows])  # (scale-invariant angular momentum, relative period)

# 1. which stable orbits from the search lie on the family?
on_curve = []
for o in found:
    if o["growth_per_period"] > 1.02 or not curve[:, 0].min() <= o["scaled_L"] <= curve[:, 0].max():
        continue
    i = int(np.argmin(np.abs(curve[:, 0] - o["scaled_L"])))
    start_gap = np.linalg.norm(np.array(o["params"]) - rows[i, 1:5])
    if start_gap > 0.02:
        continue
    phi = np.interp(o["scaled_L"], curve[::-1, 0], rows[::-1, 0]) if curve[0, 0] > curve[-1, 0] else np.interp(o["scaled_L"], curve[:, 0], rows[:, 0])
    frac = Fraction(abs(phi) / np.pi).limit_denominator(40)
    z, miss = family.member(rows, -np.pi * frac.numerator / frac.denominator)
    same = np.linalg.norm(z[:4] - np.array(o["params"]))
    on_curve.append({"params": o["params"], "period": o["period"], "scaled_L": o["scaled_L"],
                     "turn": f"{frac.numerator}/{frac.denominator}", "half_laps_to_close": o["period"] / z[4],
                     "start_differs_by": float(same), "family_miss": miss})
    print(f"found orbit T={o['period']:.4f}: turn {frac} of a half-turn per half-lap, "
          f"closes after {o['period'] / z[4]:.3f} half-laps, start matches family member to {same:.1e}")

# 2. members computed at prescribed rational angles
known = {m["turn"] for m in on_curve}
made = []
for frac in [Fraction(1, 2), Fraction(1, 4), Fraction(2, 5), Fraction(3, 5), Fraction(1, 6), Fraction(2, 7),
             Fraction(3, 7), Fraction(3, 8), Fraction(4, 9), Fraction(1, 13), Fraction(5, 8), Fraction(1, 25)]:
    if str(frac) in known:
        continue
    z, miss = family.member(rows, -np.pi * frac.numerator / frac.denominator)
    if miss > 1e-10:
        continue
    T, closure = family.closing_time(z, frac.denominator)
    if closure > 1e-8:
        continue
    o = periodic.describe(1, np.array([*z[:4], T / 2]))
    o["turn"] = str(frac)
    o["already_in_search"] = any(abs(o["scaled_period"] - f["scaled_period"]) < 1e-4 * o["scaled_period"] and
                                 abs(o["scaled_L"] - f["scaled_L"]) < 1e-4 for f in found)
    made.append(o)
    print(f"made to order, turn {frac}: period {T:.4f}, closes to {closure:.1e}, growth per period "
          f"{o['growth_per_period']:.3f}, search had it: {o['already_in_search']}")

# 3. an irrational member: never closes, yet determined for all time by one relative period
phi_irr = -np.pi * (np.sqrt(5) - 1) / 4  # golden ratio / 2 of a half-turn
z, miss = family.member(rows, phi_irr)
s0 = tb.isosceles_state(*z[:4])
axis = family.spin_axis(s0)
lap = 2 * z[4]
worst = 0.0
for laps in (10, 100, 1000):
    for part in (0.0, 0.37, 0.81):
        # one relative period of motion, rotated by the angle accumulated so far
        predicted = family.turn(tb.flow(s0, part * lap), axis, 2 * phi_irr * laps)
        worst = max(worst, np.linalg.norm(tb.flow(s0, (laps + part) * lap) - predicted))
    print(f"irrational member after {laps} laps: 'one lap + turn' vs direct simulation differ by {worst:.1e}")

json.dump({"seed": seed["params"], "points": len(rows),
           "turn_range": [float(rows[:, 0].min() / np.pi), float(rows[:, 0].max() / np.pi)],
           "spin_range": [float(curve[:, 0].min()), float(curve[:, 0].max())],
           "on_curve": on_curve, "made_to_order": [{k: m[k] for k in ("turn", "period", "params", "closure_error",
                                                                      "growth_per_period", "already_in_search")} for m in made],
           "irrational": {"turn_over_pi": float(phi_irr / np.pi), "family_miss": miss, "laps_checked": 1000,
                          "worst_gap": float(worst)}},
          open(os.path.join(tb.HERE, "results", "exp8_family.json"), "w"), indent=1)
json.dump(made, open(os.path.join(tb.HERE, "results", "family_orbits.json"), "w"), indent=1)
