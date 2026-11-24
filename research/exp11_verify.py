"""Experiment 11: independent checks on the earlier results.

  1. Validated integrator against an exact solution. Lagrange's rotating
     equilateral triangle is known in closed form; the rigorous enclosure
     must contain it, and its parameter derivatives must agree with finite
     differences of the ordinary integrator.
  2. Catalogue matching. Matched orbits should agree with the catalogue to many
     digits and unmatched ones should be far from every entry, with nothing
     in between.
  3. Stability labels. For orbits that are in the catalogue, our stable /
     unstable label should agree with the published one.
  4. A proof repeated at higher precision and series order must certify the
     same solution.
"""
import glob
import json
import os

import numpy as np
from flint import arb

import proof
import tb

RES = os.path.join(tb.HERE, "results")
out = {}

# 1. Lagrange solution: side 2, angular velocity sqrt(3/8), in the isosceles parametrisation
proof.set_precision(256)
omega = (arb(3) / 8).sqrt()
y = 1 / arb(3).sqrt()
z = [y, -omega * y, -omega, arb(0), arb(3)]
x, _, steps = proof.flow(z, order=30, tol=1e-45, rigorous=True)
start = proof.start(z)
c, s = (omega * z[4]).sin_cos()[::-1]
worst_radius, contained = 0.0, True
for body in range(6):  # three positions, then three velocities
    x0, y0, z0 = (start[3 * body + d][0] for d in range(3))
    exact = [c * x0 - s * y0, s * x0 + c * y0, z0]
    for d in range(3):
        ball = x[3 * body + d][0]
        contained &= bool(ball.overlaps(exact[d]))
        worst_radius = max(worst_radius, float(ball.rad()))
zf = np.array([float(v.mid()) for v in z])
worst_derivative = 0.0
for k in range(4):
    dz = np.zeros(5)
    dz[k] = 1e-6
    fd = (tb.flow(tb.isosceles_state(*(zf + dz)[:4]), zf[4]) - tb.flow(tb.isosceles_state(*(zf - dz)[:4]), zf[4])) / 2e-6
    ours = np.array([float(x[i][k + 1].mid()) for i in range(18)])
    worst_derivative = max(worst_derivative, float(np.abs(fd - ours).max()))
out["lagrange"] = {"steps": steps, "enclosure_contains_exact_solution": contained,
                   "largest_enclosure_radius": worst_radius, "derivative_vs_finite_difference": worst_derivative}

# 2. catalogue matching: distance of every orbit of ours to its nearest catalogue entry
published = json.load(open(os.path.join(RES, "published_3d_equal_mass.json")))
pub = []
for r in published:
    st = np.array([-1, 0, 0, 1, 0, 0, 0, 0, r["z0"], r["vx"], r["vy"], r["vz"], r["vx"], r["vy"], -r["vz"],
                   -2 * r["vx"], -2 * r["vy"], 0], dtype=float)
    E = tb.energy(st)[0]
    pub.append((r["period"] * abs(E) ** 1.5, np.linalg.norm(tb.angular_momentum(st)[0]) * abs(E) ** 0.5))
pub = np.array(pub)
ours = []
for path in sorted(glob.glob(os.path.join(RES, "spatial_orbits_seed*.json"))):
    ours += json.load(open(path))
ours += [o for o in json.load(open(os.path.join(RES, "family_orbits.json"))) if not o["already_in_search"]]
gaps = []
for o in ours:
    best = np.inf
    for k in (1, 2, 3, 4, 1 / 2, 1 / 3, 1 / 4):
        best = min(best, np.hypot((o["scaled_period"] * k - pub[:, 0]) / pub[:, 0], o["scaled_L"] - pub[:, 1]).min())
    gaps.append(best)
gaps = np.array(gaps)
out["catalogue_matching"] = {
    "orbits": len(gaps),
    "matched_largest_gap": float(gaps[gaps < 1e-5].max()), "matched": int((gaps < 1e-5).sum()),
    "unmatched_smallest_gap": float(gaps[gaps >= 1e-5].min()), "unmatched": int((gaps >= 1e-5).sum()),
    "between_1e-7_and_1e-4": int(((gaps > 1e-7) & (gaps < 1e-4)).sum()),
}

# 3. stability labels of matched orbits
agree = total = 0
for o, gap in zip(ours, gaps):
    if gap < 1e-5:
        k = int(np.argmin(np.hypot((o["scaled_period"] - pub[:, 0]) / pub[:, 0], o["scaled_L"] - pub[:, 1])))
        if abs(o["scaled_period"] - pub[k, 0]) < 1e-5 * pub[k, 0]:
            total += 1
            agree += (o["growth_per_period"] < 1.02) == published[k]["steady"]
out["stability_labels"] = {"compared": total, "agree": int(agree)}

# 4. one proof repeated at higher precision
target = next(p for p in json.load(open(os.path.join(RES, "proofs.json"))) if p["proved"] and "period 7.4208" in p["label"])
orbit = next(o for o in ours if abs(o["period"] - 7.4208) < 1e-3)
z_double = list(np.r_[orbit["params"], orbit["period"] / 2])
proof.set_precision(512)
zz, _ = proof.newton(z_double, (0, 1), order=60, tol=1e-95)
again = proof.certify(zz, (0, 1), order=60, tol=1e-95)
first = [arb(v) for v in target["centre"]]
out["proof_repeated"] = {
    "proved_again": again["proved"], "rho_first": target["rho"], "rho_again": again["rho"],
    "centres_differ_by": float(max(abs(arb(a) - b).mid() for a, b in zip(again["centre"], first))),
}

json.dump(out, open(os.path.join(RES, "exp11_verify.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
