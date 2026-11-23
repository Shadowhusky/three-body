"""Experiment 1: validation of the integrator.

Checks against the analytic two-body solution, conservation of energy, momentum
and angular momentum, time reversibility, and agreement between two very
different accuracy settings.
"""
import json
import os

import numpy as np

import tb

OUT = os.path.join(tb.HERE, "results", "exp1_validation.json")
res = {}


def kepler_exact(a, e, t, mu):
    """Two-body relative orbit starting at pericentre, from Kepler's equation."""
    n = np.sqrt(mu / a**3)
    M = n * t
    E = M
    for _ in range(100):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return np.array([a * (np.cos(E) - e), a * np.sqrt(1 - e * e) * np.sin(E)])


# 1. Two bodies (third is massless and far away), eccentricity 0.9, 100 revolutions.
m = np.array([1.0, 1.0, 0.0])
a, e, mu = 1.0, 0.9, 2.0
rp = a * (1 - e)
vp = np.sqrt(mu * (1 + e) / (a * (1 - e)))
s0 = np.zeros(18)
s0[0:3] = [-rp / 2, 0, 0]
s0[3:6] = [rp / 2, 0, 0]
s0[6:9] = [0, 1e6, 0]
s0[10], s0[13] = -vp / 2, vp / 2
period = 2 * np.pi * np.sqrt(a**3 / mu)
errs = []
for revs in (1, 10, 100):
    t = revs * period * 0.6180339887
    s1 = tb.flow(s0, t, m)
    rel = s1[3:5] - s1[0:2]
    errs.append(float(np.linalg.norm(rel - kepler_exact(a, e, t, mu))))
res["kepler_e0.9_error_after_0.6_6_62_revs"] = errs

# 2. Conservation along a violent three-body trajectory.
s0 = tb.collinear_state(0.2, 0.31)
traj = tb.sample(s0, 40.0, 4001)
E = tb.energy(traj)
L = tb.angular_momentum(traj)
P = traj[:, 9:].reshape(-1, 3, 3).sum(axis=1)
res["chaotic_run"] = {
    "energy": float(E[0]),
    "max_rel_energy_drift": float(np.max(np.abs(E - E[0])) / abs(E[0])),
    "max_angular_momentum": float(np.max(np.abs(L))),
    "max_total_momentum": float(np.max(np.abs(P))),
}

# 3. Time reversal: forward T, then backward T, on an orderly orbit and on a chaotic one.
rev = {}
for name, s0, T in (("orderly", tb.collinear_state(0.347, 0.533), 60.0),
                    ("chaotic", tb.collinear_state(0.2, 0.31), 20.0)):
    back = tb.flow(tb.flow(s0, T), -T)
    rev[name] = float(np.linalg.norm(back - s0))
res["time_reversal_error"] = rev

# 4. Self-consistency: series order 16 / tol 1e-12 against order 30 / tol 1e-19.
s0 = tb.collinear_state(0.347, 0.533)
lo = tb.flow(s0, 50.0, order=16, tol=1e-12)
hi = tb.flow(s0, 50.0, order=30, tol=1e-19)
ref = tb.flow(s0, 50.0, order=24, tol=1e-17)
res["self_consistency_T50"] = {"coarse_vs_fine": float(np.linalg.norm(lo - hi)),
                               "default_vs_fine": float(np.linalg.norm(ref - hi))}

json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps(res, indent=1))
