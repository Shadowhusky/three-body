"""Experiment 7: lifetime distribution of unstable triples.

Prediction. A body ejected with slightly negative energy e returns after one
Kepler period, P ~ |e|^(-3/2). If the distribution of e is smooth across e = 0,
then P(|e| < x) ~ x, so the survival probability should decay as t^(-2/3):
slowly enough that the mean lifetime diverges.
"""
import json
import os
import sys

import numpy as np

import tb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
T_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 3000.0

rng = np.random.default_rng(7)
q = rng.uniform(0, 0.9, size=(4 * N, 2))
q = q[3 * (q**2).sum(axis=1) < 2.4][:N]
s0 = np.array([tb.collinear_state(a, b) for a, b in q])
out, _ = tb.scan(s0, T_MAX, max_steps=6_000_000)
np.savez_compressed(os.path.join(tb.HERE, "results", "lifetimes.npz"), q=q, out=out, t_max=T_MAX)

status, t_end, growth = out[:, 0], out[:, 1], out[:, 3]
orderly = (status == 0) & (growth < 20)
doomed = ~orderly
ts = np.geomspace(3, T_MAX * 0.999, 40)
surv = np.array([np.mean((t_end[doomed] > t) | (status[doomed] != 1)) for t in ts])
fit = (ts > 100) & (surv > 0)
slope, icpt = np.polyfit(np.log(ts[fit]), np.log(surv[fit]), 1)
res = {"n": int(len(q)), "t_max": T_MAX, "orderly_fraction": float(orderly.mean()),
       "budget_exhausted": int((status == 2).sum()), "times": ts.tolist(), "surviving": surv.tolist(),
       "tail_exponent": float(slope), "predicted_exponent": -2 / 3}
json.dump(res, open(os.path.join(tb.HERE, "results", "exp7_lifetimes.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("times", "surviving")}, indent=1))
print(np.c_[ts, surv][::4])
