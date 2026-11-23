"""Experiment 3: periodic orbits in the planar, zero-angular-momentum slice.

Scan the same (p1, p2) slice as the map for orbits that nearly return to a
symmetric state, refine the candidates with Newton's method, and keep those
that close to rounding error.
"""
import json
import os
import sys

import numpy as np

import periodic
import tb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
T_HALF_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
P_MAX = 0.75

if __name__ == "__main__":
    p = (np.arange(N) + 0.5) * P_MAX / N
    P1, P2 = np.meshgrid(p, p)
    q = np.c_[P1.ravel(), P2.ravel()]
    q = q[3 * (q**2).sum(axis=1) < 2.5]  # negative total energy only
    orbits = periodic.search(0, q, T_HALF_MAX, dthresh=0.06, cell=3 * P_MAX / N, cache=f"planar{N}")
    json.dump(orbits, open(os.path.join(tb.HERE, "results", "planar_orbits.json"), "w"), indent=1)
    print(f"{len(orbits)} orbits written")
    for o in orbits:
        print(f"p=({o['params'][0]:.6f},{o['params'][1]:.6f}) T={o['period']:.4f} T*={o['scaled_period']:.4f} "
              f"growth={o['growth_per_period']:.3g} close={o['closure_error']:.1e} word={o['word'][:24]}")
