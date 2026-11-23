"""Experiment 6: three-dimensional periodic orbits.

Initial conditions are isosceles triangles fixed by a reversing symmetry: a
half-turn about the triangle's axis combined with time reversal
(tb.isosceles_state). Four parameters define one: the apex offset y and the
velocity (a, b, c) of the first body. With c != 0 the angular momentum is
nonzero and the motion is genuinely three-dimensional.
"""
import json
import os
import sys

import numpy as np

import periodic
import tb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
T_HALF_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1
DTHRESH = float(sys.argv[4]) if len(sys.argv) > 4 else 0.2

if __name__ == "__main__":
    rng = np.random.default_rng(SEED)
    q = rng.uniform([0, -0.75, 0, 0.03], [0.9, 0.75, 1.3, 0.75], size=(int(N * 1.6), 4))
    y, a, b, c = q.T
    energy = 3 * a**2 + b**2 + 3 * c**2 - 0.5 - 2 / np.sqrt(1 + 9 * y**2)
    q = q[energy < -0.05][:N]
    orbits = periodic.search(1, q, T_HALF_MAX, dthresh=DTHRESH, cell=0.04, cache=f"spatial{SEED}", max_candidates=20000)
    spatial = [o for o in orbits if o["thickness"] > 1e-3]
    out = os.path.join(tb.HERE, "results", f"spatial_orbits_seed{SEED}.json")
    json.dump(spatial, open(out, "w"), indent=1)
    print(f"{len(spatial)} spatial orbits of {len(orbits)} written")
    for o in spatial[:60]:
        print("q=(" + ",".join(f"{v:.5f}" for v in o["params"]) + f") T={o['period']:.4f} T*={o['scaled_period']:.4f} "
              f"L*={o['scaled_L']:.4f} growth={o['growth_per_period']:.3g} thick={o['thickness']:.3f} "
              f"close={o['closure_error']:.1e}")
