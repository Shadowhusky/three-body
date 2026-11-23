"""Experiment 2: stability map of the collinear, zero-angular-momentum slice.

Initial conditions: equal masses on a line at x = -1, +1, 0; the outer two share
velocity (p1, p2) and the middle one has -2(p1, p2), so total momentum and
angular momentum vanish. Each pixel is one (p1, p2).

For every pixel we integrate the orbit together with a shadow orbit displaced by
1e-8, and record whether the system broke up, when, which body escaped, and how
fast the two orbits separated (a finite-time Lyapunov indicator).
"""
import os
import sys
import time

import numpy as np

import tb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
T_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0
P_MAX = 0.75

p = (np.arange(N) + 0.5) * P_MAX / N
P1, P2 = np.meshgrid(p, p)
s0 = np.array([tb.collinear_state(a, b) for a, b in zip(P1.ravel(), P2.ravel())])
t0 = time.time()
out, sym = tb.scan(s0, T_MAX)
print(f"{len(s0)} trajectories in {time.time() - t0:.0f}s")
np.savez_compressed(os.path.join(tb.HERE, "results", "map.npz"), out=out.reshape(N, N, 9), sym=sym.reshape(N, N, -1),
                    p=p, t_max=T_MAX, p_max=P_MAX)
