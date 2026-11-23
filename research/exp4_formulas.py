"""Experiment 4: quasi-periodic fits in the regular regions.

Take initial conditions the map classifies as regular and that stay regular
over a much longer integration. Fit a sum of Fourier terms to a short span of
the orbit, then compare the extrapolated series with direct integration far
beyond the fitted span. Chaotic orbits get the same treatment as a control.
"""
import json
import os
from multiprocessing import Pool

import numpy as np

import formula
import tb

T_FIT, N_FIT, AHEAD = 600.0, 12001, 10
TERM_COUNTS = (30, 90)


def trial(args):
    p, ahead = args
    s0 = tb.collinear_state(*p)
    n_all = ahead * (N_FIT - 1) + 1
    t = np.linspace(0, T_FIT * ahead, n_all)
    path = tb.sample(s0, T_FIT * ahead, n_all)
    z = path[:, [0, 3, 6]] + 1j * path[:, [1, 4, 7]]
    out = {"params": [float(v) for v in p], "size": float(np.abs(z).max()), "fits": {}}
    for terms in TERM_COUNTS:
        rates, amps, _ = formula.fit(t[:N_FIT], z[:N_FIT], terms)
        miss = np.abs(formula.evaluate(rates, amps, t[::4]) - z[::4]).max(axis=1)
        tt = t[::4]
        out["fits"][terms] = {
            "rates": rates.tolist(),
            "amps": [[[float(a.real), float(a.imag)] for a in amps[:, b]] for b in range(3)],
            "missInFit": float(miss[tt <= T_FIT].max()),
            "missAt2x": float(miss[tt <= 2 * T_FIT].max()),
            "missAtEnd": float(miss.max()),
        }
    return out


def pick_starts(rng):
    d = np.load(os.path.join(tb.HERE, "results", "map.npz"))
    out, p = d["out"], d["p"]
    P1, P2 = np.meshgrid(p, p)
    tight = 3 * (P1**2 + P2**2) - 2.5 < -0.6  # skip weakly bound systems, which barely evolve within the survey time
    picks, verified = {}, {}
    for name, mask, t_long in (("orderly", (out[..., 0] == 0) & (out[..., 3] < 10) & tight, T_FIT * AHEAD),
                               ("chaotic", (out[..., 0] == 0) & (out[..., 3] > 60) & tight, T_FIT * 2)):
        cells = np.argwhere(mask)
        cells = cells[rng.permutation(len(cells))[:300]]
        starts = np.array([(p[i], p[j]) for j, i in cells])
        long_run, _ = tb.scan(np.array([tb.collinear_state(*q) for q in starts]), t_long, max_steps=8_000_000)
        ok = long_run[:, 0] == 0
        ok &= (long_run[:, 3] < 12) if name == "orderly" else (long_run[:, 3] > 100)
        picks[name] = starts[ok]
        verified[name] = [int(ok.sum()), len(starts), int(mask.sum())]
    return picks, verified


if __name__ == "__main__":
    rng = np.random.default_rng(4)
    picks, verified = pick_starts(rng)
    print("kept / tried / on map:", verified, flush=True)
    jobs = [(q, AHEAD) for q in picks["orderly"][:32]] + [(q, 2) for q in picks["chaotic"][:8]]
    with Pool(8) as pool:
        results = pool.map(trial, jobs, chunksize=1)
    orderly, chaotic = results[:len(jobs) - min(8, len(picks["chaotic"]))], results[-min(8, len(picks["chaotic"])):]

    best = sorted(orderly, key=lambda r: r["fits"][TERM_COUNTS[-1]]["missAtEnd"])
    shown = []
    for r in best[:5]:
        f = r["fits"][TERM_COUNTS[-1]]
        shown.append({"params": r["params"], "rates": f["rates"], "amps": f["amps"],
                      "stats": {"fitSpan": T_FIT, "testedTo": T_FIT * AHEAD, "size": r["size"],
                                "missInFit": f["missInFit"], "missAt2x": f["missAt2x"], "missAt10x": f["missAtEnd"]}})
    json.dump(shown, open(os.path.join(tb.HERE, "results", "formulas.json"), "w"))

    def table(group):
        return [{"params": r["params"], "size": r["size"],
                 **{f"terms{k}": {m: r["fits"][k][m] for m in ("missInFit", "missAt2x", "missAtEnd")} for k in TERM_COUNTS}}
                for r in group]

    summary = {"fitSpan": T_FIT, "testedTo": T_FIT * AHEAD, "verified": verified,
               "orderly": table(orderly), "chaotic": table(chaotic)}
    json.dump(summary, open(os.path.join(tb.HERE, "results", "exp4_formulas.json"), "w"), indent=1)
    for name, group in (("orderly", orderly), ("chaotic", chaotic)):
        for k in TERM_COUNTS:
            m = np.array([[r["fits"][k][x] / r["size"] for x in ("missInFit", "missAt2x", "missAtEnd")] for r in group])
            print(f"{name} n={len(m)} terms={k}: relative miss in-fit / 2x / end  median", np.round(np.median(m, axis=0), 5),
                  " best", np.round(m.min(axis=0), 6), " worst", np.round(m.max(axis=0), 3))
