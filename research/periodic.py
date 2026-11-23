"""Search for symmetric periodic orbits.

The equations of motion are reversible. An orbit that starts from a state fixed
by a reversing symmetry (sym_apply in core/tb.c) and reaches another such state
at time t is periodic with period 2t. We scan for near-returns to the symmetric
set, then refine with Newton's method until the residual is at rounding level.
"""
from multiprocessing import Pool

import numpy as np

import tb

MAKERS = {0: tb.collinear_state, 1: tb.isosceles_state}


def _residual(sym_type, z):
    s = tb.flow(MAKERS[sym_type](*z[:-1]), z[-1], order=20, tol=1e-15)
    return s - tb.sym_apply(s, sym_type)


def refine(sym_type, z0, iters=16, goal=3e-12):
    """Levenberg-Marquardt on (initial-condition parameters..., half period). Returns
    (z, residual). Gives up early when the residual is not decreasing, so poor candidates stay cheap.
    """
    z = np.array(z0, dtype=float)
    try:
        r = _residual(sym_type, z)
    except RuntimeError:
        return z, np.inf
    f = f_start = r @ r
    lam = 1e-4
    n = len(z)
    for it in range(iters):
        if f < goal * goal or (it == 5 and f > 0.05 * f_start):
            break
        J = np.empty((18, n))
        try:
            for k in range(n):
                dz = np.zeros(n)
                dz[k] = 1e-7
                J[:, k] = (_residual(sym_type, z + dz) - r) / 1e-7
        except RuntimeError:
            return z, np.inf
        g, H = J.T @ r, J.T @ J
        for _ in range(5):
            try:
                step = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-12), -g)
                zn = z + step
                if zn[-1] <= 0.05 or np.abs(step).max() > 0.1:
                    raise RuntimeError
                rn = _residual(sym_type, zn)
            except (RuntimeError, np.linalg.LinAlgError):
                lam *= 10
                continue
            if rn @ rn < f:
                z, r, f = zn, rn, rn @ rn
                lam = max(lam / 10, 1e-12)
                break
            lam *= 10
        else:
            break
    return z, float(np.sqrt(f))


def _refine_job(args):
    sym_type, z0 = args
    return refine(sym_type, z0)


def refine_many(sym_type, candidates, procs=tb.NCPU):
    with Pool(procs) as pool:
        return pool.map(_refine_job, [(sym_type, c) for c in candidates], chunksize=4)


def monodromy(s0, period):
    """Monodromy matrix (derivative of the period map) by central differences."""
    M = np.empty((18, 18))
    for k in range(18):
        d = np.zeros(18)
        d[k] = 1e-6
        M[:, k] = (tb.flow(s0 + d, period) - tb.flow(s0 - d, period)) / 2e-6
    return M


def describe(sym_type, z):
    """Summary of a converged orbit: invariants, closure error, stability, syzygy sequence."""
    s0 = MAKERS[sym_type](*z[:-1])
    period = 2 * z[-1]
    E = float(tb.energy(s0)[0])
    L = float(np.linalg.norm(tb.angular_momentum(s0)[0]))
    closure = float(np.linalg.norm(tb.flow(s0, period) - s0))
    eig = np.abs(np.linalg.eigvals(monodromy(s0, period)))
    path = tb.sample(s0, period, 1201)
    x = path[:, :9].reshape(-1, 3)
    sv = np.linalg.svd(x - x.mean(axis=0), compute_uv=False)
    out, sym = tb.scan(s0[None], period, order=24, tol=1e-16, threads=1)
    n_syz = int(out[0, 5])
    return {
        "sym": sym_type,
        "params": [float(v) for v in z[:-1]],
        "period": float(period),
        "energy": E,
        "L": L,
        "scaled_period": float(period * abs(E) ** 1.5),
        "scaled_L": float(L * abs(E) ** 0.5),
        "closure_error": closure,
        "growth_per_period": float(eig.max()),
        "min_sep": float(out[0, 4]),
        "thickness": float(sv[2] / sv[0]),
        "n_syz": n_syz,
        "word": "".join(str(c) for c in sym[0, :min(n_syz, tb.SYM_MAX)]),
        "state": [float(v) for v in s0],
    }


def _describe_job(args):
    try:
        return describe(*args)
    except (RuntimeError, np.linalg.LinAlgError):  # step budget exhausted, or a collision along the orbit
        return None


def describe_many(sym_type, zs, procs=tb.NCPU):
    with Pool(procs) as pool:
        return [d for d in pool.map(_describe_job, [(sym_type, z) for z in zs]) if d]


def dedupe(sym_type, results, miss_tol=1e-10, key_tol=2e-6):
    """Keep converged roots and remove duplicates.

    Two roots are the same orbit when their scale-invariant period and angular
    momentum agree; an orbit traversed k times is dropped in favour of the primitive one.
    """
    good = [(z, miss) for z, miss in results if miss < miss_tol]
    keyed = []
    for z, miss in good:
        s0 = MAKERS[sym_type](*z[:-1])
        E = tb.energy(s0)[0]
        if E >= 0:
            continue
        L = np.linalg.norm(tb.angular_momentum(s0)[0])
        keyed.append((2 * z[-1] * abs(E) ** 1.5, L * abs(E) ** 0.5, z))
    keyed.sort(key=lambda k: k[0])
    kept = []
    for tp, lp, z in keyed:
        dup = False
        for tq, lq, _ in kept:
            ratio = tp / tq
            if abs(lp - lq) < 1e-5 and abs(ratio - round(ratio)) < key_tol * max(1, ratio):
                dup = True
                break
        if not dup:
            kept.append((tp, lp, z))
    return [z for _, _, z in kept]


def cluster(starts, tm, dm, kp, cell, t_bin, dmax):
    """Near-returns form sheets: neighbouring initial conditions approach the symmetric
    set at similar times. Keep only the closest approach in each small cell of (parameters, time).
    """
    best = {}
    for i in np.nonzero(kp)[0]:
        for k in range(kp[i]):
            if dm[i, k] > dmax:
                continue
            key = (*np.floor(starts[i] / cell).astype(int), int(tm[i, k] / t_bin))
            if key not in best or dm[i, k] < best[key][0]:
                best[key] = (dm[i, k], (*starts[i], tm[i, k]))
    return [c for _, c in sorted(best.values())]


def search(sym_type, starts, t_half_max, dthresh, cell, cache, max_candidates=40000):
    """Full pipeline: scan for near-returns, refine, remove duplicates, describe.
    Each stage is cached under results/<cache>_*.npz so a rerun resumes.
    """
    import os
    import time

    scan_file = os.path.join(tb.HERE, "results", f"{cache}_scan.npz")
    if os.path.exists(scan_file):
        d = np.load(scan_file)
        starts, tm, dm, kp = d["starts"], d["tm"], d["dm"], d["kp"]
    else:
        s0 = np.array([MAKERS[sym_type](*q) for q in starts])
        t0 = time.time()
        tm, dm, kp = tb.symscan(s0, t_half_max, sym_type, dthresh=dthresh, kmax=6)
        print(f"scan: {len(s0)} starts in {time.time() - t0:.0f}s", flush=True)
        np.savez_compressed(scan_file, starts=starts, tm=tm, dm=dm, kp=kp)
    cands = cluster(starts, tm, dm, kp, cell, 0.3, dthresh)[:max_candidates]
    print(f"{int(kp.sum())} near-hits in {len(cands)} clusters to polish", flush=True)

    roots_file = os.path.join(tb.HERE, "results", f"{cache}_roots.npz")
    if os.path.exists(roots_file):
        d = np.load(roots_file)
        results = list(zip(d["z"], d["miss"]))
    else:
        t0 = time.time()
        results = refine_many(sym_type, cands)
        print(f"polish: {time.time() - t0:.0f}s", flush=True)
        np.savez_compressed(roots_file, z=np.array([z for z, _ in results]), miss=np.array([m for _, m in results]))
    print(f"{sum(m < 1e-10 for _, m in results)} converged", flush=True)
    roots = dedupe(sym_type, results)
    print(f"{len(roots)} distinct", flush=True)
    orbits = [o for o in describe_many(sym_type, roots) if o["closure_error"] < 1e-8]
    orbits.sort(key=lambda o: o["scaled_period"])
    return orbits
