"""Relative periodic orbits: solutions that are periodic in a rotating frame.

An orbit that starts from a symmetric isosceles state (tb.isosceles_state) and
reaches, at time t, a state that is symmetric about an axis rotated by phi
around the angular momentum vector repeats its shape every 2t, while the whole
configuration rotates by 2*phi. If phi = pi * j / k the orbit closes in the
inertial frame after k relative periods; otherwise it never closes, but one
relative period determines it for all time.
"""
import numpy as np

import tb


def spin_axis(s0):
    L = tb.angular_momentum(s0)[0]
    return L / np.linalg.norm(L)


def turn(s, axis, phi):
    """Rotate all positions and velocities about `axis` by phi."""
    k = np.asarray(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(phi) * K + (1 - np.cos(phi)) * K @ K
    return (np.asarray(s).reshape(6, 3) @ R.T).ravel()


def residual(z, phi, order=20, tol=1e-15):
    s0 = tb.isosceles_state(*z[:4])
    s = turn(tb.flow(s0, z[4], order=order, tol=tol), spin_axis(s0), -phi)
    return s - tb.sym_apply(s, 1)


def solve(z0, phi, iters=30, goal=1e-12):
    """Levenberg-Marquardt solve for (y, a, b, c, t) at fixed rotation angle phi."""
    z = np.array(z0, dtype=float)
    r = residual(z, phi)
    f, lam = r @ r, 1e-6
    for _ in range(iters):
        if f < goal * goal:
            break
        J = np.empty((18, 5))
        for k in range(5):
            dz = np.zeros(5)
            dz[k] = 1e-7
            J[:, k] = (residual(z + dz, phi) - r) / 1e-7
        g, H = J.T @ r, J.T @ J
        for _ in range(8):
            zn = z + np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-14), -g)
            rn = residual(zn, phi)
            if rn @ rn < f:
                z, r, f, lam = zn, rn, rn @ rn, max(lam / 10, 1e-14)
                break
            lam *= 10
        else:
            break
    return z, float(np.sqrt(f))


def best_phi(z):
    """Rotation angle at which a given (initial condition, t) is closest to a symmetric state."""
    phis = np.linspace(-np.pi, np.pi, 2881)
    miss = [np.linalg.norm(residual(z, p)) for p in phis]
    i = int(np.argmin(miss))
    return phis[i], miss[i]


def trace(z_start, direction, phi_end=np.pi, first_step=0.001 * np.pi, max_step=0.02 * np.pi, max_points=600):
    """Numerical continuation of the family from a known member: step the rotation
    angle and re-solve. Returns rows of (phi, y, a, b, c, t, residual).
    """
    zs, phis, rows = [np.array(z_start, dtype=float)], [0.0], []
    step = first_step
    while abs(phis[-1]) < phi_end and len(rows) < max_points:
        phi = phis[-1] + direction * step
        guess = zs[-1] if len(zs) < 2 else zs[-1] + (zs[-1] - zs[-2]) * (phi - phis[-1]) / (phis[-1] - phis[-2])
        try:
            z, miss = solve(guess, phi, iters=12)
        except (RuntimeError, np.linalg.LinAlgError):
            miss = np.inf
        if miss < 1e-10:
            zs.append(z)
            phis.append(phi)
            rows.append((phi, *z, miss))
            step = min(step * 1.5, max_step)
        else:
            step /= 3
            if step < 1e-6 * np.pi:
                break
    return np.array(rows)


def member(rows, phi):
    """Family member at a requested rotation angle, solved from the nearest traced one."""
    i = int(np.argmin(np.abs(rows[:, 0] - phi)))
    return solve(rows[i, 1:6], phi, iters=30)


def closing_time(z, max_laps):
    """Smallest multiple of the half relative period after which the orbit closes in the inertial frame."""
    s0 = tb.isosceles_state(*z[:4])
    best = (np.inf, None)
    for k in range(1, 4 * max_laps + 1):
        T = k * z[4] / 2
        miss = np.linalg.norm(tb.flow(s0, T) - s0)
        if miss < 1e-8:
            return T, miss
        best = min(best, (miss, T))
    return best[1], best[0]


def invariants(z):
    """Scale-invariant angular momentum L|E|^(1/2) and relative period 2t|E|^(3/2), plus the energy."""
    s0 = tb.isosceles_state(*z[:4])
    E = tb.energy(s0)[0]
    return np.linalg.norm(tb.angular_momentum(s0)[0]) * abs(E) ** 0.5, 2 * z[4] * abs(E) ** 1.5, E


def symmetric_epoch(q, period, samples=3000):
    """First time t > 0 at which the orbit through tb.isosceles_state(*q) is symmetric
    about an axis rotated by some angle phi around the angular momentum vector.
    Returns (t, phi, residual) with phi in (-pi/2, pi/2]; falls back to the half period.

    For a state s, the rotated-back state is linear in (1, cos phi, sin phi), so the
    best angle at each sampled time comes from a 3x3 quadratic form."""
    from scipy.optimize import least_squares

    s0 = tb.isosceles_state(*q)
    n = spin_axis(s0)
    path = tb.sample(s0, period / 2, samples + 1)[1:]
    times = np.linspace(0, period / 2, samples + 1)[1:]
    v = path.reshape(-1, 6, 3)
    along = (v @ n)[..., None] * n
    parts = [along, v - along, -np.cross(n, v)]  # coefficients of 1, cos phi, sin phi
    P = np.eye(18) - np.array([tb.sym_apply(e, 1) for e in np.eye(18)]).T
    comps = np.stack([p.reshape(-1, 18) @ P.T for p in parts], axis=1)  # (samples, 3, 18)
    gram = comps @ comps.transpose(0, 2, 1)
    phis = np.linspace(-np.pi / 2, np.pi / 2, 721)
    u = np.stack([np.ones_like(phis), np.cos(phis), np.sin(phis)])  # (3, angles)
    r2 = np.einsum("ia,sij,ja->sa", u, gram, u)
    best = r2.min(axis=1)
    dt = times[1] - times[0]
    for k in range(1, len(best) - 1):
        if best[k] < 0.05 and best[k] <= best[k - 1] and best[k] <= best[k + 1]:
            guess = [times[k], phis[int(np.argmin(r2[k]))]]
            try:
                fit = least_squares(lambda x: residual(np.r_[q, x[0]], x[1]), guess, xtol=1e-15, ftol=1e-15, gtol=1e-15,
                                    bounds=([times[k] - 2 * dt, -np.pi], [times[k] + 2 * dt, np.pi]))
            except RuntimeError:
                continue
            miss = np.linalg.norm(fit.fun)
            if miss < 1e-8:
                phi = (fit.x[1] + np.pi / 2) % np.pi - np.pi / 2
                return float(fit.x[0]), float(phi), float(miss)
    return period / 2, 0.0, float(np.linalg.norm(residual(np.r_[q, period / 2], 0.0)))


def _corrector(w_pred, tangent, iters=14):
    """Gauss-Newton onto the family, constrained to the plane through w_pred normal to the tangent."""
    def F(w):
        return np.r_[residual(w[:5], w[5]), (w - w_pred) @ tangent]

    w = w_pred.copy()
    r = F(w)
    for _ in range(iters):
        if np.linalg.norm(r) < 3e-13:
            break
        J = np.empty((19, 6))
        for k in range(6):
            d = np.zeros(6)
            d[k] = 1e-7
            J[:, k] = (F(w + d) - r) / 1e-7
        step = np.linalg.lstsq(J, -r, rcond=None)[0]
        for _ in range(6):
            trial = F(w + step)
            if np.linalg.norm(trial) < np.linalg.norm(r):
                w, r = w + step, trial
                break
            step = step / 2
        else:
            break
    return w, float(np.linalg.norm(r[:18]))


def trace_arclength(z_start, phi_start, max_points=250, ds_max=0.04, budget=1e9):
    """Pseudo-arclength continuation of the family through a known member, in both
    directions. Unlike `trace`, this follows the family around folds in phi.
    Returns rows of (y, a, b, c, t, phi, residual)."""
    import time

    t0 = time.time()
    w0 = np.r_[z_start, phi_start]
    z1, miss = solve(z_start, phi_start + 1e-3, iters=20)
    if miss > 1e-10:
        return np.array([[*w0, 0.0]])
    first_tangent = np.r_[z1, phi_start + 1e-3] - w0
    first_tangent /= np.linalg.norm(first_tangent)
    branches = []
    for direction in (1, -1):
        rows, w, tangent, ds = [], w0, direction * first_tangent, 2e-3
        while len(rows) < max_points and time.time() - t0 < budget:
            try:
                w_new, miss = _corrector(w + ds * tangent, tangent)
                ok = miss < 1e-10 and np.linalg.norm(w_new - w - ds * tangent) < 0.5 * ds
                if ok:
                    L, _, E = invariants(w_new)
                    ok = E < -1e-3 and w_new[4] > 0.05
            except (RuntimeError, np.linalg.LinAlgError):
                ok = False
            if not ok:
                ds /= 2
                if ds < 1e-6:
                    break
                continue
            tangent = (w_new - w) / np.linalg.norm(w_new - w)
            w = w_new
            rows.append([*w, miss])
            ds = min(ds * 1.4, ds_max)
            if L < 2e-3 or w[3] * w0[3] <= 0 or (len(rows) > 20 and np.linalg.norm(w - w0) < ds):
                break  # reached the collinear limit or a planar orbit, or the family closed on itself
        branches.append(rows)
    return np.array(branches[1][::-1] + [[*w0, 0.0]] + branches[0])
