"""Computer-assisted existence proofs for symmetric relative periodic orbits.

All computations use rigorous ball arithmetic (Arb, via python-flint): every
number carries an error radius that is guaranteed to contain the true value,
rounding errors included.

Statement proved for an orbit, given a numerical approximation
z0 = (y, a, b, c, t) and a rotation angle phi = pi * j / m:

  There is exactly one z* within rho of z0 such that the exact solution
  starting from tb.isosceles_state(y, a, b, c) is, at time t and after rotating
  it back by phi about its angular momentum vector, fixed by the reversing
  symmetry.

By reversibility such a solution satisfies s(u + 2t) = Rot(2 phi) s(u) for all
u, so it is periodic in the inertial frame after m relative periods (see the
lemma in `certify`).

Method. The five numbers z are the unknowns of five equations G(z) = 0. The
Krawczyk test shows that a Newton-like operator maps the box Z = z0 +/- rho
strictly into itself, which implies a unique zero in Z. It needs G(z0) and an
enclosure of the Jacobian DG over all of Z. Both come from validated
integration: a Taylor-series method whose truncation error is bounded on a
rigorous a priori enclosure of the solution.

Numbers are "duals", [value, d/dy, d/da, d/db, d/dc]: forward-mode derivatives
with respect to the initial-condition parameters.
"""
from flint import arb, arb_mat, ctx

PAIRS = ((0, 1), (0, 2), (1, 2))
ZERO = None  # set by set_precision


def set_precision(bits):
    global ZERO
    ctx.prec = bits
    ZERO = [arb(0)] * 5


def dadd(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4]]


def dsub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3], a[4] - b[4]]


def dmul(a, b):
    a0, b0 = a[0], b[0]
    return [a0 * b0, a0 * b[1] + a[1] * b0, a0 * b[2] + a[2] * b0, a0 * b[3] + a[3] * b0, a0 * b[4] + a[4] * b0]


def dscale(a, c):
    return [a[0] * c, a[1] * c, a[2] * c, a[3] * c, a[4] * c]


def jet(x, order):
    """Taylor coefficients in time of the solution through x (18 duals: positions,
    then velocities). Returns X, V with X[c][k] the k-th coefficient of coordinate c."""
    X = [[x[c]] for c in range(9)]
    V = [[x[9 + c]] for c in range(9)]
    D = [[[], [], []] for _ in PAIRS]
    S = [[] for _ in PAIRS]
    U = [[] for _ in PAIRS]
    inv_s0 = [None] * 3
    for k in range(order):
        A = [ZERO] * 9
        for p, (i, j) in enumerate(PAIRS):
            Dp = D[p]
            for d in range(3):
                Dp[d].append(dsub(X[3 * j + d][k], X[3 * i + d][k]))
            s = ZERO
            for d in range(3):
                Dd = Dp[d]
                for l in range((k + 1) // 2):
                    s = dadd(s, dscale(dmul(Dd[l], Dd[k - l]), 2))
                if k % 2 == 0:
                    s = dadd(s, dmul(Dd[k // 2], Dd[k // 2]))
            S[p].append(s)
            if k == 0:
                v = s[0]
                r = 1 / (v * v.sqrt())          # s^(-3/2); nan if the pair may collide
                g = r / v * arb(-1.5)
                U[p].append([r, g * s[1], g * s[2], g * s[3], g * s[4]])
                w = -1 / (v * v)
                inv_s0[p] = [1 / v, w * s[1], w * s[2], w * s[3], w * s[4]]
            else:
                acc = ZERO
                Sp, Up = S[p], U[p]
                for q in range(1, k + 1):
                    acc = dadd(acc, dscale(dmul(Sp[q], Up[k - q]), arb(-3 * q - 2 * (k - q)) / 2))
                U[p].append(dscale(dmul(acc, inv_s0[p]), arb(1) / k))
            Up = U[p]
            for d in range(3):
                Dd = Dp[d]
                f = ZERO
                for l in range(k + 1):
                    f = dadd(f, dmul(Dd[l], Up[k - l]))
                A[3 * i + d] = dadd(A[3 * i + d], f)
                A[3 * j + d] = dsub(A[3 * j + d], f)
        inv = arb(1) / (k + 1)
        for c in range(9):
            X[c].append(dscale(V[c][k], inv))
            V[c].append(dscale(A[c], inv))
    return X, V


def field(x):
    X, V = jet(x, 1)
    return [X[c][1] for c in range(9)] + [V[c][1] for c in range(9)]


def taylor_sum(X, V, h, upto):
    out = []
    for series in X + V:
        acc = series[upto]
        for k in range(upto - 1, -1, -1):
            acc = dadd(dscale(acc, h), series[k])
        out.append(acc)
    return out


def step_size(X, V, order, tol):
    """Heuristic only: any positive step is valid; a poor one merely fails to verify."""
    h = 1e300
    for k in (order - 1, order):
        big = max(float(abs(s[k][0].mid())) for s in X + V)
        if big > 0:
            h = min(h, (tol / big) ** (1.0 / k))
    return 0.9 * h


def enclose(x, h):
    """A priori enclosure: a box that provably contains the solution (and its parameter
    derivatives) for all times in [0, h] and all initial conditions in x, or None if
    none was verified. Test: x + [0, h] * f(B) inside B implies the solution stays in B.
    """
    span = h * arb(0.5, 0.5)
    pad = arb(2) ** (-ctx.prec + 16)

    def image(box):
        f = field(box)
        return [[xc[i] + span * fc[i] for i in range(5)] for xc, fc in zip(x, f)]

    cand = image(x)
    for _ in range(12):
        if not all(c.is_finite() for cc in cand for c in cc):
            return None
        # inflate the latest image slightly (inflating the previous box compounds and never converges)
        B = [[arb(c.mid(), c.rad() * 1.1 + pad) for c in cc] for cc in cand]
        cand = image(B)
        if all(c.is_finite() and b.contains_interior(c) for bc, cc in zip(B, cand) for b, c in zip(bc, cc)):
            return B
    return None


def rigorous_step(x, h, order):
    """Validated step: enclosure of the solution after time h (h a positive ball). The
    truncation error is the next Taylor coefficient evaluated over the a priori enclosure B.
    """
    B = enclose(x, h.upper())
    if B is None:
        return None
    X, V = jet(x, order - 1)
    XB, VB = jet(B, order)
    out = taylor_sum(X, V, h, order - 1)
    hn = h ** order
    return [dadd(o, dscale(series[order], hn)) for o, series in zip(out, XB + VB)]


def start(z):
    """tb.isosceles_state as duals, seeded with the parameter derivatives."""
    y, a, b, c = ([z[i]] + [arb(1 if j == i else 0) for j in range(4)] for i in range(4))
    const = lambda v: [arb(v)] + [arb(0)] * 4
    neg = lambda u: dscale(u, arb(-1))
    return [const(-1), y, const(0), const(1), y, const(0), const(0), dscale(y, arb(-2)), const(0),
            a, b, c, a, neg(b), c, dscale(a, arb(-2)), const(0), dscale(c, arb(-2))]


def flow(z, order, tol, rigorous):
    """State (18 duals) at time z[4] for initial-condition parameters z[:4], plus its
    time derivative. With rigorous=True and z a box, the result encloses every
    parameter choice in the box and every time in z[4].
    """
    x = start(z)
    t_goal = z[4].lower() if rigorous else z[4]
    t = arb(0)
    steps = 0
    while True:
        X, V = jet(x, order)
        h = step_size(X, V, order, tol)
        left = t_goal - t
        last = not (left > arb(h))
        hh = left if last else arb(h)
        if rigorous:
            new = rigorous_step(x, hh, order)
            tries = 0
            while new is None:
                tries += 1
                if tries > 40:
                    raise ArithmeticError("no verified enclosure")
                h, last = h / 2, False
                hh = arb(h)
                new = rigorous_step(x, hh, order)
            x = new
        else:
            x = taylor_sum(X, V, hh, order)
        t = t + hh
        steps += 1
        if last:
            break
    if rigorous:
        # cover the whole time interval z[4]: enclose the solution over [0, width]
        width = (z[4].upper() - t).upper()
        if width > 0:
            x = enclose(x, width)
            if x is None:
                raise ArithmeticError("no verified enclosure for the final time ball")
    if not all(v.is_finite() for xc in x for v in xc):
        raise ArithmeticError("enclosure blew up")
    return x, field(x), steps


def spin_axis(x0):
    L = [ZERO, ZERO, ZERO]
    for b in range(3):
        r, v = x0[3 * b:3 * b + 3], x0[9 + 3 * b:12 + 3 * b]
        L[0] = dadd(L[0], dsub(dmul(r[1], v[2]), dmul(r[2], v[1])))
        L[1] = dadd(L[1], dsub(dmul(r[2], v[0]), dmul(r[0], v[2])))
        L[2] = dadd(L[2], dsub(dmul(r[0], v[1]), dmul(r[1], v[0])))
    n2 = dadd(dadd(dmul(L[0], L[0]), dmul(L[1], L[1])), dmul(L[2], L[2]))
    r = 1 / n2[0].sqrt()
    g = r / n2[0] * arb(-0.5)
    inv_norm = [r, g * n2[1], g * n2[2], g * n2[3], g * n2[4]]
    return [dmul(c, inv_norm) for c in L]


def turn_matrix(n, sin_phi, cos_phi):
    """Rotation about unit vector n (duals) by the angle with the given sine and cosine."""
    one = [arb(1)] + [arb(0)] * 4
    K = [[ZERO, dscale(n[2], arb(-1)), n[1]], [n[2], ZERO, dscale(n[0], arb(-1))], [dscale(n[1], arb(-1)), n[0], ZERO]]
    K2 = [[dadd(dadd(dmul(K[i][0], K[0][j]), dmul(K[i][1], K[1][j])), dmul(K[i][2], K[2][j])) for j in range(3)]
          for i in range(3)]
    return [[dadd(dadd(one if i == j else ZERO, dscale(K[i][j], sin_phi)), dscale(K2[i][j], 1 - cos_phi))
             for j in range(3)] for i in range(3)]


def equations(z, turn, order, tol, rigorous):
    """G(z) (5 balls), its 5x5 Jacobian, and the first coordinate of the rotated state.
    turn = (j, m) means phi = pi * j / m; the state is rotated back by phi.
    """
    j, m = turn
    sin_phi, cos_phi = (arb.pi() * j / m).sin_cos()
    x, xdot, steps = flow(z, order, tol, rigorous)
    R = turn_matrix(spin_axis(start(z)), -sin_phi, cos_phi)

    def combos(s, M):
        def turned(vec):
            return [dadd(dadd(dmul(M[i][0], vec[0]), dmul(M[i][1], vec[1])), dmul(M[i][2], vec[2])) for i in range(3)]

        x1, x2, v1, v2 = (turned(s[k:k + 3]) for k in (0, 3, 9, 12))
        return [dadd(x1[0], x2[0]), dadd(x1[2], x2[2]), dsub(x1[1], x2[1]), dsub(v1[0], v2[0]), dadd(v1[1], v2[1])], x1[0]

    G, x1x = combos(x, R)
    # time derivative: the rotation does not depend on time, so apply it to the vector field
    plain = lambda u: [u[0]] + [arb(0)] * 4
    Gt, _ = combos([plain(v) for v in xdot], [[plain(e) for e in row] for row in R])
    values = [g[0] for g in G]
    jac = [[G[i][1], G[i][2], G[i][3], G[i][4], Gt[i][0]] for i in range(5)]
    return values, jac, x1x[0], steps


def newton(z, turn, order, tol, rounds=6, goal=None):
    """Refine a numerical approximation (no rigour needed here; the proof comes afterwards)."""
    z = [arb(v) for v in z]
    goal = goal or arb(2) ** (-ctx.prec // 2)
    size = None
    for _ in range(rounds):
        values, jac, _, _ = equations(z, turn, order, tol, rigorous=False)
        size = max(abs(v.mid()) for v in values)
        if size < goal:
            break
        J = arb_mat(5, 5, [arb(e.mid()) for row in jac for e in row])
        dz = J.solve(arb_mat(5, 1, [arb(v.mid()) for v in values]))
        z = [arb((z[i] - dz[i, 0]).mid()) for i in range(5)]
    return z, size


def certify(z0, turn, order, tol):
    """Krawczyk test on a box z0 +/- rho, with rho chosen a little larger than the
    Newton correction at z0. Returns a report dict; report["proved"] is True only
    if every rigorous check passed.

    Lemma (why five equations suffice). Write the rotated state's body-1 and
    body-2 positions as u, w and velocities as p, q. Being symmetric means six
    relations: u_x + w_x = 0, u_z + w_z = 0, u_y = w_y, p_x = q_x, p_y + q_y = 0,
    p_z = q_z. The angular momentum component along the rotated axis is exactly
    zero for these initial conditions and is conserved; given the first two
    relations it equals u_z (p_x - q_x) - u_x (p_z - q_z). So the first five
    relations plus u_x != 0 imply the sixth.
    """
    z0 = [arb(v.mid()) for v in z0]
    g0, jac0, _, steps0 = equations(z0, turn, order, tol, rigorous=True)
    C = arb_mat(5, 5, [arb(e.mid()) for row in jac0 for e in row]).inv()
    C = arb_mat(5, 5, [arb(C[i, k].mid()) for i in range(5) for k in range(5)])
    correction = C * arb_mat(5, 1, g0)
    reach = max(abs(correction[i, 0].mid()) + correction[i, 0].rad() for i in range(5))
    rho_exp = int(float((reach * 64).log() / arb(2).log())) + 1
    rho = arb(2) ** rho_exp
    Z = [arb(v.mid(), rho) for v in z0]
    _, jacZ, x1x, stepsZ = equations(Z, turn, order, tol, rigorous=True)
    DG = arb_mat(5, 5, [e for row in jacZ for e in row])
    shrink = arb_mat(5, 5, [arb(1 if i == k else 0) for i in range(5) for k in range(5)]) - C * DG
    box = arb_mat(5, 1, [arb(0, rho)] * 5)
    K = -correction + shrink * box  # relative to z0
    inside = all(arb(0, rho).contains_interior(K[i, 0]) for i in range(5))
    y, c = Z[0], Z[3]
    report = {
        "proved": bool(inside and not x1x.contains(0) and not y.contains(0) and not c.contains(0)),
        "krawczyk_inside": bool(inside),
        "image_radius_over_rho": float(max((abs(K[i, 0].mid()) + K[i, 0].rad()) for i in range(5)) / rho),
        "contraction": float(max(sum((abs(shrink[i, k].mid()) + shrink[i, k].rad()) for k in range(5)) for i in range(5))),
        "residual_at_centre": float(max(abs(v.mid()) + v.rad() for v in g0)),
        "u_x_nonzero": bool(not x1x.contains(0)),
        "leaves_the_plane": bool(not y.contains(0) and not c.contains(0)),
        "steps": [steps0, stepsZ],
        "rho": f"2^{rho_exp}",
        "bits": ctx.prec,
        "order": order,
        "centre": [v.str(40, radius=False) for v in z0],
    }
    return report
