"""Python bindings for the C integrator (core/tb.c)."""
import ctypes
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "core", "tb.c")
LIB = os.path.join(HERE, "core", "libtb.dylib")
SYM_MAX = 96
NCPU = os.cpu_count() or 8


def _build():
    if not os.path.exists(LIB) or os.path.getmtime(LIB) < os.path.getmtime(SRC):
        tmp = f"{LIB}.{os.getpid()}"
        subprocess.check_call(["clang", "-O3", "-ffp-contract=off", "-shared", "-fPIC", "-o", tmp, SRC])
        os.replace(tmp, LIB)


_build()
_lib = ctypes.CDLL(LIB)
_dp = ctypes.POINTER(ctypes.c_double)
_up = ctypes.POINTER(ctypes.c_ubyte)
_lib.tb_flow.restype = ctypes.c_long
_lib.tb_flow.argtypes = [_dp, _dp, ctypes.c_double, _dp, ctypes.c_int, ctypes.c_double, ctypes.c_long]
_lib.tb_sample.restype = ctypes.c_long
_lib.tb_sample.argtypes = [_dp, _dp, ctypes.c_double, ctypes.c_int, _dp, ctypes.c_int, ctypes.c_double, ctypes.c_long]
_lib.tb_scan_batch.restype = None
_lib.tb_scan_batch.argtypes = [_dp, ctypes.c_int, _dp, ctypes.c_double, ctypes.c_int, ctypes.c_double,
                               ctypes.c_double, ctypes.c_double, ctypes.c_long, _dp, _up]

EQUAL = np.ones(3)


def _p(a):
    return a.ctypes.data_as(_dp)


def flow(s0, t, m=EQUAL, order=24, tol=1e-17, max_steps=5_000_000):
    """State after time t (t may be negative)."""
    s0 = np.ascontiguousarray(s0, dtype=float)
    m = np.ascontiguousarray(m, dtype=float)
    s1 = np.empty(18)
    n = _lib.tb_flow(_p(s0), _p(m), float(t), _p(s1), order, tol, max_steps)
    if n < 0:
        raise RuntimeError("step budget exhausted")
    return s1


def sample(s0, t_end, n_out, m=EQUAL, order=24, tol=1e-17, max_steps=20_000_000):
    """Trajectory sampled at n_out evenly spaced times in [0, t_end]; shape (n_out, 18)."""
    s0 = np.ascontiguousarray(s0, dtype=float)
    m = np.ascontiguousarray(m, dtype=float)
    out = np.empty((n_out, 18))
    n = _lib.tb_sample(_p(s0), _p(m), float(t_end), n_out, _p(out), order, tol, max_steps)
    if n < 0:
        raise RuntimeError("step budget exhausted")
    return out


SCAN_FIELDS = ["status", "t_end", "escaper", "log_growth", "min_sep", "n_syz", "ret_min", "ret_t", "steps"]


def scan(s0s, t_max, m=EQUAL, order=16, tol=1e-12, ret_tmin=1.0, r_esc=10.0, max_steps=400_000, threads=NCPU):
    """Integrate many initial conditions and return per-orbit diagnostics.
    Returns (table[n, 9], syzygy symbols[n, SYM_MAX]); columns are named in SCAN_FIELDS.
    """
    s0s = np.ascontiguousarray(s0s, dtype=float).reshape(-1, 18)
    m = np.ascontiguousarray(m, dtype=float)
    n = len(s0s)
    outs = np.empty((n, 9))
    syms = np.zeros((n, SYM_MAX), dtype=np.uint8)
    # interleave so slow regions are spread across threads
    chunks = [np.arange(k, n, threads * 8) for k in range(min(threads * 8, n))]

    def work(idx):
        a = np.ascontiguousarray(s0s[idx])
        o = np.empty((len(idx), 9))
        y = np.zeros((len(idx), SYM_MAX), dtype=np.uint8)
        _lib.tb_scan_batch(_p(a), len(idx), _p(m), float(t_max), order, tol, ret_tmin, r_esc, max_steps, _p(o),
                           y.ctypes.data_as(_up))
        outs[idx] = o
        syms[idx] = y

    with ThreadPoolExecutor(threads) as ex:
        list(ex.map(work, chunks))
    return outs, syms


def energy(s, m=EQUAL):
    s = np.asarray(s).reshape(-1, 18)
    x = s[:, :9].reshape(-1, 3, 3)
    v = s[:, 9:].reshape(-1, 3, 3)
    ke = 0.5 * np.einsum("b,nbd,nbd->n", m, v, v)
    pe = 0
    for i, j in ((0, 1), (0, 2), (1, 2)):
        pe = pe - m[i] * m[j] / np.linalg.norm(x[:, i] - x[:, j], axis=1)
    return ke + pe


def angular_momentum(s, m=EQUAL):
    s = np.asarray(s).reshape(-1, 18)
    x = s[:, :9].reshape(-1, 3, 3)
    v = s[:, 9:].reshape(-1, 3, 3)
    return np.einsum("b,nbd->nd", m, np.cross(x, v))


def collinear_state(p1, p2, p3=0.0):
    """Equal masses on a line at x = -1, +1, 0. The outer two share velocity
    (p1, p2, +/-p3); the middle one carries the opposite momentum.
    """
    return np.array([-1, 0, 0, 1, 0, 0, 0, 0, 0,
                     p1, p2, p3, p1, p2, -p3, -2 * p1, -2 * p2, 0], dtype=float)


_ip = ctypes.POINTER(ctypes.c_int)
_lib.tb_symscan_batch.restype = None
_lib.tb_symscan_batch.argtypes = [_dp, ctypes.c_int, _dp, ctypes.c_double, ctypes.c_int, ctypes.c_int,
                                  ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long,
                                  ctypes.c_int, _dp, _dp, _ip]


def symscan(s0s, t_max, sym_type, m=EQUAL, order=16, tol=1e-12, t_min=0.3, dthresh=0.3, r_esc=10.0,
            max_steps=200_000, kmax=6, threads=NCPU):
    """For each initial condition, the times at which the orbit comes closest to a
    symmetric state, i.e. a fixed point of the reversing symmetry (sym_apply in tb.c).
    Returns (times[n, kmax], distances[n, kmax], kept[n]).
    """
    s0s = np.ascontiguousarray(s0s, dtype=float).reshape(-1, 18)
    m = np.ascontiguousarray(m, dtype=float)
    n = len(s0s)
    tm = np.zeros((n, kmax))
    dm = np.full((n, kmax), np.inf)
    kp = np.zeros(n, dtype=np.int32)
    chunks = [np.arange(k, n, threads * 8) for k in range(min(threads * 8, n))]

    def work(idx):
        a = np.ascontiguousarray(s0s[idx])
        t_ = np.zeros((len(idx), kmax))
        d_ = np.full((len(idx), kmax), np.inf)
        k_ = np.zeros(len(idx), dtype=np.int32)
        _lib.tb_symscan_batch(_p(a), len(idx), _p(m), float(t_max), sym_type, order, tol, t_min, dthresh, r_esc,
                              max_steps, kmax, _p(t_), _p(d_), k_.ctypes.data_as(_ip))
        tm[idx], dm[idx], kp[idx] = t_, d_, k_

    with ThreadPoolExecutor(threads) as ex:
        list(ex.map(work, chunks))
    return tm, dm, kp


def sym_apply(s, sym_type):
    """Image of a state under the reversing symmetry (matches sym_apply in tb.c)."""
    s = np.asarray(s, dtype=float)
    x = s[:9].reshape(3, 3)[[1, 0, 2]]
    v = s[9:].reshape(3, 3)[[1, 0, 2]]
    sg = np.array([-1.0, -1.0 if sym_type == 0 else 1.0, -1.0])
    return np.concatenate([(x * sg).ravel(), (-v * sg).ravel()])


def isosceles_state(y, a, b, c):
    """Symmetric initial condition for the half-turn symmetry: an isosceles triangle
    in the xy-plane with its apex on the y axis. y = 0, b = 0 recovers the collinear planar slice.
    """
    return np.array([-1, y, 0, 1, y, 0, 0, -2 * y, 0,
                     a, b, c, a, -b, c, -2 * a, 0, -2 * c], dtype=float)
