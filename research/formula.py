"""Frequency analysis of regular orbits.

If the motion is quasi-periodic, each body's position (as x + iy) is
z(t) = sum_k a_k exp(i w_k t). The frequencies w_k are extracted from a span of
the integrated orbit one at a time, strongest first (in the spirit of Laskar's
frequency analysis), and the amplitudes a_k are then solved for jointly.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def _peak(t, resid, window):
    n = len(t)
    dt = t[1] - t[0]
    pad = 4 * n
    spec = np.abs(np.fft.fft(resid * window[:, None], pad, axis=0)) ** 2
    k = int(np.argmax(spec.sum(axis=1)))
    w0 = 2 * np.pi * np.fft.fftfreq(pad, dt)[k]
    bin_w = 2 * np.pi / (pad * dt)

    def neg_power(w):
        e = np.exp(-1j * w * t) * window
        return -np.sum(np.abs(e @ resid) ** 2)

    return minimize_scalar(neg_power, bounds=(w0 - bin_w, w0 + bin_w), method="bounded",
                           options={"xatol": 1e-13}).x


def fit(t, z, n_terms, resolve_every=6):
    """t: evenly spaced times; z: complex positions, shape (len(t), bodies).
    Returns (frequencies[K], amplitudes[K, bodies], rms residual within the fitted span).
    """
    tc = t - t.mean()
    window = np.hanning(len(t)) ** 2
    rates = []
    resid = z.copy()
    amps = None
    for k in range(n_terms):
        w = _peak(tc, resid, window)
        rates.append(w)
        if (k + 1) % resolve_every and k + 1 < n_terms:
            # cheap update: subtract this term; all amplitudes are re-solved jointly every few terms
            e = np.exp(1j * w * tc)
            resid = resid - np.outer(e, (np.conj(e) * window) @ resid / window.sum())
        else:
            basis = np.exp(1j * np.outer(tc, rates))
            amps = np.linalg.lstsq(basis, z, rcond=None)[0]
            resid = z - basis @ amps
    rates = np.array(rates)
    amps = amps * np.exp(-1j * rates * t.mean())[:, None]  # re-reference phases to t = 0
    return rates, amps, float(np.sqrt(np.mean(np.abs(resid) ** 2)))


def evaluate(rates, amps, t):
    return np.exp(1j * np.outer(t, rates)) @ amps
