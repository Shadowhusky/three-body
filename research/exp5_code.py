"""Experiment 5: symbolic dynamics.

Each time the three bodies become collinear (a syzygy), record which body is in
the middle (1, 2 or 3). The map stores the first 96 symbols of this syzygy
sequence for each initial condition. Two questions:

  1. Does the sequence act as an address? If neighbouring initial conditions
     share only a finite prefix, and halving their distance adds a fixed number
     of shared symbols, then each symbol localises the initial condition by a
     fixed factor.
  2. Is there a grammar? Count which words of length n occur at all.
"""
import json
import os

import numpy as np

import tb

d = np.load(os.path.join(tb.HERE, "results", "map.npz"))
out, sym = d["out"], d["sym"]
n = out.shape[0]
status, growth, n_syz = out[..., 0], out[..., 3], out[..., 5].astype(int)
orderly = (status == 0) & (growth < 10)
length = np.minimum(n_syz, tb.SYM_MAX)


def shared_prefix(step):
    """Mean number of leading symbols shared by initial conditions `step` pixels apart,
    where both sequences are long enough to tell; regular and non-regular orbits separately.
    """
    a, b = sym[:, :-step], sym[:, step:]
    usable = np.minimum(length[:, :-step], length[:, step:])
    same = np.cumprod(a == b, axis=2)
    same &= np.arange(tb.SYM_MAX)[None, None, :] < usable[..., None]
    shared = same.sum(axis=2)
    capped = shared >= usable  # sequences ended before they differed
    both_orderly = orderly[:, :-step] & orderly[:, step:]
    wild = ~orderly[:, :-step] & ~orderly[:, step:] & (usable >= 12)
    return {
        "orderly_mean": float(shared[both_orderly].mean()) if both_orderly.any() else None,
        "orderly_never_differ": float(capped[both_orderly].mean()) if both_orderly.any() else None,
        "wild_mean": float(shared[wild & ~capped].mean()),
        "wild_never_differ": float(capped[wild].mean()),
    }


res = {"pixel": float(d["p_max"]) / n, "shared_prefix": {str(s): shared_prefix(s) for s in (1, 2, 4, 8, 16)}}
m1, m4 = res["shared_prefix"]["1"]["wild_mean"], res["shared_prefix"]["4"]["wild_mean"]
res["symbols_per_halving"] = (m1 - m4) / 2
res["narrowing_per_symbol"] = float(2 ** (2 / (m1 - m4))) if m1 > m4 else None

# grammar: which words occur at all
rows = sym.reshape(-1, tb.SYM_MAX)
lens = length.ravel()
repeats = 0
pairs = 0
words = {}
for k in range(1, 13):
    seen = set()
    code = np.zeros(len(rows), dtype=np.int64)
    for start in range(0, tb.SYM_MAX - k + 1, 3):
        window = rows[:, start:start + k].astype(np.int64)
        code = (window * (4 ** np.arange(k))).sum(axis=1)
        seen.update(np.unique(code[lens >= start + k]).tolist())
    words[k] = len(seen)
valid = np.arange(tb.SYM_MAX - 1)[None, :] < (lens[:, None] - 1)
res["immediate_repeat_fraction"] = float(((rows[:, :-1] == rows[:, 1:]) & valid).sum() / valid.sum())
res["distinct_words"] = words
res["all_possible"] = {k: 3**k for k in words}
res["possible_without_repeats"] = {k: 3 * 2 ** (k - 1) for k in words}
ks = np.array(sorted(words))
res["entropy_per_symbol"] = float(np.polyfit(ks[4:], np.log([words[k] for k in ks[4:]]), 1)[0])

json.dump(res, open(os.path.join(tb.HERE, "results", "exp5_code.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
