"""Export research results for the web app (app/public/data)."""
import glob
import json
import os

import numpy as np
from PIL import Image

import tb

RES = os.path.join(tb.HERE, "results")
OUT = os.path.join(tb.HERE, "..", "app", "public", "data")
BODY_RGB = np.array([[0x6E, 0xC1, 0xFF], [0xFF, 0xC1, 0x4D], [0xFF, 0x5C, 0x6C]], dtype=float)
ORDERLY_BELOW = 10.0  # ln of shadow-orbit separation growth; the histogram has its minimum here

# Published values as recalled, used only to attach familiar names to orbits we
# found independently. A wrong entry here just leaves an orbit unnamed.
PUBLISHED = [
    ("Figure-eight", 0.347111, 0.532728, 6.324449),
    ("Butterfly I", 0.30689, 0.12551, 6.2356),
    ("Butterfly II", 0.39295, 0.09758, 7.0039),
    ("Bumblebee", 0.18428, 0.58719, 63.5345),
    ("Moth I", 0.46444, 0.39606, 14.8939),
    ("Moth II", 0.43917, 0.45297, 28.6703),
    ("Butterfly III", 0.40592, 0.23016, 13.8658),
    ("Moth III", 0.38344, 0.37736, 25.8406),
    ("Goggles", 0.08330, 0.12789, 10.4668),
    ("Butterfly IV", 0.350112, 0.07934, 79.4759),
    ("Dragonfly", 0.08058, 0.58884, 21.2710),
    ("Yarn", 0.55906, 0.34919, 55.5018),
    ("Yin-yang I", 0.51394, 0.30474, 17.3284),
    ("Yin-yang I", 0.28270, 0.32721, 10.9626),
    ("Yin-yang II", 0.41682, 0.33033, 55.7898),
    ("Yin-yang II", 0.41734, 0.31310, 54.2076),
]


def published_names(flat):
    """Map orbit index -> published name. An initial condition that coincides with a
    published one wins; otherwise the scale-invariant period must agree and the name be unclaimed.
    """
    names = {}
    for name, p1, p2, period in PUBLISHED:
        for i, o in enumerate(flat):
            if np.hypot(o["params"][0] - p1, o["params"][1] - p2) < 1e-3 and abs(o["period"] - period) < 0.01 * period:
                names[i] = name
    for name, p1, p2, period in PUBLISHED:
        if name in names.values():
            continue
        scaled = period * abs(3 * (p1 * p1 + p2 * p2) - 2.5) ** 1.5
        close = [(abs(o["scaled_period"] - scaled), i) for i, o in enumerate(flat) if i not in names]
        if close and min(close)[0] < 3e-4 * scaled:
            names[min(close)[1]] = name
    return names


def export_map():
    path = os.path.join(RES, "map.npz")
    if not os.path.exists(path):
        return None
    d = np.load(path)
    out, p_max, t_max = d["out"], float(d["p_max"]), float(d["t_max"])
    n = out.shape[0]
    status, t_end, esc = out[..., 0], out[..., 1], out[..., 2].astype(int)
    growth = np.nan_to_num(out[..., 3], nan=1e3, posinf=1e3)  # undefined after a near-collision: count as chaotic
    p = d["p"]
    P1, P2 = np.meshgrid(p, p)
    bound = 3 * (P1**2 + P2**2) < 2.5
    broke = status == 1
    orderly = (status == 0) & (growth < ORDERLY_BELOW)
    restless = (status != 1) & ~orderly

    # Saturated body colour where the system breaks up quickly, fading into a haze for
    # systems still bound at the end of the survey. The haze brightens with chaos strength.
    haze = np.array([52.0, 46, 104]) + np.clip(np.log(np.maximum(growth, 10) / 10) / np.log(60), 0, 1)[..., None] * [96, 92, 96]
    undecided = np.clip(t_end / t_max, 0, 1)[..., None] ** 0.6
    rgb = haze.copy()
    for b in range(3):
        m = broke & (esc == b)
        rgb[m] = BODY_RGB[b] * (1 - undecided[m]) + haze[m] * undecided[m]
    rgb[orderly] = [246, 240, 222]
    rgb[~bound] *= 0.3
    Image.fromarray(rgb[::-1].astype(np.uint8)).save(os.path.join(OUT, "map.png"))

    cells = np.stack([np.where(broke, esc, -1), t_end, growth], axis=-1).astype(np.float32)
    cells.tofile(os.path.join(OUT, "map.bin"))
    stats = {
        "n": n, "pMax": p_max, "tMax": t_max, "orderlyBelow": ORDERLY_BELOW,
        "boundStarts": int(bound.sum()),
        "brokeUp": float((broke & bound).sum() / bound.sum()),
        "orderly": float((orderly & bound).sum() / bound.sum()),
        "restless": float((restless & bound).sum() / bound.sum()),
        "escaperShare": [float((broke & bound & (esc == b)).sum() / max(1, (broke & bound).sum())) for b in range(3)],
    }
    json.dump(stats, open(os.path.join(OUT, "map.json"), "w"), indent=1)
    return stats


def thumb(o):
    path = tb.sample(np.array(o["state"]), o["period"], 181)[:, :9].reshape(-1, 3, 3)
    flat = path.reshape(-1, 3)
    if o["sym"] == 0:
        axes = np.eye(3)[:2]
    else:
        axes = np.linalg.svd(flat - flat.mean(axis=0), full_matrices=False)[2][:2]
    xy = path @ axes.T
    xy /= np.abs(xy).max()
    return [[f"{x:.3f},{-y:.3f}" for x, y in xy[:, b]] for b in range(3)]


def export_orbits():
    orbits = []
    flat_path = os.path.join(RES, "planar_orbits.json")
    flat = json.load(open(flat_path)) if os.path.exists(flat_path) else []
    spatial = []
    for f in sorted(glob.glob(os.path.join(RES, "spatial_orbits_seed*.json"))):
        spatial += json.load(open(f))
    family_file = os.path.join(RES, "exp8_family.json")
    on_family = {round(m["period"], 4): m["turn"] for m in json.load(open(family_file))["on_curve"]} if os.path.exists(family_file) else {}
    built_file = os.path.join(RES, "family_orbits.json")
    built = json.load(open(built_file)) if os.path.exists(built_file) else []
    on_family.update({round(o["period"], 4): o["turn"] for o in built if o["already_in_search"]})
    built = [o for o in built if not o["already_in_search"]]
    checked_file = os.path.join(RES, "exp9_catalogue.json")
    checked = json.load(open(checked_file)) if os.path.exists(checked_file) else None
    in_catalogue = {round(m["period"], 4): m["published_index"] for m in checked["matched"]} if checked else {}
    proofs_file = os.path.join(RES, "proofs.json")
    proved = {p["label"] for p in json.load(open(proofs_file)) if p["proved"]} if os.path.exists(proofs_file) else set()
    seen = []
    for kind, group in (("flat", flat), ("spatial", spatial + built)):
        group = sorted(group, key=lambda o: o["scaled_period"])
        names = published_names(group) if kind == "flat" else {}
        count = 0
        for index, o in enumerate(group):
            key = (kind, round(o["scaled_period"], 4), round(o["scaled_L"], 4))
            if key in seen:
                continue
            seen.append(key)
            count += 1
            known = names.get(index)
            stable = o["growth_per_period"] < 1.02
            turn = o.get("turn") or on_family.get(round(o["period"], 4))
            if kind == "spatial":
                known = in_catalogue.get(round(o["period"], 4))
            if kind == "spatial":
                tag = "computed from the family" if "turn" in o else "found by search" if turn else ""
                if checked:
                    tag = f"{tag}, {f'catalogue orbit O{known}' if known else 'not in the 2025 catalogue'}".lstrip(", ")
                name = f"Turn {turn.replace('0/1', '0')}" if turn else f"3D orbit {count}"
                labels = (f"family turn {turn}", f"spatial orbit, period {o['period']:.4f},")  # a family member may have been proved under either label
                if any(p.startswith(labels) for p in proved):
                    tag += ", existence proved"
            else:
                tag = "known orbit" if known else ""
                name = known or f"Planar orbit {count}"
            orbits.append({
                "id": f"{kind}-{count}",
                "kind": "family" if turn else kind,
                "turn": turn,
                "name": name,
                "tag": tag,
                "params": o["params"],
                "state": o["state"],
                "period": o["period"],
                "growth": o["growth_per_period"],
                "stable": stable,
                "word": o["word"],
                "closure": o["closure_error"],
                "thumb": thumb(o),
                "featured": known == "Figure-eight",
            })
    formulas_path = os.path.join(RES, "formulas.json")
    if os.path.exists(formulas_path):
        for i, f in enumerate(json.load(open(formulas_path)), 1):
            s0 = tb.collinear_state(*f["params"])
            path = tb.sample(s0, 60.0, 361)[:, :9].reshape(-1, 3, 3)[:, :, :2]
            path = path / np.abs(path).max()
            orbits.append({
                "id": f"orderly-{i}",
                "kind": "orderly",
                "name": f"Quasi-periodic orbit {i}",
                "tag": f"predicted by a {len(f['rates'])}-term formula",
                "params": f["params"],
                "state": [float(v) for v in s0],
                "period": None,
                "stable": True,
                "growth": 1.0,
                "word": "",
                "thumb": [[f"{x:.3f},{-y:.3f}" for x, y in path[:, b]] for b in range(3)],
                "formula": {"freqs": f["rates"], "amps": f["amps"]},
                "formulaStats": f["stats"],
            })
    from fractions import Fraction
    orbits.sort(key=lambda o: Fraction(o["turn"]) if o.get("turn") else 0)  # stable: only reorders the family
    json.dump(orbits, open(os.path.join(OUT, "orbits.json"), "w"))
    return orbits


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print(export_map())
    orbits = export_orbits()
    print(len(orbits), "orbits exported:", {k: sum(o["kind"] == k for o in orbits) for k in ("flat", "family", "spatial", "orderly")})
