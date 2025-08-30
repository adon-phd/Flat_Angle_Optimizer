
# ------------------------------------------------------------------------------
# Copyright (c) 2025 adon-phd
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to use,
# copy, modify, and merge the Software for non-commercial purposes only,
# including academic and personal projects, subject to the following conditions:
#
# 1. The above copyright notice and this permission notice shall be included in
#    all copies or substantial portions of the Software.
# 2. Commercial use of the Software, in whole or in part, is strictly prohibited
#    without prior written permission from the copyright holder.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ------------------------------------------------------------------------------

#!/usr/bin/env python3
"""
WBPP rotator-angle planner + optional header tagging (2dp, joint ILP).

Highlights
- Joint optimisation across all filters: one set of centres for every filter.
- Strict ±tolerance with 2-decimal precision (no padding or epsilon tricks).
- Integer ILP (via PuLP) in centi-degrees to avoid boundary wobble.
- Greedy weighted fallback (max frames per step) if PuLP is not available.
- Optional FLATS folder; you can also supply --existing-flats (global or per-filter).
- Concise console report; optional XLSX debug output with --debug-nearest.
- WBPP tagging (ROTGRP) with ASCII-safe comments.

Dependencies
- Required: astropy
- Optional: PuLP (for ILP; falls back to greedy if not installed)
- Optional: openpyxl (for --debug-nearest XLSX output)
"""

import argparse
import os
import sys
from collections import Counter, defaultdict
from math import fmod, cos, sin, atan2, radians, degrees
from typing import List, Tuple, Dict, Optional, Union

from astropy.io import fits

from typing import Iterable


def try_merge_centres_to_existing(
    centres: List[float],
    existing: List[float],
    tol: float,
    light_angles: List[float],
    merge_radius: float,
    allow_loss: int
) -> Tuple[List[float], Dict[Tuple[float, float], int]]:
    """
    Attempt to replace any chosen centre 'c' with a nearby existing flat 'e' (within merge_radius)
    if doing so strands at most 'allow_loss' frames that were uniquely covered by 'c'.

    - All math uses your 2-dp quantisation and strict ±tol coverage.
    - 'existing' should be the angles you want to prefer (e.g., --existing-flats or any scanned flats).

    Returns:
      (new_centres, losses) where losses[(c,e)] = frames lost by that specific merge decision.
    """
    if not centres or not existing or merge_radius <= 0.0:
        return sorted(set(q(c) for c in centres)), {}

    centres = sorted(set(q(c) for c in centres))
    existing = sorted(set(q(e) for e in existing))

    # Frame weights at each 2-dp light angle
    from collections import Counter
    counts = Counter(q(a) for a in light_angles)
    uniq = sorted(counts.keys())

    def coverset(cs: Iterable[float]) -> set:
        S = set()
        for c in cs:
            for a in uniq:
                if circ_dist(a, c) <= tol:
                    S.add(a)
        return S

    # Start from full coverage with current centres
    current = set(centres)
    current_cov = coverset(current)
    losses: Dict[Tuple[float, float], int] = {}

    # Consider each centre independently (stable and predictable)
    for c in list(current):
        # find nearest existing within merge_radius
        close = [(circ_dist(c, e), e) for e in existing if circ_dist(c, e) <= merge_radius]
        if not close:
            continue
        _, e_best = min(close)

        # Angles that would become exposed if we *remove* c (unique to c)
        cov_without_c = coverset([x for x in current if x != c])
        unique_angles = current_cov - cov_without_c  # angles only c was contributing to

        # If we replace c -> e_best, what angles remain covered?
        cov_with_e = coverset(list((current - {c}) | {e_best}))
        newly_uncovered_angles = unique_angles - cov_with_e
        lost_frames = sum(counts[a] for a in newly_uncovered_angles)
        losses[(c, e_best)] = lost_frames

        if lost_frames <= allow_loss:
            # accept the merge
            current.remove(c)
            current.add(e_best)
            current_cov = cov_with_e  # update baseline to reflect the accepted merge

    return sorted(q(x) for x in current), losses

# -------------------- 2dp quantisation & helpers --------------------

ANG_PREC = 2  # fixed: two decimals everywhere

def q(a: float) -> float:
    a = fmod(a, 360.0)
    if a < 0:
        a += 360.0
    return round(a, ANG_PREC)

def fmt_angle(a: float) -> str:
    return f"{q(a):.2f}°"

def circ_dist(a: float, b: float) -> float:
    a, b = q(a), q(b)
    d = abs(a - b)
    return d if d <= 180.0 else 360.0 - d

def angle_window(c: float, tol: float) -> Tuple[float, float]:
    return q(c - tol), q(c + tol)

def fmt_window(c: float, tol: float) -> str:
    lo, hi = angle_window(c, tol)
    if hi < lo:
        return f"{fmt_angle(lo)}–360.00° & 0.00°–{fmt_angle(hi)}"
    return f"{fmt_angle(lo)}–{fmt_angle(hi)}"

def ascii_safe(s: str) -> str:
    for a, b in (("±", "+/-"), ("°", " deg "), ("–", "-"), ("—", "-")):
        s = s.replace(a, b)
    return s.encode("ascii", "ignore").decode("ascii")

# -------------------- FITS helpers --------------------

DEFAULT_EXTS = (".fits", ".fit", ".fts", ".fits.fz", ".fz")
ANGLE_KEYS_DEFAULT = ["ROTANGLE", "ROTATOR", "ROTANG", "ROT_POS", "ROTATOR_ANGLE", "CROTA2", "PA", "ROT_PA"]
FILTER_KEYS_DEFAULT = ["FILTER", "FILTNAME", "INSFLNAM", "FWHEEL", "FILTNAM", "FILTERID"]
FRAME_TYPE_KEYS = ["IMAGETYP", "FRAME", "OBSTYPE", "IMAGETYPE"]

def find_files(root: str, exts: tuple) -> List[str]:
    out = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.lower().endswith(exts):
                out.append(os.path.join(dp, fn))
    return out

def read_header(path: str) -> Optional[fits.Header]:
    try:
        with fits.open(path, memmap=True) as hdul:
            return hdul[0].header
    except Exception:
        return None

def parse_float(s) -> Optional[float]:
    try:
        return float(s)
    except Exception:
        return None

def get_rotator_angle(h: fits.Header, angle_keys: List[str]) -> Optional[float]:
    for k in angle_keys:
        if k in h:
            v = parse_float(h.get(k))
            if v is not None:
                return q(v)  # quantise at ingest
    return None

def get_filter_name(h: fits.Header, filter_keys: List[str]) -> Optional[str]:
    for k in filter_keys:
        if k in h:
            val = str(h.get(k)).strip()
            if val:
                return val.upper()
    return None

def get_frame_type(h: fits.Header) -> Optional[str]:
    for k in FRAME_TYPE_KEYS:
        if k in h:
            val = str(h.get(k)).strip().upper()
            if "LIGHT" in val or val in {"L", "LIGHTS"}:
                return "LIGHT"
            if "FLAT" in val or val in {"F", "FLATS"}:
                return "FLAT"
    return None

# -------------------- Parse --existing-flats --------------------

def parse_existing_flats(spec: Union[None, str, list]) -> Tuple[List[float], Dict[str, List[float]]]:
    """
    Accepts: None, string, or list of tokens (nargs='+')
    Tokens: "300", "315.26", "300@HA"
    Returns: (global_angles, per_filter_angles)
    """
    globals_, perf = [], defaultdict(list)
    if not spec:
        return [], {}
    if isinstance(spec, list):
        spec = " ".join(str(t) for t in spec)
    for tok in spec.replace(",", " ").split():
        if "@" in tok:
            a, f = tok.split("@", 1)
            try:
                perf[f.strip().upper()].append(q(float(a)))
            except Exception:
                pass
        else:
            try:
                globals_.append(q(float(tok)))
            except Exception:
                pass
    globals_ = sorted(set(globals_))
    for k in list(perf.keys()):
        perf[k] = sorted(set(perf[k]))
    return globals_, perf

# -------------------- Weighted greedy fallback --------------------

def greedy_cover_weighted(light_angles: List[float], tol: float) -> List[float]:
    """
    Weighted greedy cover on the circle (2dp).
    Each step picks the centre that covers the most *frames* in a 2*tol window.
    """
    if not light_angles:
        return []

    cnt = Counter(q(a) for a in light_angles)  # weights per 2dp angle
    base = sorted(cnt.keys())
    w    = [cnt[a] for a in base]
    n    = len(base)

    ext_a = base + [q(a + 360.0) for a in base]
    ext_w = w    + w

    covered = [False] * n
    remaining = sum(w)
    centres: List[float] = []

    while remaining > 0:
        best_i, best_gain, best_j = None, -1, None
        for i in range(n):
            if covered[i]:
                continue
            centre = q(ext_a[i] + tol)
            j = i
            gain = 0
            while j < i + n and circ_dist(ext_a[j] % 360.0, centre) <= tol:
                idx = j % n
                if not covered[idx]:
                    gain += ext_w[j]
                j += 1
            if gain > best_gain:
                best_i, best_gain, best_j = i, gain, j

        if best_i is None or best_gain <= 0:
            break

        centre = q(ext_a[best_i] + tol)
        for j in range(best_i, best_j):
            idx = j % n
            if not covered[idx]:
                remaining -= w[idx]
                covered[idx] = True
        centres.append(centre)

    return sorted(set(q(c) for c in centres))


def prune_redundant(centres: List[float], uniq_angles: List[float], tol_deg: float) -> List[float]:
    """Remove centres that do not uniquely cover any angle."""
    centres = sorted(set(centres))
    changed = True
    while changed:
        changed = False
        cover_map = {c: set(a for a in uniq_angles if circ_dist(a, c) <= tol_deg)
                     for c in centres}
        for c in list(centres):
            # union of others without c
            others = set().union(*(cover_map[o] for o in centres if o != c))
            if cover_map[c] <= others:
                # c adds nothing unique, drop it
                centres.remove(c)
                changed = True
                break
    return centres

# -------------------- Integer helpers (centi-degrees) --------------------

def to_cd(a: float) -> int:
    return int(round(q(a) * 100))

def from_cd(x: int) -> float:
    return q(x / 100.0)

def dist_cd(a_cd: int, b_cd: int) -> int:
    a_cd %= 36000
    b_cd %= 36000
    d = abs(a_cd - b_cd)
    return d if d <= 18000 else 36000 - d

# -------------------- Joint optimiser (across all filters) --------------------
"""
Chooses a common set of rotator centres that work across all filters.
Two modes:
- cover_all:   Make sure every light angle is covered, using the smallest
               number of centres. Reuse existing flats if possible.
- max_cover:   Given a budget of new flats, cover as many frames as 
               possible (specialized use-case).
Returns a tuple:
(chosen_centres, coverage_count_per_centre, used_ilp_flag)
"""

def optimise_joint(
    light_records: List[Tuple[str, float, str]],   # (path, angle, filter)
    existing_global: List[float],
    tol_deg: float,
    mode: str = "cover_all",
    budget: Optional[int] = None,
) -> Tuple[List[float], Dict[float, int], bool]:
    lights = [q(a) for _, a, _ in light_records]
    if not lights:
        return [], {}, True
    counts = Counter(lights)
    uniq_angles = sorted(counts.keys())

    cand = sorted(set(existing_global + uniq_angles))
    tol_cd = int(round(tol_deg * 100))
    cand_cd = [to_cd(c) for c in cand]
    existing_set = set(to_cd(x) for x in existing_global)
    uniq_cd = [to_cd(a) for a in uniq_angles]

    try:
        import pulp
    except Exception:
        # Greedy fallback
        centres = greedy_cover_weighted(lights, tol_deg)
        centres = prune_redundant(centres, uniq_angles, tol_deg)
        final = sorted(set(centres) | set(existing_global))
        cov = {ci: sum(counts[a] for a in uniq_angles if circ_dist(a, ci) <= tol_deg)
               for ci in final}
        return final, cov, False

    prob = pulp.LpProblem("joint_rotator",
                          pulp.LpMinimize if mode == "cover_all" else pulp.LpMaximize)
    x = {j: pulp.LpVariable(f"x_{j}", 0, 1, cat="Binary") for j in range(len(cand_cd))}

    if mode == "cover_all":
        # Objective: minimise new centres, then total centres
        costs = [0 if cand_cd[j] in existing_set else 1 for j in range(len(cand_cd))]
        prob += pulp.lpSum(costs[j] * x[j] for j in range(len(cand_cd))) \
              + 1e-3 * pulp.lpSum(x[j] for j in range(len(cand_cd)))

        # Coverage constraint
        for i, ai in enumerate(uniq_cd):
            J = [j for j, cj in enumerate(cand_cd) if dist_cd(ai, cj) <= tol_cd]
            if J:
                prob += pulp.lpSum(x[j] for j in J) >= 1
            else:
                raise RuntimeError(f"No candidate covers angle {from_cd(ai)} at ±{tol_deg}°")
    else:
        if budget is None:
            raise ValueError("--budget is required for --optimize max_cover")
        w = [counts[uniq_angles[i]] for i in range(len(uniq_angles))]
        prob += pulp.lpSum(w[i] * (pulp.lpSum(
            x[j] for j, cj in enumerate(cand_cd)
            if dist_cd(uniq_cd[i], cj) <= tol_cd)) >= 1 for i in range(len(uniq_cd)))
        prob += pulp.lpSum((0 if cand_cd[j] in existing_set else 1) * x[j]
                           for j in range(len(cand_cd))) <= int(budget)

    solver = pulp.PULP_CBC_CMD(msg=False)
    res = prob.solve(solver)
    if pulp.LpStatus[res] not in {"Optimal", "Feasible"}:
        centres = greedy_cover_weighted(lights, tol_deg)
        centres = prune_redundant(centres, uniq_angles, tol_deg)
        final = sorted(set(centres) | set(existing_global))
        cov = {ci: sum(counts[a] for a in uniq_angles if circ_dist(a, ci) <= tol_deg)
               for ci in final}
        return final, cov, False

    chosen = [from_cd(cand_cd[j]) for j in range(len(cand_cd)) if pulp.value(x[j]) >= 0.5]
    chosen = sorted(set(q(c) for c in chosen))
    cov = {c: sum(counts[a] for a in uniq_angles if circ_dist(a, c) <= tol_deg)
           for c in chosen}
    return chosen, cov, True

# -------------------- Per-filter needs (at chosen centres) --------------------
    """
    For each chosen centre, count per-filter LIGHTS that are NOT covered
    by any flat of the same filter within ±tol.
    """
def per_filter_needs(
    centres: List[float],
    tol: float,
    light_records: List[Tuple[str, float, str]],
    flat_records: List[Tuple[str, float, str]]
) -> Dict[float, Dict[str, int]]:
    flats_by_filter = defaultdict(list)
    for _, a_f, f_f in flat_records:
        flats_by_filter[f_f].append(q(a_f))

    needed: Dict[float, Dict[str, int]] = {}
    for c in centres:
        per_f = defaultdict(int)
        for _, a_l, f_l in light_records:
            if circ_dist(a_l, c) <= tol:
                covered = any(circ_dist(a_l, a_f) <= tol for a_f in flats_by_filter.get(f_l, []))
                if not covered:
                    per_f[f_l] += 1
        if per_f:
            needed[q(c)] = dict(sorted(per_f.items()))
    return needed

# -------------------- CSV / XLSX --------------------

def write_csv(path: str, rows: List[Dict[str, str]], per_tol_centres: Dict[float, List[float]]):
    import csv
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if rows:
            w.writerow(list(rows[0].keys()))
            for r in rows:
                w.writerow(list(r.values()))
        w.writerow([])
        w.writerow(["centres_by_tolerance"])
        for t, centres in per_tol_centres.items():
            w.writerow([f"±{t}°", ", ".join(f"{q(c):.2f}" for c in centres)])

def write_xlsx_with_debug(xlsx_path: str, summary_rows: List[Dict[str, str]], debug_rows: List[Dict]):
    try:
        from openpyxl import Workbook
    except Exception:
        # Fallback CSV
        dbg_csv = os.path.splitext(xlsx_path)[0] + "_debug.csv"
        import csv
        os.makedirs(os.path.dirname(os.path.abspath(dbg_csv)), exist_ok=True)
        with open(dbg_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["tolerance_deg", "centre_deg", "centre_window",
                        "filter", "light_angle_deg", "delta_light_to_centre_deg",
                        "nearest_flat_angle_deg", "delta_light_to_nearest_flat_deg",
                        "delta_centre_to_nearest_flat_deg"])
            for r in debug_rows:
                w.writerow([f"{r['tol']:.2f}", f"{q(r['centre']):.2f}", r['window'],
                            r['filter'], f"{q(r['light_angle']):.2f}", f"{r['d_light_centre']:.2f}",
                            (f"{q(r['nearest_flat']):.2f}" if r['nearest_flat'] is not None else ""),
                            (f"{r['d_light_flat']:.2f}" if r['d_light_flat'] is not None else ""),
                            (f"{r['d_centre_flat']:.2f}" if r['d_centre_flat'] is not None else "")])
        print(f"(openpyxl not found) Wrote debug rows to: {dbg_csv}")
        return

    wb = Workbook()
    ws = wb.active
    ws.title = "Summary"
    if summary_rows:
        ws.append(list(summary_rows[0].keys()))
        for r in summary_rows:
            row = []
            for v in r.values():
                try:
                    vv = float(v)
                    row.append(vv)
                except Exception:
                    row.append(v)
            ws.append(row)

    ws_dbg = wb.create_sheet("DebugNearest")
    ws_dbg.append(["tolerance_deg", "centre_deg", "centre_window",
                   "filter", "light_angle_deg", "delta_light_to_centre_deg",
                   "nearest_flat_angle_deg", "delta_light_to_nearest_flat_deg",
                   "delta_centre_flat_deg"])
    for r in debug_rows:
        ws_dbg.append([
            float(f"{r['tol']:.2f}"),
            float(f"{q(r['centre']):.2f}"),
            r['window'],
            r['filter'],
            float(f"{q(r['light_angle']):.2f}"),
            float(f"{r['d_light_centre']:.2f}"),
            (float(f"{q(r['nearest_flat']):.2f}") if r['nearest_flat'] is not None else None),
            (float(f"{r['d_light_flat']:.2f}") if r['d_light_flat'] is not None else None),
            (float(f"{r['d_centre_flat']:.2f}") if r['d_centre_flat'] is not None else None),
        ])
    os.makedirs(os.path.dirname(os.path.abspath(xlsx_path)), exist_ok=True)
    wb.save(xlsx_path)
    print(f"XLSX written: {xlsx_path}")

# -------------------- Debug collection (optional) --------------------

def collect_debug_rows(
    centres: List[float],
    tol_deg: float,
    light_records: List[Tuple[str, float, str]],
    flat_records: List[Tuple[str, float, str]]
) -> List[Dict]:
    rows = []
    flats_by_filter = defaultdict(list)
    for _, ang, filt in flat_records:
        flats_by_filter[filt].append(q(ang))
    for c in centres:
        window = fmt_window(c, tol_deg)
        lights_near = defaultdict(list)
        for _, ang, filt in light_records:
            da = circ_dist(ang, c)
            if da <= tol_deg:
                lights_near[filt].append((q(ang), da))
        for filt, pairs in lights_near.items():
            f_list = flats_by_filter.get(filt, [])
            nf, dcf = (None, None)
            if f_list:
                dcf, nf = min((circ_dist(c, fa), fa) for fa in f_list)
            for a, d_c in sorted(pairs):
                if f_list:
                    dlf, nfa = min((circ_dist(a, fa), fa) for fa in f_list)
                else:
                    dlf, nfa = (None, None)
                rows.append({
                    "tol": tol_deg,
                    "centre": q(c),
                    "window": window,
                    "filter": filt,
                    "light_angle": q(a),
                    "d_light_centre": d_c,
                    "nearest_flat": nfa,
                    "d_light_flat": dlf,
                    "d_centre_flat": dcf
                })
    return rows

# -------------------- Tagging helpers --------------------

# def snap_centres_to_existing_flats(centres: List[float], flats: List[float], tol: float) -> List[float]:
#     if not centres or not flats:
#         return centres[:]
#     snapped = []
#     flats = [q(f) for f in flats]
#     for c in centres:
#         d, best = min((circ_dist(c, fa), fa) for fa in flats)
#         snapped.append(best if d <= tol else q(c))
#     return sorted(set(q(a) for a in snapped))

def snap_centres_to_existing_flats(centres: List[float], flats: List[float], tol: float) -> List[float]:
    if not centres or not flats:
        return centres[:]
    snapped = []
    flats_q = [q(f) for f in flats]
    for c in centres:
        d, best = min((circ_dist(c, fa), fa) for fa in flats_q)
        snapped.append(best if d <= tol else q(c))
    # Always include all real flats
    snapped.extend(flats_q)
    return sorted(set(q(a) for a in snapped))


# def assign_to_nearest(angle: float, centres: List[float], tol: float) -> Optional[float]:
#     if not centres:
#         return None
#     best_c, best_d = None, 1e9
#     for c in centres:
#         d = circ_dist(angle, c)
#         if d < best_d:
#             best_d, best_c = d, c
#     return best_c if best_d <= tol else None

def assign_to_nearest(angle: float,
                      centres: List[float],
                      tol: float,
                      merge_r: float = 0.0) -> Optional[float]:
    """Return nearest centre if within tolerance (expanded by merge_r)."""
    if not centres:
        return None
    best_c, best_d = None, 1e9
    for c in centres:
        d = circ_dist(angle, c)
        if d < best_d:
            best_d, best_c = d, c
    # add merge_r to tolerance so tagging matches analysis
    return best_c if best_d <= tol + merge_r + 1e-6 else None

def write_rotgrp_keyword(path: str, label: str, comment: str,
                         keyword: str, overwrite: bool) -> bool:
    try:
        with fits.open(path, mode="update") as hdul:
            hdr = hdul[0].header
            if (keyword in hdr) and (not overwrite):
                return True
            hdr[keyword] = (label, ascii_safe(comment))
            hdul.flush()
        return True
    except Exception as e:
        print(f"Failed to write {keyword} for {path}: {e}", file=sys.stderr)
        return False

# -------------------- GUI --------------------

def script_dir() -> str:
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.getcwd()

def pick_directory_gui(title: str = "Select directory", initialdir: Optional[str] = None) -> Optional[str]:
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:
        return None
    try:
        root = tk.Tk()
        root.withdraw()
        root.update_idletasks()
        kwargs = {"title": title}
        if initialdir and os.path.isdir(initialdir):
            kwargs["initialdir"] = initialdir
        path = filedialog.askdirectory(**kwargs)
        root.destroy()
        return path if path else None
    except Exception:
        return None

# -------------------- Main --------------------

def main():
    default_csv = os.path.join(script_dir(), "rotator_flat_report.csv")
    default_xlsx = os.path.join(script_dir(), "rotator_flat_report.xlsx")

    p = argparse.ArgumentParser(description="WBPP rotator planner + tagging (2dp, joint ILP).")
    p.add_argument("--lights", help="Root directory of LIGHT frames")
    p.add_argument("--flats", help="Root directory of FLAT frames (optional)")
    p.add_argument("--gui", action="store_true", help="Open GUI pickers if paths not provided")
    p.add_argument("--picker-default", help="Default initial folder for GUI pickers")

    p.add_argument("--extensions", nargs="+", default=list(DEFAULT_EXTS))
    p.add_argument("--angle-keys", nargs="+", default=ANGLE_KEYS_DEFAULT)
    p.add_argument("--filter-keys", nargs="+", default=FILTER_KEYS_DEFAULT)
    p.add_argument("--tolerances", nargs="+", type=float, default=[1, 2, 5, 10, 15])

    p.add_argument("--existing-flats", "--existing_flats", nargs="+", metavar="ANGLE[@FILTER]",
                   help="Angles available without a folder; e.g. 300,315.26 or 300@HA 300@OIII")

    # Optimiser mode — kept for parity with earlier versions
    p.add_argument("--optimize", choices=["cover_all", "max_cover"], default="cover_all",
                   help="Optimise to cover all lights or maximise coverage under a budget")
    p.add_argument("--budget", type=int, help="Budget of NEW centres (required for max_cover)")

    # Display / debug
    p.add_argument("--csv", default=default_csv, help="CSV path (summary)")
    p.add_argument("--xlsx", default=default_xlsx, help="XLSX path (when --debug-nearest)")
    p.add_argument("--debug-nearest", action="store_true", help="Console debug + XLSX DebugNearest sheet")

    # Frame-type fallback
    p.add_argument("--assume-lights", action="store_true", help="Treat files under --lights as LIGHT if missing type")
    p.add_argument("--assume-flats", action="store_true", help="Treat files under --flats as FLAT if missing type")

    # Tagging
    p.add_argument("--noninteractive", action="store_true", help="Skip prompt and tag using --tolerance")
    p.add_argument("--tolerance", type=float, help="Tolerance for tagging (deg)")
    p.add_argument("--keyword", default="ROTGRP", help="Header keyword to write")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing keyword")
    p.add_argument("--dry-run", action="store_true", help="Compute only; do not modify files")


    p.add_argument(
        "--centre-merge-radius",
        type=float,
        default=0.0,
        help="After optimisation, collapse centres to a nearby EXISTING flat within this radius (deg). 0 = off.",
    )
    p.add_argument(
        "--merge-allow-loss",
        type=int,
        default=0,
        help="Max frames allowed to become uncovered by a single centre merge (default 0).",
    )

    args = p.parse_args()

    if args.picker_default:
        args.picker_default = os.path.expanduser(os.path.expandvars(args.picker_default))

    # Parse existing flats
    existing_global, existing_per_filter = parse_existing_flats(args.existing_flats)
    if existing_global or existing_per_filter:
        pretty_pf = ", ".join(f"{f}:{', '.join(f'{a:.2f}' for a in angs)}" for f, angs in existing_per_filter.items())
        print(f"Existing flat angles (global): {', '.join(f'{a:.2f}' for a in existing_global) or '(none)'}")
        print(f"Existing flat angles (per-filter): {pretty_pf or '(none)'}")

    # Resolve roots
    lights_root = args.lights
    flats_root  = args.flats

    if args.gui or not lights_root:
        print("Select LIGHTS root…")
        picked = pick_directory_gui("Select LIGHTS root directory", args.picker_default)
        if picked:
            lights_root = picked
            print(f"LIGHTS root: {lights_root}")
        elif not args.lights:
            print("No LIGHTS directory selected; exiting.", file=sys.stderr)
            sys.exit(1)

    have_virtual_flats = bool(existing_global or existing_per_filter)
    if args.gui and not flats_root and not have_virtual_flats:
        print("Select FLATS root (optional)…")
        picked = pick_directory_gui("Select FLATS root directory (optional)", args.picker_default)
        if picked:
            flats_root = picked
            print(f"FLATS root:  {flats_root}")
        else:
            print("No FLATS directory selected; proceeding without flats.")

    # Gather files
    exts = tuple(s.lower() for s in args.extensions)
    lights_files = find_files(lights_root, exts) if lights_root else []
    flats_files  = find_files(flats_root,  exts) if flats_root else []

    if not lights_files:
        print("No matching LIGHT files found.", file=sys.stderr)
        sys.exit(1)

    angle_keys  = [k.upper() for k in args.angle_keys]
    filter_keys = [k.upper() for k in args.filter_keys]

    headers_ok = with_angle = 0
    lights: List[float] = []
    flats_all: List[float] = []
    light_records: List[Tuple[str, float, str]] = []
    flat_records:  List[Tuple[str, float, str]] = []
    file_records:  List[Tuple[str, float, str, str]] = []  # (path, angle, type, filter)

    def process(files, assume_type: Optional[str]):
        nonlocal headers_ok, with_angle, lights, flats_all, light_records, flat_records, file_records
        for fp in files:
            h = read_header(fp)
            if h is None:
                continue
            headers_ok += 1
            ang = get_rotator_angle(h, angle_keys)
            if ang is None:
                continue
            filt = get_filter_name(h, filter_keys) or "UNKNOWN"
            with_angle += 1
            ftype = get_frame_type(h)
            if ftype is None and assume_type:
                ftype = assume_type
            if ftype == "LIGHT":
                lights.append(ang)
                light_records.append((fp, ang, filt))
            elif ftype == "FLAT":
                flats_all.append(ang)
                flat_records.append((fp, ang, filt))
            file_records.append((fp, ang, ftype or "", filt))

    process(lights_files, "LIGHT" if args.assume_lights else None)
    process(flats_files,  "FLAT"  if args.assume_flats  else None)

    # Inject virtual flats from --existing-flats (global and per-filter)
    for ang in existing_global:
        flats_all.append(ang)
        for _, _, filt in light_records:
            flat_records.append(("<virtual>", ang, filt))
    for filt, angs in existing_per_filter.items():
        for ang in angs:
            flats_all.append(ang)
            flat_records.append(("<virtual>", ang, filt))

    if not light_records:
        print("No frames with a valid rotator angle were found.", file=sys.stderr)
        sys.exit(1)

    # Summary
    light_unique = sorted(set(lights))
    flat_unique  = sorted(set(flats_all))
    print("\n=== Scan Summary ===")
    print(f"Files scanned (lights):     {len(lights_files)}")
    print(f"Files scanned (flats):      {len(flats_files)}")
    print(f"Readable FITS headers:      {headers_ok}")
    print(f"Frames with angle:          {with_angle}")
    print(f"LIGHT frames:               {len(lights)}")
    print(f"FLAT frames:                {len(flats_all)}")
    print(f"Unique LIGHT angles:        {len(light_unique)}")
    print(f"Unique FLAT angles:         {len(flat_unique)}")

    # ---- Optimise jointly per tolerance ----
    per_tol_centres: Dict[float, List[float]] = {}
    summary_rows: List[Dict[str, str]] = []
    debug_rows_all: List[Dict] = []

    for t in args.tolerances:
        centres, cov_map, used_ilp = optimise_joint(
            light_records=light_records,
            existing_global=existing_global,
            tol_deg=t,
            mode=args.optimize,
            budget=args.budget
        )
        per_tol_centres[t] = centres


        # Optional fuzzy collapse towards existing flats (uses provided --existing-flats first;
        # if none were provided, fall back to any scanned flats)
        existing_for_merge = sorted(set(existing_global)) or sorted(set(f for f in flats_all))
        merge_r = float(args.centre_merge_radius or 0.0)
        allow_loss = int(args.merge_allow_loss or 0)

        if existing_for_merge and merge_r > 0.0:
            merged, losses = try_merge_centres_to_existing(
                centres=centres,
                existing=existing_for_merge,
                tol=t,
                light_angles=lights,   # all frames for frame-weighted loss
                merge_radius=merge_r,
                allow_loss=allow_loss
            )

            if merged != centres:
                print(f"Fuzzy merge: collapsed {len(centres)} → {len(merged)} centres "
                    f"(radius {merge_r:.2f}°, allow_loss {allow_loss}).")
                for (c_old, e_new), lost in losses.items():
                    if c_old in centres and e_new in merged and c_old not in merged:
                        print(f"  merged {c_old:.2f}° → {e_new:.2f}° (lost {lost} frames)")
                centres = merged

            # Preserve real flat centres if they exist within the dataset
            if flats_all:
                for f in flats_all:
                    f_q = q(f)
                    # Is there at least one light within tolerance of this flat?
                    if any(circ_dist(l, f_q) <= t for l in lights):
                        # If optimiser dropped it, restore it
                        if all(circ_dist(f_q, c) > 1e-6 for c in centres):
                            centres.append(f_q)
                centres = sorted(set(q(c) for c in centres))

        # Coverage accounting (union across all chosen centres)
        counts_all = Counter(lights)
        covered_angles = set()
        for c in centres:
            for a in light_unique:
                if circ_dist(a, c) <= t:
                    covered_angles.add(a)
        covered_frames = sum(counts_all[a] for a in covered_angles)
        uncovered_angles = [a for a in light_unique if a not in covered_angles]
        uncovered_frames = sum(counts_all[a] for a in uncovered_angles)

        print(f"\n--- Tolerance ±{t}° ---")
        print("Chosen centres (optimised): " + (", ".join(f"{q(c):.2f}" for c in centres) if centres else "(none)"))
        print(f"Covered unique angles:      {len(covered_angles)}")
        print(f"Covered light frames:       {covered_frames}")
        print(f"UNCOVERED unique angles:    {len(uncovered_angles)}")
        print(f"UNCOVERED light frames:     {uncovered_frames}")
        if uncovered_angles:
            print("Residual light angles:")
            print("  " + ", ".join(f"{q(a):.2f}" for a in uncovered_angles))

        # Per-filter needs at these joint centres
        needs = per_filter_needs(centres, t, light_records, flat_records)
        if needs:
            print("\nPer-filter requirements at chosen centres:")
            for c in centres:
                if c not in needs:
                    continue
                win = fmt_window(c, t)
                items = ", ".join(f"{filt} (needs flats; {cnt} lights)" for filt, cnt in needs[c].items())
                print(f"  Angle {fmt_angle(c)} (covers {win}) → {items}")
        else:
            print("All required filters/angles are covered at these centres.")

        # Rows for CSV summary
        for c in centres:
            win = fmt_window(c, t)
            covered_here = sum(counts_all[a] for a in light_unique if circ_dist(a, c) <= t)
            summary_rows.append({
                "tolerance_deg": f"{t}",
                "centre_deg": f"{q(c):.2f}",
                "coverage_window": win,
                "covered_light_frames": f"{covered_here}",
            })

        # Debug rows (optional)
        if args.debug_nearest and centres:
            debug_rows_all.extend(collect_debug_rows(centres, t, light_records, flat_records))

    # Write CSV / XLSX
    try:
        write_csv(args.csv, summary_rows, per_tol_centres)
        print(f"\nCSV written: {args.csv}")
    except Exception as e:
        print(f"\nFailed to write CSV: {e}", file=sys.stderr)

    if args.debug_nearest and debug_rows_all:
        try:
            write_xlsx_with_debug(args.xlsx, summary_rows, debug_rows_all)
        except Exception as e:
            print(f"Failed to write XLSX: {e}", file=sys.stderr)

    # ----- Tagging (single tolerance) -----
    if args.noninteractive:
        if args.tolerance is None:
            print("\n--noninteractive requires --tolerance.", file=sys.stderr)
            sys.exit(2)
        chosen_tol = float(args.tolerance)
    else:
        print("\nNext: tag headers with a grouping keyword for WBPP.")
        print("Enter a tolerance to tag with these centres, or 'q' to quit without tagging.")
        try:
            user_in = input("Tolerance (deg) or 'q': ").strip()
        except EOFError:
            user_in = "q"
        if user_in.lower() in {"q", "quit", "exit"}:
            print("Exiting after writing the report. No headers modified.")
            return
        chosen_tol = float(user_in)

    # Recompute centres jointly for the chosen tolerance (match report)
    centres_final, _, _ = optimise_joint(
        light_records=light_records,
        existing_global=existing_global,
        tol_deg=chosen_tol,
        mode=args.optimize,
        budget=args.budget
    )


    # Apply the same fuzzy merge policy for tagging
    existing_for_merge = sorted(set(existing_global)) or sorted(set(f for f in flats_all))
    merge_r = float(args.centre_merge_radius or 0.0)
    allow_loss = int(args.merge_allow_loss or 0)
    if existing_for_merge and merge_r > 0.0:
        centres_final, _ = try_merge_centres_to_existing(
            centres=centres_final,
            existing=existing_for_merge,
            tol=chosen_tol,
            light_angles=lights,
            merge_radius=merge_r,
            allow_loss=allow_loss
        )

    # Optional: snap to any real/provided flats if within tol (kept as before)
    if flats_all:
        centres_final = snap_centres_to_existing_flats(centres_final, flats_all, chosen_tol)


    print(f"\nTagging with tolerance ±{chosen_tol:.1f}°")
    print("Final centres used for tagging:")
    print("  " + (", ".join(f"{q(c):.2f}" for c in centres_final) if centres_final else "(none)"))
    print("Coverage windows:")
    for c in centres_final:
        print(f"  {fmt_angle(c)} → {fmt_window(c, chosen_tol)}")

    # Coverage per centre (quick)
    for c in centres_final:
        n_l = sum(1 for _, a, _ in light_records if circ_dist(a, c) <= chosen_tol)
        n_f = sum(1 for _, a, _ in flat_records  if circ_dist(a, c) <= chosen_tol)
        print(f"  {fmt_angle(c)}  ->  lights: {n_l:4d} | flats: {n_f:4d}")

    # Write headers
    wrote, assigned_total, skipped = 0, 0, 0
    for fp, ang, ftype, filt in file_records:
        # centre = assign_to_nearest(ang, centres_final, chosen_tol)
        centre = assign_to_nearest(ang, centres_final,
                           chosen_tol,
                           merge_r=float(args.centre_merge_radius or 0.0))

        if centre is None:
            skipped += 1
            continue
        assigned_total += 1
        label = f"ROT_{q(centre):.1f}"
        comment = f"Rotator group {q(centre):.1f} deg (+/-{chosen_tol:.1f} deg)"
        if args.dry_run:
            continue
        if write_rotgrp_keyword(fp, label, comment, args.keyword, args.overwrite):
            wrote += 1

    print("\n=== Tagging Summary ===")
    print(f"Frames with angle:          {len(file_records)}")
    print(f"Assigned to centres:        {assigned_total}")
    print(f"Skipped (no centre match):  {skipped}")
    if args.dry_run:
        print("Dry run: no headers modified.")
    else:
        print(f"Wrote {args.keyword}:             {wrote}")
    print("\nWBPP → Grouping → Keyword → ROTGRP")

if __name__ == "__main__":
    main()
