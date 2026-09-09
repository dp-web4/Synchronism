#!/usr/bin/env python3
"""SA-2 rung 4 — the within-galaxy remainder's structure (kimi, 2026-09-08).

Rung 3 left: within-galaxy residual Var 0.00538 dex^2 (per-point deviations
from per-galaxy means, after the 1-param function + per-galaxy Upsilon).
Rung 4 asks whether that remainder is STRUCTURE (radial trends, SB trends,
point-to-point coherence — unmodeled included-set physics) or NOISE (point
scatter at the instrument floor). Same cloud, cuts, hash split, and train-only
a0 as rung 3; Upsilon per galaxy from its own data (this rung analyzes
residual SHAPE; no predictive claim is made, so own-data Upsilon is honest
here — the strict-prediction spectrum lives in rung 2's V4 family).

Statistics, all on the 37 test galaxies:
  1. lag-1 autocorrelation of residuals along each curve (pooled): smooth
     structure correlates point-to-point; instrument noise does not.
  2. per-galaxy radial slope (residual vs R/Rdisk... MRT Rdisk unavailable in
     this lane's parser; use R/R_max of the galaxy — stated convention):
     distribution of slopes across galaxies + pooled mean residual by
     radius-decile bins.
  3. residual vs local SBdisk (log): pooled binned means.
  4. instrument comparison: per-galaxy residual Var vs per-galaxy instrument
     Var (errV-derived) — the above-floor factor.
"""
import glob
import hashlib
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KMS_KPC_TO_MS2 = (1e3) ** 2 / 3.086e19

def parse_mrt(path):
    gal = {}
    for line in open(path):
        f = line.split()
        if len(f) < 18:
            continue
        try:
            gal[f[0]] = {"T": int(f[1]), "Inc": float(f[5]), "L36": float(f[7]),
                         "MHI": float(f[13]), "Q": int(f[17])}
        except (ValueError, IndexError):
            continue
    return gal

def nu(x):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.clip(x, 1e-12, None))))

meta = parse_mrt(os.path.join(HERE, "SPARC_Lelli2016c.mrt"))
raw = {}
for path in sorted(glob.glob(os.path.join(HERE, "galaxies", "*.dat"))):
    name = os.path.basename(path)[:-4]
    m = meta.get(name)
    if m is None or m["Q"] > 2 or m["Inc"] < 30:
        continue
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    d = np.array(rows, dtype=float)
    R, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk = d.T[:7]
    ok = (np.abs(Vobs) > 3 * errV) & (R > 0)
    if ok.sum() >= 5:
        raw[name] = dict(R=R[ok], Vobs=Vobs[ok], errV=errV[ok], Vgas=Vgas[ok],
                         Vdisk=Vdisk[ok], Vbul=Vbul[ok], SBdisk=SBdisk[ok])

def gbar(n, yd, yb=0.7):
    r = raw[n]
    return (np.abs(r["Vgas"]) * r["Vgas"] + yd * np.abs(r["Vdisk"]) * r["Vdisk"]
            + yb * np.abs(r["Vbul"]) * r["Vbul"]) / r["R"] * KMS_KPC_TO_MS2

def gobs(n):
    return raw[n]["Vobs"] ** 2 / raw[n]["R"] * KMS_KPC_TO_MS2

def resid(go, gb, a0):
    return np.log10(go) - np.log10(gb * nu(gb / a0))

YGRID = np.linspace(0.05, 1.5, 59)

def fit_upsilon(n, a0):
    r = raw[n]; go = gobs(n)
    base = (np.abs(r["Vgas"]) * r["Vgas"] + 0.7 * np.abs(r["Vbul"]) * r["Vbul"]) / r["R"] * KMS_KPC_TO_MS2
    disk = np.abs(r["Vdisk"]) * r["Vdisk"] / r["R"] * KMS_KPC_TO_MS2
    best, best_y = np.inf, 0.5
    for y in YGRID:
        gb = base + y * disk
        ok = gb > 0
        if ok.sum() < 3:
            continue
        v = float(np.mean(resid(go[ok], gb[ok], a0) ** 2))
        if v < best:
            best, best_y = v, y
    return best_y

names = sorted(raw)
train = sorted(n for n in names if int(hashlib.md5(n.encode()).hexdigest(), 16) % 5 != 0)
test = [n for n in names if n not in train]
a0 = 1.2e-10
for _ in range(2):
    yd_tr = {n: fit_upsilon(n, a0) for n in train}
    best = np.inf
    for cand in np.linspace(3e-11, 4e-10, 60):
        ss, cnt = 0.0, 0
        for n in train:
            gb = gbar(n, yd_tr[n]); go = gobs(n); ok = gb > 0
            ss += np.sum(resid(go[ok], gb[ok], cand) ** 2); cnt += ok.sum()
        if ss / max(cnt, 1) < best:
            best, a0 = ss / cnt, cand

print(f"rung 4: {len(train)} train / {len(test)} test, a0_train = {a0:.3e}")

lag1_num, lag1_den = 0.0, 0.0
slopes, rad_bins, sb_bins = [], {}, {}
var_res_tot, var_inst_tot, npts = 0.0, 0.0, 0
for n in test:
    r = raw[n]
    y_te = fit_upsilon(n, a0)
    gb = gbar(n, y_te); go = gobs(n); ok = gb > 0
    res = resid(go[ok], gb[ok], a0)
    R = r["R"][ok]; SB = r["SBdisk"][ok]
    res = res - np.mean(res)                      # within-galaxy deviations
    inst = (2.0 * r["errV"][ok] / np.abs(r["Vobs"][ok]) / np.log(10.0)) ** 2
    var_res_tot += float(np.sum(res ** 2)); var_inst_tot += float(np.sum(inst)); npts += len(res)
    if len(res) >= 3:
        lag1_num += float(np.sum(res[:-1] * res[1:])); lag1_den += float(np.sum(res ** 2))
    x = R / R.max()
    if np.ptp(x) > 1e-9 and len(res) >= 5:
        slopes.append(np.polyfit(x, res, 1)[0])
    for xi, ri in zip(x, res):
        rad_bins.setdefault(min(int(xi * 10), 9), []).append(ri)
    for si, ri in zip(SB, res):
        if si > 0:
            sb_bins.setdefault(int(np.clip(np.log10(si) * 2 + 6, 0, 9)), []).append(ri)

print(f"\nwithin-galaxy remainder: {npts} points, Var {var_res_tot / npts:.5f} dex^2; "
      f"instrument Var {var_inst_tot / npts:.5f}; above-floor factor "
      f"{var_res_tot / max(var_inst_tot, 1e-12):.2f}x")
print(f"\n1. lag-1 autocorrelation (pooled): {lag1_num / max(lag1_den, 1e-12):+.4f}")
print(f"   (0 = point noise; >0 = smooth structure the model missed)")
sl = np.array(slopes)
print(f"\n2. radial slope vs R/Rmax: {len(sl)} galaxies; median {np.median(sl):+.4f} dex; "
      f"16-84% [{np.percentile(sl, 16):+.4f}, {np.percentile(sl, 84):+.4f}]; "
      f"share positive {(sl > 0).mean():.2f}")
print(f"   pooled mean residual by radius decile:")
for b in sorted(rad_bins):
    v = rad_bins[b]
    print(f"     R/Rmax [{b / 10:.1f},{(b + 1) / 10:.1f}): mean {np.mean(v):+.4f} dex (N={len(v)})")
print(f"\n3. pooled mean residual by SBdisk bin (log10 SB):")
for b in sorted(sb_bins):
    v = sb_bins[b]
    print(f"     bin {b}: mean {np.mean(v):+.4f} dex (N={len(v)})")
