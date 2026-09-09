#!/usr/bin/env python3
"""SA-2 rung 3 — the chained variance budget with the 6-var correction series,
on ONE point cloud with ONE cut set, with the identifiability column
(kimi, 2026-09-08).

Registered design (sub-arc map rung 3 + codex-review adoptions):
  - Same cloud/cuts as rung 2 corrected: Q<=2, Inc>=30, |Vobs|>3 sigma, SPARC
    Vgas as tabulated (He already included), 153 galaxies.
  - SAME galaxy split as V4 (md5(name)%5): 116 train / 37 test.
  - a0 fitted on TRAIN only; Upsilon_d per galaxy from that galaxy's own data.
  - STAGE 1 (the 1-param function): residual r = log g_obs - log[g_bar*nu].
  - STAGE 2 (the correction series): the Session-484 feature set —
    logV=log10(Vflat), logL=log10(L36/1e9), f_gas=1.33*MHI/(1.33*MHI+Ud*L36),
    c_V=<V_inner>/<V_outer> (halves by R), plus logV*c_V and logL*f_gas —
    OLS on TRAIN per-galaxy mean residuals, evaluated on TEST galaxies.
    NOTE: Upsilon is fitted BEFORE the series here (rung-2 order), so the
    series prices only what Upsilon left behind — sequential, never
    double-counted. The original 0.938 was measured with fixed Upsilon; the
    chained number WILL be smaller, and that is the point.
  - BUDGET (test points): anomaly Var = instrument + 1-param priced +
    between-galaxy priced by series + within-galaxy + REMAINDER.
  - IDENTIFIABILITY column (codex T5): the remainder is a COMPATIBLE
    INTERVAL [0, remainder] for excluded-sector signal — this channel cannot
    distinguish unmodeled included-set structure from a sector. Printed as a
    column, not a footnote.
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
            gal[f[0]] = {"T": int(f[1]), "D": float(f[2]), "Inc": float(f[5]),
                         "L36": float(f[7]), "MHI": float(f[13]),
                         "Vflat": float(f[15]), "Q": int(f[17])}
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
    R, Vobs, errV, Vgas, Vdisk, Vbul = d.T[:6]
    ok = (np.abs(Vobs) > 3 * errV) & (R > 0)
    R, Vobs, errV, Vgas, Vdisk, Vbul = (a[ok] for a in (R, Vobs, errV, Vgas, Vdisk, Vbul))
    if len(R) >= 5:
        raw[name] = dict(R=R, Vobs=Vobs, errV=errV, Vgas=Vgas, Vdisk=Vdisk, Vbul=Vbul)

def gbar_parts(n, yd, yb=0.7):
    r = raw[n]
    return (np.abs(r["Vgas"]) * r["Vgas"] / r["R"] * KMS_KPC_TO_MS2
            + yd * np.abs(r["Vdisk"]) * r["Vdisk"] / r["R"] * KMS_KPC_TO_MS2
            + yb * np.abs(r["Vbul"]) * r["Vbul"] / r["R"] * KMS_KPC_TO_MS2)

def gobs(n):
    r = raw[n]
    return r["Vobs"] ** 2 / r["R"] * KMS_KPC_TO_MS2

def resid(go, gb, a0):
    return np.log10(go) - np.log10(gb * nu(gb / a0))

YGRID = np.linspace(0.05, 1.5, 59)

def fit_upsilon(go, n, a0, yb=0.7):
    r = raw[n]
    disk = np.abs(r["Vdisk"]) * r["Vdisk"] / r["R"] * KMS_KPC_TO_MS2
    base = (np.abs(r["Vgas"]) * r["Vgas"] / r["R"] * KMS_KPC_TO_MS2
            + yb * np.abs(r["Vbul"]) * r["Vbul"] / r["R"] * KMS_KPC_TO_MS2)
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
print(f"galaxies: {len(names)} total ({len(raw)} with >=5 pts), "
      f"{len(train)} train / {len(test)} test")

# a0 on train (with per-galaxy Upsilon), two alternations
a0 = 1.2e-10
yd = {}
for _ in range(2):
    yd = {n: fit_upsilon(gobs(n), n, a0) for n in train}
    best, best_a = np.inf, a0
    for cand in np.linspace(3e-11, 4e-10, 60):
        ss, cnt = 0.0, 0
        for n in train:
            gb = gbar_parts(n, yd[n]); go = gobs(n); ok = gb > 0
            ss += np.sum(resid(go[ok], gb[ok], cand) ** 2); cnt += ok.sum()
        if ss / max(cnt, 1) < best:
            best, best_a = ss / cnt, cand
    a0 = best_a
for n in test:  # test galaxies: Upsilon from their own data at train a0
    yd[n] = fit_upsilon(gobs(n), n, a0)
print(f"a0_train = {a0:.3e} m/s^2; Upsilon_d median {np.median(list(yd.values())):.3f}")

# per-galaxy residuals, features
def features(n):
    r = raw[n]; m = meta[n]
    mid = len(r["R"]) // 2
    cv = np.mean(np.abs(r["Vobs"][:mid])) / max(np.mean(np.abs(r["Vobs"][mid:])), 1e-9)
    mgas = 1.33 * m["MHI"]
    fg = mgas / (mgas + yd[n] * m["L36"])
    lv, ll = np.log10(max(m["Vflat"], 1e-3)), np.log10(max(m["L36"], 1e-6))
    return np.array([1.0, lv, ll, cv, fg, lv * cv, ll * fg])

per_gal = {}
for n in names:
    gb = gbar_parts(n, yd[n]); go = gobs(n); ok = gb > 0
    res = resid(go[ok], gb[ok], a0)
    per_gal[n] = dict(mean=float(np.mean(res)), res=res,
                      errV=raw[n]["errV"][ok], vobs=raw[n]["Vobs"][ok],
                      X=features(n))

Xtr = np.array([per_gal[n]["X"] for n in train])
ytr = np.array([per_gal[n]["mean"] for n in train])
beta, *_ = np.linalg.lstsq(Xtr, ytr, rcond=None)
# train LOO for the series itself
loo_pred = np.empty_like(ytr)
for i in range(len(train)):
    m_ = np.ones(len(train), bool); m_[i] = False
    b, *_ = np.linalg.lstsq(Xtr[m_], ytr[m_], rcond=None)
    loo_pred[i] = Xtr[i] @ b
loo_r2 = 1 - np.var(ytr - loo_pred) / np.var(ytr)
print(f"series on train per-galaxy mean residuals: R2 {1 - np.var(ytr - Xtr @ beta) / np.var(ytr):.4f}, "
      f"LOO R2 {loo_r2:.4f} (session-484 measured 0.945/0.938 with FIXED Upsilon;")

# ---- the chained budget on TEST points
anom, res_all, inst = [], [], []
m_te, mhat_te = [], []
for n in test:
    g = per_gal[n]
    gb = gbar_parts(n, yd[n]); go = gobs(n); ok = gb > 0
    anom.extend(np.log10(go[ok]) - np.log10(gb[ok]))
    res_all.extend(g["res"])
    inst.extend((2.0 * g["errV"] / np.abs(g["vobs"]) / np.log(10.0)) ** 2)
    m_te.append(g["mean"]); mhat_te.append(float(g["X"] @ beta))
anom, res_all, inst = map(np.array, (anom, res_all, inst))
m_te, mhat_te = np.array(m_te), np.array(mhat_te)

var_anom = np.var(anom)
var_res = np.var(res_all)
var_inst = float(np.mean(inst))
between = np.var(m_te)
between_after = np.var(m_te - mhat_te)
within = var_res - between          # law of total variance (approx: means unweighted)
series_priced = between - between_after
remainder = var_res - var_inst - series_priced
print(f"\n=== CHAINED BUDGET (37 test galaxies, {len(anom)} points) ===")
print(f"anomaly Var                       {var_anom:.5f} dex^2")
print(f"1-param function+Ud priced        {var_anom - var_res:.5f} ({(var_anom - var_res) / var_anom * 100:.1f}%)")
print(f"instrument floor (errV)           {var_inst:.5f} ({var_inst / var_anom * 100:.1f}%)")
print(f"between-galaxy Var (before series){between:.5f}; series prices {series_priced:.5f} "
      f"-> after {between_after:.5f}")
print(f"within-galaxy Var                 {within:.5f}")
print(f"REMAINDER (res - inst - series)   {remainder:.5f} ({remainder / var_anom * 100:.1f}%)")
print(f"E_corr after series = 1-(res-inst-series)/anom = {1 - remainder / var_anom:.4f}")
print(f"series out-of-sample: Var(test means) {between:.5f} -> {between_after:.5f} "
      f"({'prices' if series_priced > 0 else 'DEGRADES'} "
      f"{abs(series_priced) / max(between, 1e-12) * 100:.1f}% of between-galaxy Var)")
print(f"\nIDENTIFIABILITY (codex T5): the remainder {remainder:.5f} is a COMPATIBLE")
print(f"INTERVAL [0, {remainder:.5f}] for excluded-sector signal. This channel")
print(f"(rotation curves + 3.6um photometry) cannot distinguish unmodeled")
print(f"included-set structure from a sector. The interval is the column; the")
print(f"answer requires a new channel or a restricted model class.")
print(f"\ntest residual mean {np.mean(res_all):+.4f} dex (bias column, codex T6)")
