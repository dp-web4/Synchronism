#!/usr/bin/env python3
"""SA-2 rung 2 + SA-1b first cut — the epicycle index under honest cuts,
per-galaxy Upsilon with LOO and GALAXY holdout, and the residual-scatter
scaling (kimi, 2026-09-08; errata from codex review same day incorporated).

RUNG 2 (honest E): quality cuts Q<=2 and Inc>=30 deg from SPARC_Lelli2016c.mrt,
point cut |Vobs|>3 sigma. SPARC's tabulated Vgas ALREADY includes the x1.33 He
factor (Lelli+2016 section 3.3) — no gas correction is applied anywhere in this
file (the first commit multiplied by 1.33 again: a double-count, caught in
codex's review 2026-09-08). Variants:
  V1 fiducial:   Upsilon_d=0.5, Upsilon_b=0.7, a0=1.2e-10 (literature)
  V2 fitted:     Upsilon_d per galaxy (grid, 0.05..1.5), Upsilon_b=0.7,
                 a0 global fit, two alternations — the standard practice
  V3 point-LOO:  V2 but each point's residual uses Upsilon_d refit WITHOUT that
                 point. Measures within-galaxy interpolation stability with a
                 shared global a0 — NOT new-galaxy prediction (codex review).
  V4 GALAXY holdout: the named task — predict a NEW galaxy's rotation curve
                 from its baryonic data + the global law. a0 fitted on 116
                 train galaxies only; Upsilon_d of each of the 37 held-out
                 galaxies fitted from that galaxy's own points (a legitimate
                 per-galaxy baryonic property); E over held-out points.
E_raw / E_corr as in day zero (instrument floor from errV), with residual
means reported alongside (E is variance-only and bias-blind — codex T6).

SA-1b FIRST CUT (the exponent question, galactic): after V2, per-galaxy
residual std sigma_g vs the included-DOF proxy M_star = Upsilon_d * L[3.6].
A -1/2 log-log slope = Gaussian counting of independent units; 0 = NOT counting
of independent units. Codex's narrowing applies: slope ~0 does not by itself
establish non-Gaussianity, causal structure, or escape from K2 — dependence,
distribution shape, and physical novelty are distinct questions; the slope
says only that the scatter does not behave like averaging of independent
identical units of this proxy.
"""
import glob
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
            gal[f[0]] = {
                "T": int(f[1]), "D": float(f[2]), "Inc": float(f[5]),
                "L36": float(f[7]), "SBeff": float(f[10]), "MHI": float(f[13]),
                "Q": int(f[17]),
            }
        except (ValueError, IndexError):
            continue
    return gal

def nu(x):
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.clip(x, 1e-12, None))))

def load_points(meta, yd_map=None, yb=0.7, he=1.0):   # he: SPARC Vgas already includes the 1.33 He factor (Lelli+2016 §3.3 — erratum 2026-09-08, codex review; the x1.33 in the first commit of this script was a double-count)
    """yd_map: galaxy -> Upsilon_d; default 0.5. Returns per-point dict of arrays."""
    out = {k: [] for k in ("g_obs", "g_bar", "errV", "vobs", "gal")}
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
        if not len(R):
            continue
        yd = (yd_map or {}).get(name, 0.5)
        vbar2 = he * np.abs(Vgas) * Vgas + yd * np.abs(Vdisk) * Vdisk + yb * np.abs(Vbul) * Vbul
        g_obs = Vobs ** 2 / R * KMS_KPC_TO_MS2
        g_bar = vbar2 / R * KMS_KPC_TO_MS2
        ok2 = g_bar > 0
        for k, arr in (("g_obs", g_obs[ok2]), ("g_bar", g_bar[ok2]),
                       ("errV", errV[ok2]), ("vobs", Vobs[ok2])):
            out[k].extend(arr)
        out["gal"].extend([name] * ok2.sum())
    return {k: np.array(v) for k, v in out.items()}

def resid(g_obs, g_bar, a0):
    return np.log10(g_obs) - np.log10(g_bar * nu(g_bar / a0))

def index_E(anomaly, residual, errV, vobs):
    var_inst = float(np.mean((2.0 * (errV / np.abs(vobs)) / np.log(10.0)) ** 2))
    e_raw = 1 - np.var(residual) / np.var(anomaly)
    e_corr = 1 - (np.var(residual) - var_inst) / np.var(anomaly)
    return e_raw, e_corr, np.sqrt(var_inst)

YGRID = np.linspace(0.05, 1.5, 59)

def fit_upsilon(g_obs, g_bar_parts, a0, exclude=None):
    """1-D grid search for Upsilon_d minimizing residual variance (per galaxy)."""
    he_part, disk_part, bulge_part = g_bar_parts
    best, best_y = np.inf, 0.5
    for y in YGRID:
        gb = he_part + y * disk_part + bulge_part
        if exclude is not None:
            gb, go = gb[~exclude], g_obs[~exclude]
        else:
            go = g_obs
        ok = gb > 0
        if ok.sum() < 3:
            continue
        r = resid(go[ok], gb[ok], a0)
        v = float(np.mean(r ** 2))   # RMS, not variance: the per-galaxy mean counts
        if v < best:
            best, best_y = v, y
    return best_y

meta = parse_mrt(os.path.join(HERE, "SPARC_Lelli2016c.mrt"))
print(f"galaxies in MRT: {len(meta)}; after Q<=2 & Inc>=30: "
      f"{sum(1 for m in meta.values() if m['Q'] <= 2 and m['Inc'] >= 30)}")

# ---------------- V1: fiducial
p = load_points(meta)
a0_fid = 1.2e-10
anom = np.log10(p["g_obs"]) - np.log10(p["g_bar"])
res = resid(p["g_obs"], p["g_bar"], a0_fid)
e1 = index_E(anom, res, p["errV"], p["vobs"])
print(f"\nV1 fiducial (Ud=0.5, Ub=0.7, Vgas as tabulated incl. He, a0=1.2e-10): "
      f"N={len(anom)} pts, {len(set(p['gal']))} galaxies")
print(f"  anomaly std {np.std(anom):.4f} | residual std {np.std(res):.4f} "
      f"| instrument {e1[2]:.4f} | E_raw {e1[0]:.4f} | E_corr {e1[1]:.4f}")

# ---------------- V2: per-galaxy Upsilon_d + global a0
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
    if len(R) >= 3:
        raw[name] = dict(R=R, Vobs=Vobs, errV=errV, Vgas=Vgas, Vdisk=Vdisk, Vbul=Vbul)

def parts(name, yb=0.7, he=1.0):
    r = raw[name]
    return (he * np.abs(r["Vgas"]) * r["Vgas"] / r["R"] * KMS_KPC_TO_MS2,
            np.abs(r["Vdisk"]) * r["Vdisk"] / r["R"] * KMS_KPC_TO_MS2,
            yb * np.abs(r["Vbul"]) * r["Vbul"] / r["R"] * KMS_KPC_TO_MS2)

# ---------------- V2 + V3 in one pass over raw (consistent point set)
a0 = a0_fid
yd_map = {}
for _ in range(2):
    yd_map = {n: fit_upsilon(raw[n]["Vobs"] ** 2 / raw[n]["R"] * KMS_KPC_TO_MS2,
                             parts(n), a0) for n in raw}
    best, best_a = np.inf, a0
    for cand in np.linspace(3e-11, 4e-10, 60):
        tot = 0.0; cnt = 0
        for n in raw:
            hp, dp, bp = parts(n)
            gb = hp + yd_map[n] * dp + bp
            go = raw[n]["Vobs"] ** 2 / raw[n]["R"] * KMS_KPC_TO_MS2
            ok = gb > 0
            tot += np.sum(resid(go[ok], gb[ok], cand) ** 2); cnt += ok.sum()
        if tot / max(cnt, 1) < best:
            best, best_a = tot / cnt, cand
    a0 = best_a

anom2, res2, res3, err2, vob2 = [], [], [], [], []
for n in raw:
    r = raw[n]
    go_all = r["Vobs"] ** 2 / r["R"] * KMS_KPC_TO_MS2
    hp, dp, bp = parts(n)
    gb_all = hp + yd_map[n] * dp + bp
    ok = gb_all > 0
    for j in np.where(ok)[0]:
        anom2.append(np.log10(go_all[j]) - np.log10(gb_all[j]))
        res2.append(float(resid(np.array([go_all[j]]), np.array([gb_all[j]]), a0)[0]))
        exc = np.zeros(len(go_all), bool); exc[j] = True
        y_loo = fit_upsilon(go_all, (hp, dp, bp), a0, exclude=exc)
        gb_loo = hp[j] + y_loo * dp[j] + bp[j]
        res3.append(np.log10(go_all[j]) - np.log10(gb_loo * nu(gb_loo / a0)))
        err2.append(r["errV"][j]); vob2.append(r["Vobs"][j])
anom2, res2, res3, err2, vob2 = map(np.array, (anom2, res2, res3, err2, vob2))
e2 = index_E(anom2, res2, err2, vob2)
yvals = np.array(list(yd_map.values()))
print(f"\nV2 fitted (Ud per galaxy, a0 global): a0 = {a0:.3e} m/s^2; "
      f"Ud median {np.median(yvals):.3f}, 16-84% [{np.percentile(yvals,16):.3f}, "
      f"{np.percentile(yvals,84):.3f}]")
print(f"  anomaly std {np.std(anom2):.4f} | residual std {np.std(res2):.4f} "
      f"mean {np.mean(res2):+.4f} | instrument {e2[2]:.4f} "
      f"| E_raw {e2[0]:.4f} | E_corr {e2[1]:.4f}")

# ---------------- V3: point-LOO Upsilon (residuals built in the same pass above)
e3 = index_E(anom2, res3, err2, vob2)
print(f"\nV3 point-LOO (Ud refit per held-out point; measures within-galaxy")
print(f"  interpolation stability with shared global a0 — NOT new-galaxy prediction):")
print(f"  residual std {np.std(res3):.4f} mean {np.mean(res3):+.4f} "
      f"| E_raw {e3[0]:.4f} | E_corr {e3[1]:.4f}")

# ---------------- V4 family: GALAXY-level holdout, four increasingly strict
# variants (V4 relabeled + V4a/V4b/V4c added after codex's V4-leakage flag,
# 2026-09-08). The named strict task: predict a NEW galaxy's rotation curve
# using ONLY its baryonic/photometric data + the train-fitted law.
# V4  (relabeled): test-galaxy rotation-curve calibration — Upsilon_d fitted
#      from the held-out galaxy's OWN ROTATION CURVE (target-dependent; kept
#      for comparison, no longer called prediction).
# V4a (strict prior): Upsilon_d = 0.5 fixed for all test galaxies — no test
#      information of any kind in the evaluation.
# V4b (strict mapping): Upsilon_d = f(galaxy properties), f fitted on train
#      only (Session-484 feature set, ridge), applied to test photometry.
# V4c (few-shot, LABELED): first half of each test galaxy's points (by R)
#      calibrates Upsilon_d, the OTHER half is evaluated — few-shot
#      within-galaxy prediction, never substituted for V4a/V4b.
# Upstream-dependence repair (codex point 2): the train-side Upsilon/a0
# alternation restarts from the fiducial a0 STRICTLY within train galaxies.
import hashlib
names = sorted(raw)
train = {n for n in names if int(hashlib.md5(n.encode()).hexdigest(), 16) % 5 != 0}
test = [n for n in names if n not in train]

a0_tr = 1.2e-10
for _ in range(2):
    yd_tr = {n: fit_upsilon(raw[n]["Vobs"] ** 2 / raw[n]["R"] * KMS_KPC_TO_MS2,
                            parts(n), a0_tr) for n in sorted(train)}
    best = np.inf
    for cand in np.linspace(3e-11, 4e-10, 60):
        tot = 0.0; cnt = 0
        for n in sorted(train):
            hp, dp, bp = parts(n)
            gb = hp + yd_tr[n] * dp + bp
            go = raw[n]["Vobs"] ** 2 / raw[n]["R"] * KMS_KPC_TO_MS2
            ok = gb > 0
            tot += np.sum(resid(go[ok], gb[ok], cand) ** 2); cnt += ok.sum()
        if tot / max(cnt, 1) < best:
            best, a0_tr = tot / cnt, cand

def galaxy_features(n, yd_val=None):
    """PHOTOMETRIC + HI features only — nothing rotation-curve-derived
    (Vflat, c_V are RC quantities and would leak the target)."""
    m = meta[n]
    mgas = 1.33 * m["MHI"]
    fg = mgas / (mgas + 0.5 * m["L36"])   # reference Upsilon=0.5: photometric
    ll = np.log10(max(m["L36"], 1e-6))
    lsb = np.log10(max(m["SBeff"], 1e-6))
    return np.array([1.0, ll, lsb, fg, ll * fg, m["T"] / 10.0])

def eval_test(yd_of, split_half=False):
    an, rs, er, vo = [], [], [], []
    for n in test:
        r = raw[n]
        go_all = r["Vobs"] ** 2 / r["R"] * KMS_KPC_TO_MS2
        hp, dp, bp = parts(n)
        if split_half:
            half = len(go_all) // 2
            y_te = fit_upsilon(go_all[:half], (hp[:half], dp[:half], bp[:half]), a0_tr)
            idxs = range(half, len(go_all))
        else:
            y_te = yd_of(n)
            idxs = range(len(go_all))
        gb_all = hp + y_te * dp + bp
        for j in idxs:
            if gb_all[j] <= 0:
                continue
            an.append(np.log10(go_all[j]) - np.log10(gb_all[j]))
            rs.append(np.log10(go_all[j]) - np.log10(gb_all[j] * nu(gb_all[j] / a0_tr)))
            er.append(r["errV"][j]); vo.append(r["Vobs"][j])
    an, rs, er, vo = map(np.array, (an, rs, er, vo))
    e = index_E(an, rs, er, vo)
    return e, np.std(rs), np.mean(rs), len(rs)

e4, s4, m4, n4 = eval_test(lambda n: fit_upsilon(
    raw[n]["Vobs"] ** 2 / raw[n]["R"] * KMS_KPC_TO_MS2, parts(n), a0_tr))
print(f"\nV4 RELABELED — test-galaxy ROTATION-CURVE CALIBRATION ({len(train)} train /"
      f" {len(test)} test; a0_train {a0_tr:.3e}; NOT the strict prediction task):")
print(f"  residual std {s4:.4f} mean {m4:+.4f} | E_raw {e4[0]:.4f} | E_corr {e4[1]:.4f}")

e4a, s4a, m4a, n4a = eval_test(lambda n: 0.5)
print(f"\nV4a STRICT PRIOR (Upsilon=0.5 fixed; zero test information):")
print(f"  residual std {s4a:.4f} mean {m4a:+.4f} | E_raw {e4a[0]:.4f} | E_corr {e4a[1]:.4f}")

Xtr = np.array([galaxy_features(n, yd_tr[n]) for n in sorted(train)])
ytr = np.array([yd_tr[n] for n in sorted(train)])
lam = 1e-3 * np.trace(Xtr.T @ Xtr) / Xtr.shape[1]
beta_y = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1]), Xtr.T @ ytr)
yd_pred = {n: float(np.clip(galaxy_features(n, 0.5) @ beta_y, 0.05, 1.5)) for n in test}
e4b, s4b, m4b, n4b = eval_test(lambda n: yd_pred[n])
mapped = np.array(list(yd_pred.values()))
print(f"\nV4b STRICT MAPPING (Upsilon = f(properties), f ridge on train only; "
      f"mapped Upsilon median {np.median(mapped):.3f}):")
print(f"  residual std {s4b:.4f} mean {m4b:+.4f} | E_raw {e4b[0]:.4f} | E_corr {e4b[1]:.4f}")
print(f"  (mapping quality: train R2 for Upsilon = "
      f"{1 - np.var(ytr - Xtr @ beta_y) / np.var(ytr):.3f})")

e4c, s4c, m4c, n4c = eval_test(None, split_half=True)
print(f"\nV4c FEW-SHOT (first half by R calibrates Upsilon, second half scored; "
      f"{n4c} scored points; LABELED within-galaxy few-shot):")
print(f"  residual std {s4c:.4f} mean {m4c:+.4f} | E_raw {e4c[0]:.4f} | E_corr {e4c[1]:.4f}")

# ---------------- SA-1b first cut: sigma_g vs M_star
print("\nSA-1b first cut: per-galaxy residual std vs M_star = Ud * L[3.6]")
ms, sg = [], []
i0 = 0
for n in raw:
    cnt = int((parts(n)[0] + yd_map[n] * parts(n)[1] + parts(n)[2] > 0).sum())
    seg = res2[i0:i0 + cnt]; i0 += cnt
    if cnt < 5:
        continue
    ms.append(yd_map[n] * meta[n]["L36"])     # 1e9 solMass
    sg.append(np.std(seg))
ms, sg = np.array(ms), np.array(sg)
okm = (ms > 0) & (sg > 0)
slope, intercept = np.polyfit(np.log10(ms[okm]), np.log10(sg[okm]), 1)
print(f"  galaxies: {okm.sum()}; log-log slope sigma_g vs M_star: {slope:+.3f}")
print(f"  (-0.5 = Gaussian counting of independent units; 0 = structure, not counting)")
# split by M_star terciles for a look at the shape
qs = np.quantile(np.log10(ms[okm]), [1/3, 2/3])
lm = np.log10(ms[okm])
for lab, m_ in (("low M*", lm < qs[0]), ("mid M*", (lm >= qs[0]) & (lm < qs[1])),
                ("high M*", lm >= qs[1])):
    print(f"  {lab}: median sigma_g {np.median(sg[okm][m_]):.4f} dex, N={m_.sum()}")
