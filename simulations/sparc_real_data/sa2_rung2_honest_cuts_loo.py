#!/usr/bin/env python3
"""SA-2 rung 2 + SA-1b first cut — the epicycle index under honest cuts,
per-galaxy Upsilon with LOO, and the residual-scatter scaling (kimi, 2026-09-08).

RUNG 2 (honest E): quality cuts Q<=2 and Inc>=30 deg from SPARC_Lelli2016c.mrt,
He x1.33 gas correction, point cut |Vobs|>3 sigma. Three variants:
  V1 fiducial:   Upsilon_d=0.5, Upsilon_b=0.7, a0=1.2e-10 (literature)
  V2 fitted:     Upsilon_d per galaxy (grid, 0.05..1.5), Upsilon_b=0.7,
                 a0 global fit, two alternations — the standard practice
  V3 LOO:        V2 but each point's residual uses Upsilon_d refit WITHOUT that
                 point (a0 kept global; one global param over ~3000 pts)
E_raw / E_corr as in day zero (instrument floor from errV).

SA-1b FIRST CUT (the exponent question, galactic): after V2, per-galaxy
residual std sigma_g vs the included-DOF proxy M_star = Upsilon_d * L[3.6].
A -1/2 log-log slope = Gaussian counting of independent units; 0 = the residual
is STRUCTURE (series terms or sector), not counting fluctuation. Registered
expectation per the program: nearer 0 than -1/2 — but measured, not assumed.
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

def load_points(meta, yd_map=None, yb=0.7, he=1.33):
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
print(f"\nV1 fiducial (Ud=0.5, Ub=0.7, He x1.33, a0=1.2e-10): "
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

def parts(name, yb=0.7, he=1.33):
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
      f"| instrument {e2[2]:.4f} | E_raw {e2[0]:.4f} | E_corr {e2[1]:.4f}")

# ---------------- V3: LOO Upsilon (residuals built in the same pass above)
e3 = index_E(anom2, res3, err2, vob2)
print(f"\nV3 LOO (Ud refit per held-out point):")
print(f"  residual std {np.std(res3):.4f} | E_raw {e3[0]:.4f} | E_corr {e3[1]:.4f}")

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
