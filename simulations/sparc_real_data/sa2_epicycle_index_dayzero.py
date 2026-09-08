#!/usr/bin/env python3
"""SA-2 day-zero baseline — the epicycle index of the galactic horizon (kimi, 2026-09-08).

The MRH-validity program's closure test, minimal form. The "anomaly" at the
galactic scale is A = g_obs - g_bar (Newtonian-baryons residual — the signal the
excluded-variable posit, dark matter, claims). The epicycle question: how much
of A is a FUNCTION OF THE INCLUDED SET alone? The RAR says: almost all of it —
g_obs is a tight one-parameter function of g_bar (an included variable). What
remains after that function is what either included-set structure (the 6-var
correction series) or an excluded sector (a halo) must price.

Index (registered in the SA-2 map doc):
  E = 1 - Var(log g_obs - log[g_bar * nu(g_bar/a0)]) / Var(log g_obs - log g_bar)
E -> 1: the anomaly is epicycle-class — fully a function of the included set,
         the excluded sector is invisible to the instrument beyond its
         correlation with it. E < 1: an unexplained remainder exists.

Data: simulations/sparc_real_data/galaxies/*.dat (175 SPARC galaxies,
Lelli-McGaugh-Schombert 2016 format). Cuts (stated, not tuned): |Vobs| > 3*errV,
R > 0, g_bar > 0, fiducial M/L = 1 for disk and bulge (the files' Vdisk/Vbul are
already at M/L=1, 3.6um; gas x1.33 for He is NOT applied — recorded as a known
simplification affecting absolute a0 fits, not the index shape).
"""
import glob
import os
import numpy as np

A0 = 1.2e-10          # m/s^2, McGaugh fiducial
KMS_KPC_TO_MS2 = (1e3) ** 2 / 3.086e19   # (km/s)^2 / kpc -> m/s^2

def nu(x):            # McGaugh interpolating function, x = g_bar/a0
    return 1.0 / (1.0 - np.exp(-np.sqrt(np.clip(x, 1e-12, None))))

g_obs_all, g_bar_all, err_all, vobs_all = [], [], [], []
per_galaxy = {}
for path in sorted(glob.glob(os.path.join(os.path.dirname(__file__),
                                          "galaxies", "*.dat"))):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    d = np.array(rows, dtype=float)
    R, Vobs, errV, Vgas, Vdisk, Vbul = d.T[:6]
    ok = (np.abs(Vobs) > 3 * errV) & (R > 0)
    R, Vobs, errV, Vgas, Vdisk, Vbul = (a[ok] for a in (R, Vobs, errV, Vgas, Vdisk, Vbul))
    if len(R) == 0:
        continue
    vbar2 = np.abs(Vgas) * Vgas + np.abs(Vdisk) * Vdisk + np.abs(Vbul) * Vbul
    g_obs = Vobs ** 2 / R * KMS_KPC_TO_MS2
    g_bar = vbar2 / R * KMS_KPC_TO_MS2
    ok2 = g_bar > 0
    g_obs_all.extend(g_obs[ok2]); g_bar_all.extend(g_bar[ok2]); err_all.extend(errV[ok2])
    vobs_all.extend(Vobs[ok2])
    per_galaxy[os.path.basename(path)[:-4]] = len(g_obs[ok2])

g_obs = np.array(g_obs_all); g_bar = np.array(g_bar_all)
x = g_bar / A0
anomaly = np.log10(g_obs) - np.log10(g_bar)            # what the halo posit claims
residual = np.log10(g_obs) - np.log10(g_bar * nu(x))   # what remains after the RAR

E = 1.0 - np.var(residual) / np.var(anomaly)
# instrument floor: log10 g_obs carries 2*errV/Vobs/ln(10) of measurement scatter;
# E raw understates horizon-pricing by that variance. E_corr removes the floor.
err_frac = 2.0 * (np.array(err_all) / np.abs(vobs_all)) / np.log(10.0)
var_inst = float(np.mean(err_frac ** 2))
E_corr = 1.0 - (np.var(residual) - var_inst) / np.var(anomaly)
print(f"galaxies used: {len(per_galaxy)}   points: {len(g_obs)}")
print(f"anomaly  log10(g_obs/g_bar):        std {np.std(anomaly):.4f} dex  "
      f"mean {np.mean(anomaly):+.4f} dex")
print(f"residual log10(g_obs/(g_bar*nu)):   std {np.std(residual):.4f} dex  "
      f"mean {np.mean(residual):+.4f} dex")
print(f"instrument floor (from errV):       std {np.sqrt(var_inst):.4f} dex")
print(f"EPICYCLE INDEX  E_raw  = 1 - Var(residual)/Var(anomaly)            = {E:.4f}")
print(f"EPICYCLE INDEX  E_corr = 1 - (Var(res)-Var(inst))/Var(anomaly)     = {E_corr:.4f}")
print()
# variance budget by acceleration decade (where does the remainder live?)
print(f"{'g_bar/a0 range':>18} {'N':>6} {'anom std':>9} {'resid std':>9} {'E_local':>8}")
for lo, hi in ((0.01, 0.1), (0.1, 1), (1, 10), (10, 100)):
    m = (x >= lo) & (x < hi)
    if m.sum() > 10:
        e_loc = 1 - np.var(residual[m]) / np.var(anomaly[m])
        print(f"{f'[{lo}, {hi})':>18} {m.sum():>6} {np.std(anomaly[m]):>9.4f} "
              f"{np.std(residual[m]):>9.4f} {e_loc:>8.4f}")
