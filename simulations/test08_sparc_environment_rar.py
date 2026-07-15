#!/usr/bin/env python3
"""
TEST-08 (registered run): SPARC Environment Catalog
====================================================
Date: 2026-07-14
Registered in: Research/EXPERIMENTAL_TEST_CATALOG.md (TEST-08, lines ~128-141)

Prediction:    RAR scatter correlates with galactic environment (cluster vs
               field vs void), not just internal structure.
Expected:      Environment explains >20% of RAR scatter.
Falsification: Environment correlation < 0.3 (r^2 < 0.09).
MOND (no EFE): zero environment dependence.

PRE-DECLARED before computing any correlation (see exploration note):
  PRIMARY statistic = Pearson r between per-galaxy RAR offset (mean dex
  residual from the pooled RAR fit) and distance-corrected log density
  from redshift-space cylinder counts in Cosmicflows-4
  (projected radius 2 Mpc, |dV| < 500 km/s, self-excluded).
  r^2 < 0.09  -> TEST-08 REFUTED as registered.
  r^2 >= 0.20 -> registered expectation met.
  0.09 <= r^2 < 0.20 -> alive but under-strength (per catalog wording).
  Secondaries (context, not verdict): 5 Mpc sphere count, 5th-NN 3D density,
  Virgo-centric distance, cluster/field/void tercile ANOVA, and partial
  correlation controlling for internal structure (log L36, log SBeff, log D).

Data:
  simulations/sparc_real_data/SPARC_Lelli2016c.mrt      (175 galaxies, table1)
  simulations/sparc_real_data/MassModels_Lelli2016c.mrt (3391 points, table2)
  simulations/sparc_real_data/sparc_ned_coordinates.json (NED, fetched 2026-07-14)
  data/cf4/table2.dat                                    (Cosmicflows-4, 55877 gal)

RAR construction follows McGaugh, Lelli & Schombert 2016 (PRL 117, 201101):
  M/L[3.6] = 0.5 (disk), 0.7 (bulge); quality cuts Q<=2, i>=30 deg,
  e_Vobs/Vobs < 0.1; g_bar > 0; g = V^2/R.
"""

import json
import os
import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
SPARC_T1 = os.path.join(HERE, 'sparc_real_data', 'SPARC_Lelli2016c.mrt')
SPARC_T2 = os.path.join(HERE, 'sparc_real_data', 'MassModels_Lelli2016c.mrt')
COORDS = os.path.join(HERE, 'sparc_real_data', 'sparc_ned_coordinates.json')
CF4 = os.path.join(HERE, '..', 'data', 'cf4', 'table2.dat')

KM2_S2_PER_KPC_TO_M_S2 = 1.0e6 / 3.0857e19   # (km/s)^2/kpc -> m/s^2
H0 = 74.6                                     # CF4 zero-point, km/s/Mpc
ML_DISK, ML_BUL = 0.5, 0.7

rng_note = []

# ---------------------------------------------------------------- SPARC table1
def parse_table1():
    gals = {}
    lines = open(SPARC_T1).readlines()
    last_dash = max(i for i, l in enumerate(lines) if l.startswith('-----'))
    for l in lines[last_dash + 1:]:
        if len(l) < 95:
            continue
        name = l[0:11].strip()
        f = l[11:].split()
        # T D e_D f_D Inc e_Inc L36 e_L36 Reff SBeff Rdisk SBdisk MHI RHI
        # Vflat e_Vflat Q Ref...
        gals[name] = dict(
            T=int(f[0]), D=float(f[1]),
            inc=float(f[4]), L36=float(f[6]),
            SBeff=float(f[9]), Vflat=float(f[14]),
            Q=int(f[16]),
        )
    return gals

# ------------------------------------------------------------- SPARC table2
def parse_table2():
    pts = {}
    lines = open(SPARC_T2).readlines()
    last_dash = max(i for i, l in enumerate(lines) if l.startswith('-----'))
    for l in lines[last_dash + 1:]:
        if len(l) < 60:
            continue
        name = l[0:11].strip()
        f = l[11:].split()
        if len(f) < 9:
            continue
        try:
            # D R Vobs e_Vobs Vgas Vdisk Vbul SBdisk SBbul
            row = dict(R=float(f[1]), Vobs=float(f[2]),
                       eV=float(f[3]), Vgas=float(f[4]),
                       Vdisk=float(f[5]), Vbul=float(f[6]))
        except ValueError:
            continue
        pts.setdefault(name, []).append(row)
    return pts

# ---------------------------------------------------------------------- CF4
def parse_cf4():
    pgc, vcmb, dm, ra, dec = [], [], [], [], []
    with open(CF4) as f:
        for l in f:
            toks = l.split()
            if len(toks) < 12:
                continue
            try:
                # head: PGC 1PGC T17 Vcmb DM e_DM ... ; tail (always present):
                # RAdeg DEdeg GLON GLAT SGL SGB
                head = l[:41].split()
                vcmb.append(float(head[3]))
                dm.append(float(head[4]))
                ra.append(float(toks[-6]))
                dec.append(float(toks[-5]))
                pgc.append(int(head[0]))
            except (ValueError, IndexError):
                continue
    vcmb, dm = np.array(vcmb), np.array(dm)
    ra, dec = np.radians(np.array(ra)), np.radians(np.array(dec))
    dist = 10 ** ((dm - 25.0) / 5.0)  # Mpc
    return dict(pgc=np.array(pgc), vcmb=vcmb, dist=dist, ra=ra, dec=dec)

def unit_vec(ra, dec):
    return np.stack([np.cos(dec) * np.cos(ra),
                     np.cos(dec) * np.sin(ra),
                     np.sin(dec)], axis=-1)

# ------------------------------------------------------------------ RAR side
t1 = parse_table1()
t2 = parse_table2()
coords = json.load(open(COORDS))

g_bar_all, g_obs_all, gal_of_pt = [], [], []
for name, rows in t2.items():
    meta = t1.get(name)
    if meta is None or meta['Q'] > 2 or meta['inc'] < 30:
        continue
    for r in rows:
        if r['R'] <= 0 or r['Vobs'] <= 0 or r['eV'] / r['Vobs'] >= 0.1:
            continue
        gb = (r['Vgas'] * abs(r['Vgas'])
              + ML_DISK * r['Vdisk'] * abs(r['Vdisk'])
              + ML_BUL * r['Vbul'] * abs(r['Vbul'])) / r['R']
        if gb <= 0:
            continue
        g_bar_all.append(gb * KM2_S2_PER_KPC_TO_M_S2)
        g_obs_all.append(r['Vobs'] ** 2 / r['R'] * KM2_S2_PER_KPC_TO_M_S2)
        gal_of_pt.append(name)

g_bar_all = np.array(g_bar_all)
g_obs_all = np.array(g_obs_all)
gal_of_pt = np.array(gal_of_pt)
print(f'RAR points after cuts: {len(g_bar_all)} '
      f'({len(set(gal_of_pt))} galaxies)')

def rar_fn(gbar, gdag):
    return gbar / (1.0 - np.exp(-np.sqrt(gbar / gdag)))

# fit g-dagger on pooled points (log-space least squares)
def loss(log_gdag):
    pred = rar_fn(g_bar_all, 10 ** log_gdag)
    return np.sum((np.log10(g_obs_all) - np.log10(pred)) ** 2)

res = optimize.minimize_scalar(loss, bounds=(-11, -9), method='bounded')
gdag = 10 ** res.x
print(f'fitted g-dagger = {gdag:.3e} m/s^2 (McGaugh+16: 1.20e-10)')

resid = np.log10(g_obs_all) - np.log10(rar_fn(g_bar_all, gdag))
print(f'pooled point scatter = {np.std(resid):.3f} dex '
      f'(McGaugh+16: ~0.11 dex)')

# per-galaxy offsets (>=3 points)
gal_names, gal_offset, gal_npts = [], [], []
for name in sorted(set(gal_of_pt)):
    m = gal_of_pt == name
    if m.sum() < 3:
        continue
    gal_names.append(name)
    gal_offset.append(float(np.mean(resid[m])))
    gal_npts.append(int(m.sum()))
gal_offset = np.array(gal_offset)
print(f'per-galaxy offsets: {len(gal_names)} galaxies, '
      f'offset std = {np.std(gal_offset):.3f} dex')

# --------------------------------------------------------------- environment
cf4 = parse_cf4()
cf4_xyz = unit_vec(cf4['ra'], cf4['dec']) * cf4['dist'][:, None]
cf4_uv = unit_vec(cf4['ra'], cf4['dec'])

env = {}
for name in gal_names:
    c = coords[name]
    D = t1[name]['D']
    ra, dec = np.radians(c['ra']), np.radians(c['dec'])
    uv = unit_vec(np.array(ra), np.array(dec))
    xyz = uv * D
    vsys = H0 * D

    # angular separation -> projected distance at the SPARC galaxy's D
    cosang = np.clip(cf4_uv @ uv, -1, 1)
    ang = np.arccos(cosang)
    proj = ang * D                        # Mpc (small-angle; fine at these scales)
    dv = np.abs(cf4['vcmb'] - vsys)

    selfm = (proj < 0.1) & (dv < 200)     # exclude the galaxy itself
    cyl2 = int(np.sum((proj < 2.0) & (dv < 500) & ~selfm))

    d3 = np.linalg.norm(cf4_xyz - xyz, axis=1)
    sph5 = int(np.sum((d3 < 5.0) & ~selfm))
    d3s = np.sort(d3[~selfm])
    d5nn = float(d3s[4])                  # distance to 5th NN, Mpc
    rho5 = 5.0 / (4.0 / 3.0 * np.pi * d5nn ** 3)

    # Virgo-centric distance (12h30m, +12d18m, 16.5 Mpc)
    virgo = unit_vec(np.radians(187.7), np.radians(12.34)) * 16.5
    dvirgo = float(np.linalg.norm(xyz - virgo))

    env[name] = dict(cyl2=cyl2, sph5=sph5, rho5=rho5, dvirgo=dvirgo, D=D)

# ------------------------------------------------------ correlate (primary)
D_arr = np.array([env[n]['D'] for n in gal_names])
logD = np.log10(D_arr)

def dist_correct(metric):
    """residual of metric after linear regression on log D (CF4 is
    flux-limited: raw counts fall with distance; this removes that)."""
    A = np.stack([np.ones_like(logD), logD], axis=1)
    coef, *_ = np.linalg.lstsq(A, metric, rcond=None)
    return metric - A @ coef

log_cyl2 = np.log10(1 + np.array([env[n]['cyl2'] for n in gal_names]))
log_sph5 = np.log10(1 + np.array([env[n]['sph5'] for n in gal_names]))
log_rho5 = np.log10(np.array([env[n]['rho5'] for n in gal_names]))
dvirgo = np.array([env[n]['dvirgo'] for n in gal_names])

prim_env = dist_correct(log_cyl2)

r_p, p_p = stats.pearsonr(gal_offset, prim_env)
rs_p, ps_p = stats.spearmanr(gal_offset, prim_env)
print('\n================ PRIMARY (pre-declared) ================')
print(f'RAR offset vs dist-corrected log(1+N_cyl[2Mpc,500km/s]):')
print(f'  Pearson  r = {r_p:+.3f}  (r^2 = {r_p**2:.4f})  p = {p_p:.3g}')
print(f'  Spearman r = {rs_p:+.3f}  p = {ps_p:.3g}')
print(f'  N = {len(gal_names)} galaxies')
verdict = ('REFUTED (r^2 < 0.09)' if r_p ** 2 < 0.09 else
           'EXPECTATION MET (r^2 >= 0.20)' if r_p ** 2 >= 0.20 else
           'ALIVE BUT UNDER-STRENGTH (0.09 <= r^2 < 0.20)')
print(f'  REGISTERED VERDICT: {verdict}')

print('\n---------------- secondaries (context) ----------------')
for label, metric in [
        ('log(1+N_sphere[5Mpc])  dist-corr', dist_correct(log_sph5)),
        ('log rho_5NN (3D)       dist-corr', dist_correct(log_rho5)),
        ('log(1+N_cyl) RAW (no dist corr) ', log_cyl2),
        ('Virgo-centric distance [Mpc]    ', dvirgo),
        ('log D (distance itself)         ', logD)]:
    r, p = stats.pearsonr(gal_offset, metric)
    rs, _ = stats.spearmanr(gal_offset, metric)
    print(f'  {label}: r={r:+.3f} r^2={r**2:.4f} p={p:.3g} (spearman {rs:+.3f})')

# tercile classification: cluster / field / void by primary density
terc = np.quantile(prim_env, [1 / 3, 2 / 3])
lo = gal_offset[prim_env <= terc[0]]
mid = gal_offset[(prim_env > terc[0]) & (prim_env <= terc[1])]
hi = gal_offset[prim_env > terc[1]]
F, pF = stats.f_oneway(lo, mid, hi)
print(f'\n  terciles (void/field/cluster): '
      f'{np.mean(lo):+.4f} / {np.mean(mid):+.4f} / {np.mean(hi):+.4f} dex '
      f'(ANOVA F={F:.2f}, p={pF:.3g})')

# partial correlation: control internal structure (log L36, log SBeff, log D)
logL = np.log10(np.array([max(t1[n]['L36'], 1e-4) for n in gal_names]))
logSB = np.log10(np.array([max(t1[n]['SBeff'], 1e-2) for n in gal_names]))

def partial_r(x, y, controls):
    A = np.stack([np.ones_like(x)] + controls, axis=1)
    rx = x - A @ np.linalg.lstsq(A, x, rcond=None)[0]
    ry = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    return stats.pearsonr(rx, ry)

rp, pp = partial_r(gal_offset, prim_env, [logL, logSB, logD])
print(f'  partial r (control logL36, logSBeff, logD): '
      f'r={rp:+.3f} r^2={rp**2:.4f} p={pp:.3g}')

# save per-galaxy table for the record
out = {n: dict(offset_dex=float(gal_offset[i]), n_pts=gal_npts[i], **env[n])
       for i, n in enumerate(gal_names)}
outfile = os.path.join(HERE, 'test08_per_galaxy_results.json')
json.dump(out, open(outfile, 'w'), indent=1)
print(f'\nper-galaxy table -> {outfile}')
