#!/usr/bin/env python3
# Publisher 2026-09-02 -- IN-REPO REPRODUCTION of synchronism-site explorer/findings/scripts/eps0_mass_relation_last_escape.py
# (site explorer session 2026-09-01). At scan time the site's committed output (280d92c, 08:16:00 PDT) stopped at eps0 = 0.3551
# of the sweep; the full R0-R5 output finished at 08:17:39 and sat UNCOMMITTED in the site working tree, with no finding
# inscribed. "An uncommitted output is not a result": this file is the site script verbatim except for three path lines
# (the solver modules and the 08-30 E2 cache are read from the site checkout; the .npy is written here) and runs on the
# in-repo real SPARC through the site loader. Printed values are compared line by line against the site's working-tree
# output (publisher_20260902_last_escape_compare.txt). Seeded rng (20260901) => permutation p and bootstrap CIs should match exactly.
r"""
IS eps0(M_bar) A RELATION OF THE THEORY, OR IS IT WHAT MOND INDUCES?
====================================================================
Explorer session 2026-09-01.  Topic: eps0-mass-relation-the-last-escape.md

2026-08-30 measured the boost ceiling per galaxy under the framework's own field
equation  div(C(rho) grad Phi) = 4 pi G rho  and found the fitted eps0 tracks
baryonic mass (rho_s = +0.758 with rho_c FROZEN; +0.512 with rho_c co-fit),
while MOND's a0 measured identically does not (+0.07).  The seeded question is
whether eps0(M_bar) is tight enough to be an extra stated relation of the theory.

WAKE reframed it.  If MOND is empirically right, the outer boost is nu(g_bar/a0),
so ANY per-galaxy amplitude fit must come out as a function of the galaxy's outer
acceleration -- which is a function of (M_bar, R).  A mass correlation is what
MOND PREDICTS a density-keyed amplitude fit will show.  So the decisive test is
not "is eps0(M) tight?" but "is eps0 anything other than what MOND would induce?"

That is answerable at zero extra cost: the PDE solve does not depend on the
target, so every solve is scored twice -- against the DATA and against MOND's
PREDICTED curve for the same galaxy.  eps0_MOND,i is then "the eps0 this galaxy
would measure if MOND were exactly true".

PRE-REGISTERED (logs/2026-09-01.md, written before this script ran):
  R0  uncensored eps0 grid 0.005..0.98, rho_c fixed (A) and co-fit (B, from the
      08-30 E2 grid on disk).
  R1  OLS log eps0 = a + k log(M/1e10); residual sigma.  "relation" iff sigma <
      sigma(log a0 per galaxy) [pipeline noise floor for a near-universal
      constant]; "tight" iff sigma < 0.20 dex.
  R2  same for a0.  If a0 gains a comparable |k'|, R1 is pipeline artifact.
  R3  parameter-matched global refit: class with eps0 = A M^k (2 params, rho_c
      fixed) vs MOND with a0 = A' M^k' (2 params), via per-galaxy chi2(param)
      interpolation, verified by real solves.  Supermodel-with-real-difference
      iff class chi2/N <= 1.2 x MOND's at matched count.
  R4  eps0_MOND,i vs eps0_data,i.  If MOND predicts eps0_data better than M_bar
      does, the relation is MOND-induced.  Partial regression decides whether
      M_bar carries information beyond eps0_MOND.
  R5  R1 residuals vs g_bar,last, R_last/R_d, Sigma_eff: what the relation proxies.
  null  permutation over galaxy labels for every rho_s.

NUISANCES: Upsilon_disk = 0.5 FIXED both sides; Upsilon_bul = 0.7; distance and
inclination NOT marginalised (distance moves M_bar ~ D^2 and eps0_fit ~ D
together; R2 is the guard).  chi2 weights galaxies by point count.
"""
import os
import sys
import time
import numpy as np
from scipy.stats import spearmanr, wilcoxon

SITE = "/mnt/c/exe/projects/ai-agents/synchronism-site/explorer/findings/scripts"  # site solver + 08-30 E2 cache, imported not copied
HERE = SITE
OUT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l2_sparc_core as K                      # noqa: E402
import l2_field_equation_on_sparc as D         # noqa: E402

RHO_C_FIXED = 3.5e-06
GAMMA = 0.489
EPS = np.logspace(np.log10(0.005), np.log10(0.98), 27)
A0S = np.logspace(np.log10(0.03e-10), np.log10(30e-10), 61)
A0_LIT = 1.2e-10
rng = np.random.default_rng(20260901)
SMOKE = "--smoke" in sys.argv


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78, flush=True)


def score_to(g, v2pred, target):
    v = np.sqrt(np.clip(v2pred, 0, None))
    o = g.ok & (v > 0)
    r = (target[o] - v[o]) / g.d["eVobs"][o]
    return float(np.sum(r**2))


def refine(grid_log, chi):
    """Sub-grid minimum by a parabola through the 3 points around the argmin
    (in log parameter).  Returns (log_param_min, at_edge)."""
    j = int(np.argmin(chi))
    if j == 0 or j == len(chi) - 1:
        return grid_log[j], True
    x0, x1, x2 = grid_log[j - 1:j + 2]
    y0, y1, y2 = chi[j - 1:j + 2]
    denom = (y0 - 2 * y1 + y2)
    if denom <= 0:
        return grid_log[j], False
    dx = 0.5 * (y0 - y2) / denom * (x1 - x0)
    return x1 + np.clip(dx, -(x1 - x0), (x2 - x1)), False


def ols(x, y):
    """y = a + k x.  Returns a, k, residual std (ddof=2), residuals."""
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    return coef[0], coef[1], float(np.sqrt(np.sum(res**2) / max(len(y) - 2, 1))), res


def robust_sigma(res):
    return float(1.4826 * np.median(np.abs(res - np.median(res))))


def perm_p(x, y, n=20000):
    r0 = spearmanr(x, y)[0]
    cnt = 0
    for _ in range(n):
        if abs(spearmanr(x, rng.permutation(y))[0]) >= abs(r0):
            cnt += 1
    return r0, (cnt + 1) / (n + 1)


t0 = time.time()
gals = K.load_sparc()
G = [D.Gal(d) for d in gals.values()]
if SMOKE:
    G = G[:8]; EPS = EPS[::6]; A0S = A0S[::10]
NG = len(G)
NN = np.array([int(g.ok.sum()) for g in G])
print(f"built {NG} galaxies in {time.time()-t0:.0f}s;  N = {NN.sum()} points", flush=True)

# MOND's predicted curve for each galaxy (simple mu, a0 = 1.2e-10), as pseudo-data
V_MOND = [np.sqrt(K.mond_simple(g.gbar_sparc) * g.d["R"]) for g in G]

# galaxy observables
logM = np.array([np.log10(g.Mtot) for g in G]) - 10.0          # log(M_bar / 1e10)
logVf = np.array([np.log10(max(g.d["props"]["Vflat"], 10.0)) for g in G])
logRd = np.array([np.log10(max(g.d["props"]["Rdisk"], 0.05)) for g in G])
Rlast = np.array([g.d["R"][g.ok][-1] for g in G])
g_last = np.array([g.gbar_sparc[g.ok][-1] for g in G])          # (km/s)^2/kpc
B_last = np.array([g.d["Vobs"][g.ok][-1]**2 / g.vbar2[g.ok][-1] for g in G])
Sig_eff = np.array([g.Mtot / (2 * np.pi * (1e3 * max(g.d["props"]["Rdisk"], 0.05))**2) for g in G])  # Msun/pc^2
rho_mid = np.array([np.median(g.rho_mid[g.ok] / K.KPC3) for g in G])
y_last = g_last / K.A0_KPC
mu_last = 1.0 / (0.5 + np.sqrt(0.25 + 1.0 / y_last))            # 1/nu_simple(y_last) = mu at the last point

# ------------------------------------------------------------------ solves
hdr("0. SOLVES: eps0 grid, rho_c fixed; every solve scored against DATA and against MOND")
chi_data = np.zeros((len(EPS), NG))
chi_mond = np.zeros((len(EPS), NG))
for a, e0 in enumerate(EPS):
    Cf = K.C_framework(GAMMA, RHO_C_FIXED, e0)
    for i, g in enumerate(G):
        v2 = g.solve(Cf, e0)[0] * g.vbar2
        chi_data[a, i] = score_to(g, v2, g.d["Vobs"])
        chi_mond[a, i] = score_to(g, v2, V_MOND[i])
    print(f"   eps0 = {e0:.4f}  chi2/N vs data {chi_data[a].sum()/NN.sum():8.2f}"
          f"   vs MOND-curve {chi_mond[a].sum()/NN.sum():8.2f}   [{time.time()-t0:.0f}s]", flush=True)

chi_a0 = np.zeros((len(A0S), NG))
old = K.A0_KPC
for a, a0 in enumerate(A0S):
    K.A0_KPC = a0 / 3.24078e-14
    for i, g in enumerate(G):
        chi_a0[a, i] = score_to(g, K.mond_simple(g.gbar_sparc) * g.d["R"], g.d["Vobs"])
K.A0_KPC = old

lE, lA = np.log10(EPS), np.log10(A0S)
e_data = np.zeros(NG); e_mond = np.zeros(NG); a_gal = np.zeros(NG)
edge_e = np.zeros(NG, bool); edge_m = np.zeros(NG, bool); edge_a = np.zeros(NG, bool)
for i in range(NG):
    e_data[i], edge_e[i] = refine(lE, chi_data[:, i])
    e_mond[i], edge_m[i] = refine(lE, chi_mond[:, i])
    a_gal[i], edge_a[i] = refine(lA, chi_a0[:, i])

# treatment B: co-fit (eps0, rho_c) per galaxy from the 08-30 E2 grid
try:
    co = np.load(os.path.join(HERE, "epsilon0_per_galaxy_fw.npy"))
    e_cofit = np.log10(co[0][:NG]); rc_cofit = co[1][:NG]
    have_B = (len(co[0]) == NG)
except Exception:
    have_B = False

# ------------------------------------------------------------------ R0
hdr("R0. UNCENSORED PER-GALAXY CONSTANTS")
print(f"   grid: eps0 {EPS[0]:.4f}..{EPS[-1]:.3f} ({len(EPS)} pts);  a0 {A0S[0]:.2e}..{A0S[-1]:.2e} ({len(A0S)} pts)")
for nm, v, ed in (("eps0 (data)", e_data, edge_e), ("eps0 (MOND-curve)", e_mond, edge_m), ("a0", a_gal, edge_a)):
    p16, p50, p84 = np.percentile(v, [16, 50, 84])
    print(f"   {nm:<18s} median {10**p50:9.3e}   16-84% {p84-p16:5.3f} dex   std {np.std(v):5.3f} dex   at edge {100*ed.mean():4.0f}%")
if have_B:
    p16, p50, p84 = np.percentile(e_cofit, [16, 50, 84])
    print(f"   {'eps0 (co-fit B)':<18s} median {10**p50:9.3e}   16-84% {p84-p16:5.3f} dex   std {np.std(e_cofit):5.3f} dex"
          f"   rho_c at top edge {100*np.mean(rc_cofit >= 30):3.0f}%")
ie = int(chi_data.sum(1).argmin()); ia = int(chi_a0.sum(1).argmin())
print(f"\n   universal: class eps0 = {EPS[ie]:.3f} chi2/N = {chi_data[ie].sum()/NN.sum():.2f};"
      f"  MOND a0 = {A0S[ia]:.3e} chi2/N = {chi_a0[ia].sum()/NN.sum():.2f}")
print(f"   per-galaxy 1-param: class chi2/N = {chi_data.min(0).sum()/NN.sum():.2f};  MOND {chi_a0.min(0).sum()/NN.sum():.2f}")
print(f"   sanity: eps0 fitted to MOND's own curve, universal best = {EPS[int(chi_mond.sum(1).argmin())]:.3f}"
      f"  (chi2/N {chi_mond.sum(1).min()/NN.sum():.2f}: how well the class can imitate MOND at all)")

# ------------------------------------------------------------------ R1 / R2
hdr("R1/R2. THE MASS RELATION -- eps0 vs a0, same galaxies, same pipeline")
SIG_A0 = float(np.std(a_gal))
print(f"   decision thresholds (pre-registered): 'relation' iff resid sigma < sigma(log a0) = {SIG_A0:.3f} dex;"
      f"  'tight' iff < 0.200 dex\n")
print(f"   {'constant':<18s} {'k (dex/dex)':>12s} {'a @1e10':>9s} {'resid sigma':>12s} {'robust':>8s} {'raw sigma':>10s} {'rho_s':>7s} {'perm p':>9s}")
fits = {}
rows = [("eps0 (data)", e_data), ("eps0 (MOND-curve)", e_mond), ("a0", a_gal)]
if have_B:
    rows.append(("eps0 (co-fit B)", e_cofit))
for nm, v in rows:
    a, k, s, res = ols(logM, v)
    r, p = perm_p(logM, v, 5000 if SMOKE else 20000)
    fits[nm] = (a, k, s, res)
    print(f"   {nm:<18s} {k:+12.3f} {10**a:9.3e} {s:12.3f} {robust_sigma(res):8.3f} {np.std(v):10.3f} {r:+7.3f} {p:9.2e}")
print("\n   uncensored only (drop galaxies whose eps0 sits at a grid edge):")
m = ~edge_e
a, k, s, res = ols(logM[m], e_data[m])
print(f"   eps0 (data), n={m.sum()}: k = {k:+.3f}, resid sigma = {s:.3f} dex, robust {robust_sigma(res):.3f}")

# ------------------------------------------------------------------ R4
hdr("R4. IS eps0_data ANYTHING OTHER THAN WHAT MOND INDUCES?")
d = e_data - e_mond
print(f"   log(eps0_data / eps0_MOND): median {np.median(d):+.3f}, sigma {np.std(d):.3f} dex, robust {robust_sigma(d):.3f} dex")
r, p = perm_p(e_mond, e_data, 5000 if SMOKE else 20000)
print(f"   Spearman(eps0_data, eps0_MOND) = {r:+.3f}  (perm p = {p:.2e})")
a, k, s, res = ols(e_mond, e_data)
print(f"   OLS eps0_data on eps0_MOND: slope {k:+.3f}, resid sigma {s:.3f} dex  (vs {fits['eps0 (data)'][2]:.3f} dex about M_bar)")
r2, p2 = perm_p(logM, e_mond, 5000 if SMOKE else 20000)
print(f"   Spearman(eps0_MOND, log M) = {r2:+.3f} (perm p = {p2:.2e})  <- the mass relation MOND *induces* in this pipeline")
# partial regression: does M_bar add information beyond eps0_MOND, and vice versa?
X = np.vstack([np.ones(NG), e_mond, logM]).T
coef, *_ = np.linalg.lstsq(X, e_data, rcond=None)
res2 = e_data - X @ coef
s2 = np.sqrt(np.sum(res2**2) / (NG - 3))
# bootstrap CIs
bs = []
for _ in range(2000):
    idx = rng.integers(0, NG, NG)
    c_, *_ = np.linalg.lstsq(X[idx], e_data[idx], rcond=None)
    bs.append(c_)
bs = np.array(bs)
print(f"\n   joint OLS  log eps0_data = a + b log eps0_MOND + c log M :")
print(f"      b = {coef[1]:+.3f}  [95% {np.percentile(bs[:,1],2.5):+.3f}, {np.percentile(bs[:,1],97.5):+.3f}]")
print(f"      c = {coef[2]:+.3f}  [95% {np.percentile(bs[:,2],2.5):+.3f}, {np.percentile(bs[:,2],97.5):+.3f}]   (M_bar's slope AFTER MOND's induced eps0 is removed)")
print(f"      resid sigma {s2:.3f} dex   (M only: {fits['eps0 (data)'][2]:.3f};  eps0_MOND only: {s:.3f})")
# the analytic version: eps0 ~ mu(y_last)?
r3 = spearmanr(np.log10(mu_last), e_data)[0]
r4 = spearmanr(np.log10(mu_last), e_mond)[0]
print(f"\n   analytic proxy mu_simple(g_bar,last/a0): Spearman with eps0_data {r3:+.3f}, with eps0_MOND {r4:+.3f}")

# ------------------------------------------------------------------ R5
hdr("R5. WHAT THE RESIDUALS ABOUT eps0(M) STILL KNOW")
res_e = fits["eps0 (data)"][3]; res_a = fits["a0"][3]
print(f"   {'observable':<16s} {'rho_s(res eps0|M)':>18s} {'rho_s(res a0|M)':>16s}   {'rho_s(raw eps0)':>15s} {'rho_s(raw a0)':>13s}")
for nm, v in (("log g_bar,last", np.log10(g_last)), ("log B_req,last", np.log10(B_last)),
              ("log R_last/R_d", np.log10(Rlast) - logRd), ("log Sigma_eff", np.log10(Sig_eff)),
              ("log rho_mid", np.log10(rho_mid)), ("log V_flat", logVf), ("log R_d", logRd),
              ("log eps0_MOND", e_mond)):
    print(f"   {nm:<16s} {spearmanr(v, res_e)[0]:+18.3f} {spearmanr(v, res_a)[0]:+16.3f}   "
          f"{spearmanr(v, e_data)[0]:+15.3f} {spearmanr(v, a_gal)[0]:+13.3f}")

# ------------------------------------------------------------------ R3
hdr("R3. PARAMETER-MATCHED GLOBAL FIT: class with eps0 = A (M/1e10)^k  vs  MOND with a0 = A' (M/1e10)^k'")


def interp_chi(lgrid, chi, lp):
    """Per-galaxy chi2 at log-param lp by linear interpolation on the grid (clamped)."""
    lp = np.clip(lp, lgrid[0], lgrid[-1])
    tot = 0.0
    for i in range(chi.shape[1]):
        tot += np.interp(lp[i], lgrid, chi[:, i])
    return tot


KS = np.linspace(-0.6, 1.4, 81)
best = {}
for nm, lgrid, chi in (("class eps0", lE, chi_data), ("MOND a0", lA, chi_a0)):
    AS = np.linspace(lgrid[0], lgrid[-1], 71)
    tab = np.zeros((len(AS), len(KS)))
    for ia_, A in enumerate(AS):
        for ik, k in enumerate(KS):
            tab[ia_, ik] = interp_chi(lgrid, chi, A + k * logM)
    ia_, ik = np.unravel_index(tab.argmin(), tab.shape)
    k0 = int(np.argmin(np.abs(KS)))
    best[nm] = (AS[ia_], KS[ik], tab[ia_, ik] / NN.sum(), tab[:, k0].min() / NN.sum())
    print(f"   {nm:<10s}: best A = {10**AS[ia_]:.3e}, k = {KS[ik]:+.3f}  ->  chi2/N = {tab[ia_, ik]/NN.sum():8.2f}"
          f"   (k = 0 universal: {tab[:, k0].min()/NN.sum():8.2f})")
ratio = best["class eps0"][2] / best["MOND a0"][2]
print(f"\n   matched 2-parameter ratio class/MOND = {ratio:.2f}x   (pre-registered pass: <= 1.20)")
print(f"   what the relation buys each side: class {best['class eps0'][3]/best['class eps0'][2]:.2f}x,"
      f"  MOND {best['MOND a0'][3]/best['MOND a0'][2]:.2f}x")

# verify the class relation with real solves (interpolation check) and get B-profiles
print("\n   verification with real solves at the best class relation:")
A, k = best["class eps0"][0], best["class eps0"][1]
tot = 0.0; shape = []; shape_m = []
for i, g in enumerate(G):
    e0 = float(np.clip(10**(A + k * logM[i]), EPS[0], EPS[-1]))
    B, _ = g.solve(K.C_framework(GAMMA, RHO_C_FIXED, e0), e0)
    tot += score_to(g, B * g.vbar2, g.d["Vobs"])
    ok = g.ok
    Bok = B[ok]; nu = (g.d["R"][ok] * K.mond_simple(g.gbar_sparc[ok])) / g.vbar2[ok]
    shape.append(np.log10(Bok[-1] / Bok[len(Bok) // 2])); shape_m.append(np.log10(nu[-1] / nu[len(nu) // 2]))
print(f"   real-solve chi2/N = {tot/NN.sum():8.2f}  (interpolated {best['class eps0'][2]:8.2f})")
print(f"   boost SHAPE at the relation: median log(B_last/B_mid) class {np.median(shape):+.3f} dex vs MOND nu {np.median(shape_m):+.3f} dex"
      f"  -- with rho_c = {RHO_C_FIXED:.1e} the class's boost is {'nearly uniform' if abs(np.median(shape)) < 0.05 else 'radially varying'} in the disc")

# MOND at the best relation, exact
A2, k2 = best["MOND a0"][0], best["MOND a0"][1]
tot2 = 0.0
for i, g in enumerate(G):
    K.A0_KPC = 10**(A2 + k2 * logM[i]) / 3.24078e-14
    tot2 += score_to(g, K.mond_simple(g.gbar_sparc) * g.d["R"], g.d["Vobs"])
K.A0_KPC = old
print(f"   MOND at a0 = A' M^k' exact: chi2/N = {tot2/NN.sum():8.2f}  (interpolated {best['MOND a0'][2]:8.2f})")

# ------------------------------------------------------------------ galaxy-level
hdr("GALAXY-LEVEL (N_eff = number of galaxies) for the matched 2-parameter fits")
pg_c = np.array([np.interp(np.clip(A + k * logM[i], lE[0], lE[-1]), lE, chi_data[:, i]) for i in range(NG)]) / NN
pg_m = np.array([np.interp(np.clip(A2 + k2 * logM[i], lA[0], lA[-1]), lA, chi_a0[:, i]) for i in range(NG)]) / NN
print(f"   median reduced chi2: class {np.median(pg_c):.2f}  MOND {np.median(pg_m):.2f};  class wins {100*np.mean(pg_c < pg_m):.0f}%;"
      f"  Wilcoxon p = {wilcoxon(pg_c, pg_m).pvalue:.2e}")

np.save(os.path.join(OUT, "publisher_20260902_last_escape_reproduction.npy"),
        np.vstack([e_data, e_mond, a_gal, logM, logVf, logRd, np.log10(g_last), np.log10(B_last),
                   np.log10(Rlast), np.log10(Sig_eff), edge_e, edge_m, edge_a, NN]))
print(f"\ntotal {time.time()-t0:.0f}s")
