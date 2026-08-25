#!/usr/bin/env python3
"""
Independent Publisher-lane verification of the synchronism-site explorer finding
"the argument of C: three functions" (2026-08-24).

Four claims are checked, each re-implemented from scratch (the site's scripts were
NOT imported or read for the numerics):

  R1  deep-limit slopes: explicit-in-g gives dlog g_obs/dlog g_bar -> 0 (v ~ sqrt(r),
      never flat); implicit-in-g gives 1/2 (MOND's flat law).
  R2  at gamma = 1/2, tanh(0.5*ln(1+x)) == x/(x+2) == mu_simple(x/2), and the
      algebraic-mu EFE response factor mu(x)[1+dln mu/dln x] coincides with
      MOND-simple's at x/2.
  R3  head-to-head fit on real SPARC: C keyed to g_obs (implicit) vs C keyed to
      local baryonic density rho. Same points, same likelihood, 3 free parameters
      each, so BIC penalties cancel and dBIC = d(-2 lnL).
  R4  galaxy-level bootstrap on gamma (the site claims the published +/-0.02 is
      correlated-N inflation).

Data: simulations/sparc_real_data/MassModels_Lelli2016c.mrt (real SPARC, Lelli+2016)
      simulations/sparc_real_data/SPARC_Lelli2016c.mrt      (Rdisk, Vflat, Q, Inc)
"""
import numpy as np
from scipy.optimize import brentq, minimize
import sympy as sp

rng = np.random.default_rng(20260825)

KPC   = 3.0856775814913673e19   # m
KMS   = 1.0e3                   # m/s
MSUN  = 1.98892e30              # kg
G     = 6.67430e-11
LSUN_PC2_TO_MSUN_M2 = None      # set per upsilon below

# ----------------------------------------------------------------- data load
def _data_lines(path):
    """Return the lines after the LAST '---' separator of an AAS .mrt header."""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    last = max(i for i, l in enumerate(lines) if l.startswith("-----"))
    return lines[last + 1:]

# NOTE: the committed .mrt files do NOT honour the byte ranges in their own
# byte-by-byte headers (fields are shifted ~1-2 chars).  Verified by inspection
# 2026-08-25; parse on whitespace instead.  No field in either table contains a
# space (galaxy names and reference lists are space-free).
def load_galaxy_table(path):
    out = {}
    for ln in _data_lines(path):
        f = ln.split()
        if len(f) < 18:
            continue
        try:
            out[f[0]] = dict(inc=float(f[5]), Rdisk=float(f[11]),
                             SBdisk=float(f[12]), Vflat=float(f[15]), Q=int(f[17]))
        except ValueError:
            continue
    return out

def load_mass_models(path):
    rows = []
    for ln in _data_lines(path):
        f = ln.split()
        if len(f) < 10:
            continue
        try:
            rows.append(dict(name=f[0], R=float(f[2]), Vobs=float(f[3]), eV=float(f[4]),
                             Vgas=float(f[5]), Vdisk=float(f[6]), Vbul=float(f[7]),
                             SBd=float(f[8]), SBb=float(f[9])))
        except ValueError:
            continue
    return rows

BASE = "simulations/sparc_real_data/"
gal  = load_galaxy_table(BASE + "SPARC_Lelli2016c.mrt")
rows = load_mass_models(BASE + "MassModels_Lelli2016c.mrt")
print(f"loaded: {len(gal)} galaxies in Table1, {len(rows)} rotation-curve points")

def build_sample(ups_d=0.5, ups_b=0.7, h_mode="const0.3"):
    """Standard SPARC RAR cuts (Lelli+2016 / McGaugh+2016 style), stated explicitly:
       Q < 3, inclination > 30 deg, e_Vobs/Vobs < 0.10, R > 0, Vobs > 0, g_bar > 0."""
    name, gb, go, rho, w = [], [], [], [], []
    for r in rows:
        g = gal.get(r["name"])
        if g is None or g["Q"] >= 3 or g["inc"] <= 30:
            continue
        if r["Vobs"] <= 0 or r["R"] <= 0 or r["eV"] <= 0:
            continue
        if r["eV"] / r["Vobs"] >= 0.10:
            continue
        Rm = r["R"] * KPC
        vb2 = (r["Vgas"] * abs(r["Vgas"]) + ups_d * r["Vdisk"]**2 + ups_b * r["Vbul"]**2) * KMS**2
        if vb2 <= 0:
            continue
        gbar = vb2 / Rm
        gobs = (r["Vobs"] * KMS)**2 / Rm
        # local baryonic volume density rho = Sigma / (2h), Sigma from [3.6] surface brightness
        sigma = (ups_d * r["SBd"] + ups_b * r["SBb"]) * MSUN / (KPC/1000.0)**2   # kg/m^2
        if sigma <= 0:
            continue
        if h_mode == "const0.3":
            h = 0.3 * KPC
        elif h_mode == "Rd/5":
            h = max(g["Rdisk"], 0.1) / 5.0 * KPC
        else:
            raise ValueError(h_mode)
        name.append(r["name"]); gb.append(gbar); go.append(gobs)
        rho.append(sigma / (2.0 * h))
        w.append(2.0 * r["eV"] / (r["Vobs"] * np.log(10.0)))   # sigma on log10 g_obs
    return (np.array(name), np.array(gb), np.array(go), np.array(rho), np.array(w))

names, gbar, gobs, rho, sobs = build_sample()
print(f"sample: N = {len(gbar)} points, {len(set(names))} galaxies "
      f"(cuts: Q<3, inc>30, e_Vobs/Vobs<0.10)")

# ------------------------------------------------------- R1 deep-limit slopes
print("\n=== R1  deep-limit slopes (dlog g_obs / dlog g_bar, g_bar << a0) ===")
A0 = 1.2e-10
def C_tanhlog(x, gamma):
    return np.tanh(gamma * np.log1p(x))

def gobs_explicit(gb, a0, gamma):
    """C evaluated at the BARYONIC field: g_obs = g_bar / C(g_bar/a0)."""
    return gb / C_tanhlog(gb / a0, gamma)

def gobs_implicit(gb, a0, gamma):
    """C at the TOTAL field: solve x*tanh(g*ln(1+x)) = g_bar/a0 for x = g_obs/a0.
    Vectorised bisection in log x (f is strictly increasing in x)."""
    t = np.atleast_1d(np.asarray(gb, dtype=float) / a0)
    f = lambda x: x * np.tanh(gamma * np.log1p(x)) - t
    lo = np.maximum(t, 1e-30)                      # x >= t since tanh(.) <= 1
    hi = np.maximum(2.0 * t + 1.0, 2.0)
    for _ in range(200):                           # widen until bracketed
        bad = f(hi) < 0
        if not np.any(bad):
            break
        hi = np.where(bad, hi * 4.0, hi)
    llo, lhi = np.log(lo), np.log(hi)
    for _ in range(120):                           # 120 halvings -> machine precision
        lmid = 0.5 * (llo + lhi)
        neg = f(np.exp(lmid)) < 0
        llo = np.where(neg, lmid, llo)
        lhi = np.where(neg, lhi, lmid)
    return np.exp(0.5 * (llo + lhi)) * a0

for gamma in (0.5, 2.0):
    gb_lo = np.logspace(-13, -12, 60)
    for label, fn in (("explicit in g_bar", gobs_explicit), ("implicit in g_obs", gobs_implicit)):
        y = fn(gb_lo, A0, gamma)
        s = np.polyfit(np.log10(gb_lo), np.log10(y), 1)[0]
        print(f"  gamma={gamma:>4}  {label:<18} slope = {s:+.4f}   "
              f"deep-limit v(r) ~ r^{(1-2*s)/2:+.3f}")

# ---------------------------------------------------------- R2 gamma=1/2 EFE
print("\n=== R2  gamma = 1/2 identity, extended to the algebraic-mu EFE ===")
x = sp.symbols('x', positive=True)
lhs = sp.tanh(sp.Rational(1,2) * sp.log(1 + x))
print("  sympy simplify(tanh(0.5*ln(1+x)) - x/(x+2)) =",
      sp.simplify(sp.expand(sp.simplify(lhs)) - x/(x+2)))
# numeric over six decades
xs = np.logspace(-3, 3, 4001)
print("  max |tanh(0.5 ln(1+x)) - mu_simple(x/2)| =",
      f"{np.max(np.abs(np.tanh(0.5*np.log1p(xs)) - (xs/2)/(1+xs/2))):.3e}")
def efe_factor(mu, xs, dx=1e-6):
    m  = mu(xs)
    dl = (np.log(mu(xs*(1+dx))) - np.log(mu(xs*(1-dx)))) / (2*dx)
    return m * (1.0 + dl)
f_tanh = efe_factor(lambda u: np.tanh(0.5*np.log1p(u)), xs)
f_simp = efe_factor(lambda u: (u/2)/(1+u/2), xs)
print(f"  max |mu[1+L]_tanhlog(x) - mu[1+L]_simple(x/2)| = {np.max(np.abs(f_tanh-f_simp)):.3e}")
# and at gamma != 1/2 the same comparison, to show epsilon = 2gamma-1 is the deformation
for gamma in (0.4, 0.5093, 0.6, 2.0):
    fg = efe_factor(lambda u: np.tanh(gamma*np.log1p(u)), xs)
    print(f"  gamma={gamma:<7} eps=2g-1={2*gamma-1:+.4f}  max|diff vs simple(x/2)| = "
          f"{np.max(np.abs(fg-f_simp)):.3e}")

# ------------------------------------------------------- R3 head-to-head fit
print("\n=== R3  head-to-head fit, C(g_obs) vs C(rho), real SPARC ===")
LG_OBS = np.log10(gobs)

def nll_Cg(p, gb, lgo, s):
    lga0, gamma, lsig = p
    if not (-12 < lga0 < -8 and 0.01 < gamma < 6 and -4 < lsig < 1):
        return 1e12
    a0 = 10**lga0
    try:
        pred = np.log10(gobs_implicit(gb, a0, gamma))
    except Exception:
        return 1e12
    v = s**2 + (10**lsig)**2
    return 0.5*np.sum((lgo - pred)**2/v + np.log(2*np.pi*v))

def nll_Crho(p, gb, rh, lgo, s):
    lgrc, gamma, lsig = p
    if not (-30 < lgrc < -15 and 0.001 < gamma < 8 and -4 < lsig < 1):
        return 1e12
    rc = 10**lgrc
    C = np.tanh(gamma*np.log1p(rh/rc))
    if np.any(C <= 1e-12):
        return 1e12
    pred = np.log10(gb / C)
    v = s**2 + (10**lsig)**2
    return 0.5*np.sum((lgo - pred)**2/v + np.log(2*np.pi*v))

def best(fn, starts, args):
    bl, bp = np.inf, None
    for s0 in starts:
        r = minimize(fn, s0, args=args, method="Nelder-Mead",
                     options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-8))
        if r.fun < bl:
            bl, bp = r.fun, r.x
    return bl, bp

starts_g  = [[-10.0, 0.5, -1.0], [-10.3, 2.0, -0.9], [-9.9, 0.3, -1.1]]
starts_r  = [[-24.0, 0.5, -0.7], [-25.0, 2.0, -0.6], [-23.0, 0.1, -0.8], [-26.0, 0.05, -0.6]]

Lg, Pg = best(nll_Cg,   starts_g, (gbar, LG_OBS, sobs))
Lr, Pr = best(nll_Crho, starts_r, (gbar, rho, LG_OBS, sobs))
print(f"  C_g   : a0 = {10**Pg[0]:.3e}  gamma = {Pg[1]:.4f}  sigma_int = {10**Pg[2]:.4f} dex   -2lnL = {2*Lg:.1f}")
print(f"  C_rho : rho_c = {10**Pr[0]:.3e}  gamma = {Pr[1]:.4f}  sigma_int = {10**Pr[2]:.4f} dex   -2lnL = {2*Lr:.1f}")
dbic = 2*Lr - 2*Lg
ngal = len(set(names))
print(f"  dBIC (C_rho vs C_g) = {dbic:+.1f}   [same N={len(gbar)}, same 3 params -> BIC penalties cancel]")
print(f"  N_eff deflation (N_gal/N = {ngal}/{len(gbar)}): dBIC_eff = {dbic*ngal/len(gbar):+.1f}")

# estimator sweep
print("\n  estimator sweep (rho = Sigma/2h):")
for hm in ("const0.3", "Rd/5"):
    for ud in (0.5, 0.7):
        nm2, gb2, go2, rh2, s2 = build_sample(ups_d=ud, h_mode=hm)
        lg2 = np.log10(go2)
        l1, _ = best(nll_Cg,   starts_g, (gb2, lg2, s2))
        l2, p2 = best(nll_Crho, starts_r, (gb2, rh2, lg2, s2))
        d = 2*l2 - 2*l1
        print(f"    h={hm:<9} ups_d={ud}  N={len(gb2):<5} dBIC = {d:+9.1f}  "
              f"(deflated {d*len(set(nm2))/len(gb2):+7.1f})  gamma_rho = {p2[1]:.4f}  sig_int = {10**p2[2]:.4f}")

# ------------------------------------------------------------ R4 bootstrap
print("\n=== R4  galaxy-level bootstrap on gamma (C_g) ===")
ug = np.array(sorted(set(names)))
idx = {g: np.where(names == g)[0] for g in ug}
draws = []
for b in range(150):
    pick = rng.choice(ug, size=len(ug), replace=True)
    sel = np.concatenate([idx[g] for g in pick])
    try:
        lb, pb = best(nll_Cg, [Pg], (gbar[sel], LG_OBS[sel], sobs[sel]))
        if lb < 1e11:
            draws.append(pb[1])
    except Exception:
        pass
draws = np.array(draws)
print(f"  bootstrap draws that converged: {len(draws)}")
print(f"  gamma = {np.mean(draws):.4f} +/- {np.std(draws, ddof=1):.4f}   "
      f"(median {np.median(draws):.4f}, 16-84% [{np.percentile(draws,16):.4f}, {np.percentile(draws,84):.4f}])")
eps = 2*draws - 1
print(f"  eps = 2*gamma-1 = {np.mean(eps):+.4f} +/- {np.std(eps, ddof=1):.4f}  "
      f"-> {abs(np.mean(eps))/np.std(eps, ddof=1):.2f} sigma from 0")
