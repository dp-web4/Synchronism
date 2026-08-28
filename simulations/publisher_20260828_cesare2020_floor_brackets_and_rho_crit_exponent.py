#!/usr/bin/env python3
"""
Publisher 2026-08-28 — verify before propagating.

Three things are checked here, each against a source rather than a relay.

(A) Cesare, Diaferio, Matsakos & Angus 2020 (A&A 637 A70; arXiv:2003.07377), read in full
    today (the paper flagged UNREAD on 2026-08-27).  Their numbers, transcribed from the text:
      per-galaxy (Table 2, 30 DiskMass galaxies; Sect. 5, Fig. 8):
          eps_0 = 0.56 +/- 0.19,  Q = 0.92 +/- 0.24,  log10(rho_c / g cm^-3) = -25.30 +/- 0.70
          (mean per-galaxy uncertainties 0.16, 0.71, 1.22)
      universal combination (Sect. 5, Fig. 9, MCMC over the whole sample at fixed Upsilon, h_z):
          eps_0 = 0.661 (+0.007/-0.007),  Q = 1.79 (+0.14/-0.26),  log10 rho_c = -24.54 (+0.08/-0.07)
      Sect. 6: RG "tends to underestimate the relation (20) at low g_bar" (abstract: 0.1-0.3 dex).
    Question: where does the framework's DERIVED floor eps_0 = Omega_m = 0.315 sit relative to
    the RG literature's FITTED floors, once Cesare+2020 is added to Matsakos & Diaferio 2016?

(B) The site explorer's 2026-08-27 measurement rho_crit ∝ V^(s-2), s = dlogSigma_c/dlogV on
    SPARC Table 1 (N=129, Q<=2, V_flat>0).  Reproduced here from the in-repo copy of the table,
    with the identity p = 2 - s/2 checked on the same galaxies.

(C) The 4*pi in the whitepaper's Sect. 6.4 row "x <= 0.019 beta_J^2, no by ~40x".  The primary
    (Research/Session53_Theoretical_Foundations.md, line 62) computes A = 1/(G alpha^2 R_0^2)
    with NO 4*pi; the compilation chain (session67_B_derivation.py:372) writes A = 4*pi/(...).
    x = rho_gal/rho_crit under each convention, and the local-vs-mean density ratio of an
    exponential disc, computed rather than quoted.

Exit 0 means every assertion held.
"""
import math, os, re, sys
import numpy as np

G_PC = 4.30091e-3          # pc (km/s)^2 / Msun
MSUN_G = 1.98847e33
PC_CM = 3.0857e18
MSUN_PC3_TO_G_CM3 = MSUN_G / PC_CM**3

def dex(x): return math.log10(x)

print("=" * 88)
print("(A) Cesare+2020 read in full: the RG literature's fitted floors, all of them")
print("=" * 88)
print(f"  1 Msun/pc^3 = {MSUN_PC3_TO_G_CM3:.3e} g/cm^3")
rho_c_univ = 10 ** -24.54 / MSUN_PC3_TO_G_CM3
rho_c_gal = 10 ** -25.30 / MSUN_PC3_TO_G_CM3
rho_c_md = 10 ** -24.0 / MSUN_PC3_TO_G_CM3
site_rho = 0.161
print(f"  Cesare+2020 universal rho_c = 10^-24.54 g/cm^3 = {rho_c_univ:.2e} Msun/pc^3")
print(f"  Cesare+2020 per-galaxy mean = 10^-25.30 g/cm^3 = {rho_c_gal:.2e} Msun/pc^3 (scatter 0.70 dex)")
print(f"  M&D 2016 galaxy example     = 10^-24    g/cm^3 = {rho_c_md:.2e} Msun/pc^3")
print(f"  site 08-27 Jeans-knee median = {site_rho} Msun/pc^3 = {site_rho*MSUN_PC3_TO_G_CM3:.2e} g/cm^3"
      f" = 10^{dex(site_rho*MSUN_PC3_TO_G_CM3):.2f}")
print(f"  ratio site / Cesare-universal = {site_rho/rho_c_univ:.0f}x ({dex(site_rho/rho_c_univ):.2f} dex);"
      f" / per-galaxy mean = {site_rho/rho_c_gal:.0f}x")
print("  NOTE: different constructions (a Jeans knee at R_half vs a LOCAL permittivity parameter);"
      " the ratio is reported, not adjudicated.")

# the fitted floors
floors = [
    ("M&D 2016  NGC 6946 (spiral)", 0.25),
    ("M&D 2016  NGC 1560 (spiral)", 0.20),
    ("M&D 2016  model G0", 1/6),
    ("M&D 2016  A1991 (cluster)", 0.045),
    ("M&D 2016  A1795 (cluster)", 0.065),
    ("Cesare+2020 per-galaxy mean (30 DMS)", 0.56),
    ("Cesare+2020 per-galaxy mean -1sd", 0.56 - 0.19),
    ("Cesare+2020 per-galaxy mean +1sd", 0.56 + 0.19),
    ("Cesare+2020 UNIVERSAL fit", 0.661),
]
fw = 0.315
print(f"\n  {'system':40s} {'eps_0':>7s} {'1/eps_0':>8s}  vs framework 1/Omega_m = {1/fw:.2f}")
below, above = [], []
for name, e in floors:
    tag = "ceiling BELOW demand" if 1/e > 1/fw else "ceiling ABOVE demand"
    (below if 1/e > 1/fw else above).append(name)
    print(f"  {name:40s} {e:7.3f} {1/e:8.2f}  {tag}")
lo = min(e for _, e in floors); hi = max(e for _, e in floors)
print(f"\n  RG-literature fitted range of eps_0: {lo:.3f} .. {hi:.3f}  (spread {hi/lo:.1f}x; M&D alone was 5.6x)")
assert lo < fw < hi, "0.315 is NOT bracketed?"
print(f"  ==> 0.315 is BRACKETED by RG's own fits: {len(below)} systems demand MORE boost than 3.17,"
      f" {len(above)} (all DiskMass) demand LESS.")
print("  ==> 'the ceiling sits below RG's ENTIRE fitted range' (Publisher 08-27) was true of M&D 2016 only.")
# formal exclusion of 0.315 by the universal fit
z = (0.661 - fw) / 0.007
print(f"  Universal DMS fit excludes 0.315 at {z:.0f} sigma FORMAL (stat. only; 68% CI +/-0.007, fixed Upsilon/h_z)")
print(f"  Per-galaxy DMS distribution puts 0.315 at {(0.56-fw)/0.19:.1f} sd below the mean;"
      f" M&D spirals put it {fw/0.25:.2f}-{fw/0.20:.2f}x ABOVE their fits.")
print("  DILEMMA: if eps_0 is universal, DMS's 0.661+/-0.007 refutes 0.315; if it is not universal,"
      " a floor DERIVED from Omega_m is not the same object as RG's per-system fit.")

# optional: pull the 30 per-galaxy eps_0 from the local text extraction, if present
txt = "/tmp/cesare2020.txt"
if os.path.exists(txt):
    lines = open(txt, encoding="utf-8", errors="replace").read().splitlines()
    vals = []
    for ln in lines[820:1010]:
        toks = ln.split()
        if len(toks) == 3 and re.fullmatch(r"[0-9]\.[0-9]{2}\+[0-9.]+", toks[2]) and re.fullmatch(r"[0-9]\.[0-9]{2}", toks[1]):
            vals.append(float(toks[2].split("+")[0]))
    if vals:
        v = np.array(vals)
        print(f"\n  Table 2 eps_0 parsed from local text: N={len(v)}, min {v.min():.2f}, max {v.max():.2f},"
              f" mean {v.mean():.2f}, sd {v.std(ddof=1):.2f}  (paper: 0.56 +/- 0.19)")
        print(f"  galaxies with eps_0 <= 0.40: {(v<=0.40).sum()};  <= 0.315: {(v<=fw).sum()}")
else:
    print("\n  (local text extraction not present; per-galaxy parse skipped — paper summary values used)")

print("\n" + "=" * 88)
print("(B) rho_crit ∝ V^(s-2): reproduce the site's slopes on the in-repo SPARC Table 1")
print("=" * 88)
here = os.path.dirname(os.path.abspath(__file__))
mrt = os.path.join(here, "sparc_real_data", "SPARC_Lelli2016c.mrt")
rows = []
with open(mrt, encoding="utf-8", errors="replace") as f:
    # whitespace-separated in this copy (the byte layout in the header is not what the rows carry);
    # column ORDER is the header's: Galaxy T D e_D f_D Inc e_Inc L e_L Reff SBeff Rdisk SBdisk MHI RHI Vflat e_Vflat Q [Ref]
    for ln in f:
        t = ln.split()
        if len(t) < 18 or not t[1].isdigit() or not t[17].isdigit():
            continue
        try:
            name = t[0]; T = int(t[1]); D = float(t[2])
            Reff = float(t[9]); Rdisk = float(t[11]); SBdisk = float(t[12])
            Vflat = float(t[15]); eV = float(t[16]); Q = int(t[17])
        except ValueError:
            continue
        rows.append((name, T, D, Reff, Rdisk, SBdisk, Vflat, eV, Q))
print(f"  parsed {len(rows)} galaxies")
assert len(rows) == 175, "expected 175 SPARC galaxies"
sel = [r for r in rows if r[8] <= 2 and r[6] > 0 and r[5] > 0 and r[4] > 0]
N = len(sel)
print(f"  Q<=2, V_flat>0: N = {N}   (site: 129)")
assert N == 129
lV = np.log10([r[6] for r in sel]); lS = np.log10([r[5] for r in sel]); lR = np.log10([r[4] for r in sel])
lRe = np.log10([r[3] for r in sel]); eV = np.array([r[7] for r in sel]); V = np.array([r[6] for r in sel])

def ols(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    b, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ b; n = len(x)
    s2 = ((y - yhat) ** 2).sum() / (n - 2)
    se = math.sqrt(s2 / ((x - x.mean()) ** 2).sum())
    r = np.corrcoef(x, y)[0, 1]
    return b[0], se, r

s, s_se, r_s = ols(lV, lS)
s_inv, _, _ = ols(lS, lV); s_inv = 1 / s_inv
p, p_se, r_p = ols(lV, lR)
p_inv, _, _ = ols(lR, lV); p_inv = 1 / p_inv
pe, pe_se, _ = ols(lV, lRe)
print(f"  s = dlogSigma_c/dlogV  (OLS Sigma|V) = {s:.3f} +/- {s_se:.3f}   inverse {s_inv:.3f}   r = {r_s:.3f}"
      f"     (site: 1.837 +/- 0.176, inv 4.438, r 0.643)")
print(f"  p = dlogR_disk/dlogV   (OLS R|V)     = {p:.3f} +/- {p_se:.3f}   inverse {p_inv:.3f}   r = {r_p:.3f}"
      f"     (site: 1.038 +/- 0.098, inv 1.931, r 0.733)")
print(f"  p (R_eff)                            = {pe:.3f} +/- {pe_se:.3f}                     (site: 0.914 +/- 0.095)")
assert abs(s - 1.837) < 0.02 and abs(p - 1.038) < 0.02, "site slopes NOT reproduced"
print(f"  identity p = 2 - s/2 : predicted {2 - s/2:.3f} vs measured {p:.3f}  -> {abs(p-(2-s/2))/math.hypot(p_se, s_se/2):.2f} sigma")
# attenuation from e_Vflat
sig_meas = np.median(eV / V) / math.log(10)
sig_obs = lV.std(ddof=1)
atten = 1 - (sig_meas / sig_obs) ** 2
print(f"  attenuation: sigma(logV) meas-error {sig_meas:.4f} vs observed {sig_obs:.4f} -> factor {atten:.4f}"
      f" -> corrected s = {s/atten:.3f}  (site: 0.9930, 1.850)")
expo = s - 2; expo_se = s_se
print(f"  rho_crit exponent s-2 = {expo:+.3f} +/- {expo_se:.3f}   (site: -0.15 +/- 0.18)")
print(f"    framework's +2 excluded at {(2-expo)/expo_se:.1f} sigma;  MOND-tracking -2 excluded at {(expo+2)/expo_se:.1f} sigma")
print(f"    via p: 2-2p = {2-2*p:+.3f} +/- {2*p_se:.3f}")
assert (2 - expo) / expo_se > 10 and (expo + 2) / expo_se > 9
# the Jeans-knee density itself, S53 primitive, alpha = 1.1, R_half = 1.678 R_disk
alpha = 1.1
Rh = 1.678 * np.array([r[4] for r in sel]) * 1e3   # pc
rho_crit = V ** 2 / (G_PC * alpha ** 2 * Rh ** 2)
med = np.median(rho_crit); sc = np.std(np.log10(rho_crit), ddof=1)
print(f"  S53 primitive rho_crit = V^2/(G alpha^2 R_half^2), R_half = 1.678 R_disk, alpha=1.1:"
      f" median {med:.3f} Msun/pc^3, scatter {sc:.2f} dex   (site: 0.161, 0.45 dex)")
Rh2 = np.array([r[3] for r in sel]) * 1e3
rho2 = V ** 2 / (G_PC * alpha ** 2 * Rh2 ** 2)
print(f"    with R_half = R_eff instead: median {np.median(rho2):.3f} Msun/pc^3")
# the V^2 law the ledger uses, for comparison, at the same galaxies
print(f"    ledger law 0.029*V^2 at the sample median V = {np.median(V):.0f} km/s: {0.029*np.median(V)**2:.0f} Msun/pc^3"
      f"  ({0.029*np.median(V)**2/med:.0f}x the Jeans-knee median)")
print(f"    S53 law    0.028*V^0.5 at the same V: {0.028*np.median(V)**0.5:.3f} Msun/pc^3")

print("\n" + "=" * 88)
print("(C) the 4*pi in Sect. 6.4's 'x <= 0.019 beta_J^2 — no, by ~40x'")
print("=" * 88)
for a in (1.0, 1.1):
    x_s53 = 3 * a ** 2 / (4 * math.pi)            # A = 1/(G a^2 R^2), rho_gal = 3V^2/(4 pi G R^2)
    x_4pi = 3 * a ** 2 / (16 * math.pi ** 2)      # A = 4 pi/(G a^2 R^2)
    print(f"  alpha = {a}:  S53 form x = 3a^2/(4pi) = {x_s53:.3f} (shortfall {1/x_s53:.1f}x);"
          f"  4pi form x = 3a^2/(16pi^2) = {x_4pi:.4f} (shortfall {1/x_4pi:.0f}x);  ratio {x_s53/x_4pi:.2f} = 4pi")
assert abs(3 * 1.1 ** 2 / (16 * math.pi ** 2) - 0.023) < 0.001   # the whitepaper's 0.023 at beta_J = 1.1
# local midplane density vs mean density inside R_half, exponential disc, sech^2 vertical profile
print("  exponential disc, Sigma = Sigma_0 exp(-R/R_d), midplane rho = Sigma/(2 h_z):")
for Rd in (2.0, 3.0, 4.5):
    for hz in (0.3, Rd / 7.3):
        Rh = 1.678 * Rd
        mean = (math.pi * Rd ** 2) / (4 / 3 * math.pi * Rh ** 3)          # Sigma_0 = 1 units: M(<R_half)=pi Sigma0 Rd^2
        loc0 = 1 / (2 * hz)
        locRh = math.exp(-1.678) / (2 * hz)
        print(f"    R_d = {Rd:3.1f} kpc, h_z = {hz:.2f} kpc: local/mean at centre {loc0/mean:5.1f}x,  at R_half {locRh/mean:4.1f}x")
print("  ==> the site's 25-35x is the CENTRAL midplane ratio at h_z = 0.3 kpc; at R_half it is ~5x."
      " The 3.5x S53-form shortfall is inside the central ratio and comparable to the R_half ratio;"
      " the 40x 4pi-form shortfall is outside both. The verdict is set by the 4pi and the estimator.")

print("\nALL ASSERTIONS HELD")
