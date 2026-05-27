#!/usr/bin/env python3
"""
Session 675: Provenance check of TEST-17 (scale-dependent speed of light).
Continues the frontier per-test verification S674 recommended; applies S632's
method (does the stated number follow from the stated derivation?) by COMPUTATION.

Stated formula (Literature_Review_Novelty_2025-11-06.md line 84):
    c_eff(kappa) = c_0 * (1 + alpha * ln(kappa / ell_P)),   alpha ~ 1e-5

Catalog TEST-17 specific predictions (EXPERIMENTAL_TEST_CATALOG.md):
    atomic   kappa=1e-10 m  ->  -17 km/s
    solar    kappa=1e12  m  ->  +33 km/s
    galactic kappa=1e20  m  ->  +39 km/s
    headline: ~60 km/s difference between atomic-derived c and astrophysical c

Question: do -17/+33/+39 follow from the formula? Two checks:
  (1) plug kappa into the formula -> compare sign & magnitude
  (2) collinearity: ANY law c_eff = c_0(1+alpha ln(kappa/ref)) makes Delta_c LINEAR
      in ln(kappa). Are the three catalog points collinear in ln(kappa)? A log law
      forces a CONSTANT slope (km/s per e-fold). Check it.
"""
import numpy as np

c0 = 299792.458          # km/s
ellP = 1.616e-35         # m
alpha = 1e-5             # stated

kappa = {"atomic": 1e-10, "solar": 1e12, "galactic": 1e20}
catalog = {"atomic": -17.0, "solar": +33.0, "galactic": +39.0}

print("=" * 74)
print("(1) Does the stated formula reproduce the catalog numbers?")
print("=" * 74)
print(f"   c_eff(k) = c0*(1 + alpha*ln(k/ellP)), c0={c0:.0f} km/s, alpha={alpha:g}")
print(f"   {'scale':<9} {'kappa(m)':>9} {'ln(k/ellP)':>11} {'formula dC (km/s)':>18} {'catalog':>9}")
for s,k in kappa.items():
    L = np.log(k/ellP)
    dC = c0 * alpha * L
    print(f"   {s:<9} {k:>9.0e} {L:>11.1f} {dC:>18.1f} {catalog[s]:>9.0f}")
print("   => formula gives all POSITIVE (k>>ellP => ln>0) and ~170-390 km/s.")
print("      Catalog has a NEGATIVE atomic value and magnitudes ~10x smaller.")
print("      The numbers do NOT come from the stated formula with alpha=1e-5.")

print()
print("=" * 74)
print("(2) Are the three catalog points collinear in ln(kappa)? (a log law requires it)")
print("=" * 74)
lnk = np.array([np.log(kappa[s]) for s in ["atomic","solar","galactic"]])
dc  = np.array([catalog[s]       for s in ["atomic","solar","galactic"]])
# slope from atomic->solar, then predict galactic; compare to catalog
slope_as = (dc[1]-dc[0])/(lnk[1]-lnk[0])
pred_gal = dc[0] + slope_as*(lnk[2]-lnk[0])
slope_sg = (dc[2]-dc[1])/(lnk[2]-lnk[1])
print(f"   points (ln k, dC): "
      f"({lnk[0]:.1f},{dc[0]:.0f}) ({lnk[1]:.1f},{dc[1]:.0f}) ({lnk[2]:.1f},{dc[2]:.0f})")
print(f"   slope atomic->solar   = {slope_as:.3f} km/s per e-fold")
print(f"   slope solar->galactic = {slope_sg:.3f} km/s per e-fold")
print(f"   a single log law needs ONE constant slope; these differ by "
      f"{slope_as/slope_sg:.1f}x.")
print(f"   line through (atomic,solar) predicts galactic = {pred_gal:.1f} km/s; "
      f"catalog says {dc[2]:.0f}.")
print("   => the three points are NOT collinear in ln(kappa). No single logarithmic")
print("      c_eff(kappa) law produces -17/+33/+39. The numbers are picked, not derived.")

print()
print("=" * 74)
print("(3) Conceptual + empirical problems (independent of the above)")
print("=" * 74)
dc_over_c = 60.0/c0
print(f"   claimed atomic-vs-astrophysical difference ~60 km/s => Dc/c ~ {dc_over_c:.1e}")
print("   - c has been SI-DEFINED exactly since 1983 (299792458 m/s); it is not")
print("     'measured' at a scale. The 'atomic c vs astrophysical c' framing is")
print("     confused: what is testable is consistency of c-dependent relations, and")
print("     those agree far better than 2e-4.")
print("   - A scale/position-dependent c at the 2e-4 level violates local Lorentz")
print("     invariance, constrained to <~1e-15 (clock-comparison, Michelson-Morley")
print("     class). Excluded by ~11 orders of magnitude.")
print("   - Contradicts the framework's OWN substrate analysis: S667 (continuum")
print("     transfer rule is parabolic => infinite propagation speed, no finite c to")
print("     vary) and S641 (Lorentz invariance is an unfilled kinematic-layer gap).")

print()
print("=" * 74)
print("VERDICT: TEST-17 amplitude is NOT derived (same failure mode as S632's 500 Mpc)")
print("=" * 74)
print("   - numbers do not follow from the stated formula (wrong sign, ~10x magnitude)")
print("   - not collinear in ln(kappa): no single log law yields them")
print("   - claimed 2e-4 variation is excluded by Lorentz tests by ~11 OOM, and the")
print("     'scale-measured c' framing conflicts with c being SI-defined")
print("   - contradicts S667 (parabolic substrate) and S641 (Lorentz gap)")
print("   Frontier provenance now settled for TEST-07 (S632), TEST-12 (self-flagged),")
print("   TEST-17 (here): all NOT derived. 6 of 9 frontier tests still genuinely unverified.")
