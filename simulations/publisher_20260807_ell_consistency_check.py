#!/usr/bin/env python3
"""
Publisher 2026-08-07 — re-execution check on OQ-Coarsening (whitepaper §6.4, added 2026-08-06).

Five questions, all answered by computation, none by re-reading:

  Q1. Is A a function of the PRODUCT (alpha*R0) only?  If so the 08-06 "two substitutions"
      framing and the 08-05 proposal's "317 pc" are the same number re-factorized, and the
      644x is degenerate rather than mis-substituted.
  Q2. Does whitepaper 6.4 point 3 (x prop ell^2; NGC 3198 C: 0.0003 -> 0.86) reproduce, and
      what coarse-graining convention does it assume?
  Q3. Under the two OTHER self-consistent conventions -- (a) ell smooths BOTH rho and rho_crit,
      (b) ell smooths rho ONLY, rho_crit set by galaxy size -- what does x(ell) do?
  Q4. Origin of the 1.515x density discrepancy the sibling repo flagged as unexplained
      (site/whitepaper rho(0)=0.934 vs plotter-exact 0.6164 for NGC 3198).

Run: python3 simulations/publisher_20260807_ell_consistency_check.py
"""

import numpy as np

# ---- constants (G in the archive's working units: M_sun, pc, km/s) ----
G_SI = 6.674e-11          # m^3 kg^-1 s^-2
M_sun = 1.989e30          # kg
pc = 3.086e16             # m
km = 1e3                  # m
# G in pc * (km/s)^2 / M_sun
G = G_SI * M_sun / (pc * km**2)

print("=" * 78)
print("Q1 -- is A a function of the product (alpha * R0) only?")
print("=" * 78)

def A_of(alpha, R0_kpc, fourpi=False):
    """A = [4pi] / (alpha^2 G R0^2), R0 in kpc/(km/s)^0.75 -> convert to pc units."""
    R0_pc = R0_kpc * 1e3
    return (4 * np.pi if fourpi else 1.0) / (alpha**2 * G * R0_pc**2)

A_code = A_of(4.5, 0.07)          # session66_A_gap_investigation.py inputs
A_stated = A_of(1.0, 8.0)         # the whitepaper/site's stated-formula inputs
print(f"  A(alpha=4.5, R0=0.07)  = {A_code:.5g}   <- simulations/session66_A_gap_investigation.py")
print(f"  A(beta_J=1, R0=8 kpc)  = {A_stated:.5g}   <- stated formula")
print(f"  ratio                  = {A_code / A_stated:.4f}")
print(f"  closed form (8/0.07)^2 / 4.5^2 = {(8/0.07)**2 / 4.5**2:.4f}   [Publisher 08-06]")

# product invariance: hold alpha*R0 fixed, vary the factorization
prod = 4.5 * 0.07
print(f"\n  product alpha*R0 = {prod:.4f} kpc.  Refactorizations at fixed product:")
for a in (0.5, 1.0, 1.1, 4.5, 12.0):
    print(f"    alpha={a:5.2f}  R0={prod/a:8.5f}  ->  A = {A_of(a, prod/a):.6g}")
print("  => A depends ONLY on the product. Any factorization reproduces A exactly.")

# the 08-05 proposal's inversion
ell_prop = 0.3174   # kpc, the proposal's 317 pc at beta_J = 1
print(f"\n  08-05 proposal:  beta_J=1, R0={ell_prop} kpc  -> product = {1.0*ell_prop:.4f} kpc")
print(f"  Session 66 documented: alpha=4.5, R0=0.07     -> product = {prod:.4f} kpc")
print(f"  agreement on the invariant: {abs(1.0*ell_prop - prod)/prod*100:.2f}%")
print("  => the 317 pc is NOT manufactured; it is the SAME invariant, re-factorized.")
print("     (Publisher 08-06 called it manufactured and quoted R0=70.5 pc as 'the implied")
print("      length, 4.25x off'. R0 alone is not the invariant. Correction below.)")

print()
print("=" * 78)
print("Q2 -- reproduce whitepaper 6.4 point 3, and name its convention")
print("=" * 78)

def C(x, gamma=2.0):
    return np.tanh(gamma * np.log1p(x))

def x_from_C(c, gamma=2.0):
    return np.expm1(np.arctanh(c) / gamma)

print(f"  6.4 states C = 0.0003 at ell=100 pc and C = 0.86 at ell=8 kpc (gamma=2).")
print(f"    C=0.0003 -> x = {x_from_C(0.0003):.6g}")
print(f"    C=0.86   -> x = {x_from_C(0.86):.6g}")
print(f"    x ratio  = {x_from_C(0.86)/x_from_C(0.0003):.5g}")
print(f"    (8.0/0.1)^2 = {(8.0/0.1)**2:.5g}   -> the stated x prop ell^2 reproduces.")
print("  Convention assumed: ell enters rho_crit ONLY (via A prop 1/ell^2), rho held fixed.")
print("  This is one-sided: it coarse-grains the threshold and not the field it compares to.")

print()
print("=" * 78)
print("Q3 -- the two self-consistent conventions")
print("=" * 78)

# --- (a) ell smooths BOTH: top-hat sphere of radius ell enclosing M(ell) ---
# rho_ell     = 3 M / (4 pi ell^3)
# rho_crit(l) = 4 pi V^2 / (beta^2 G ell^2)
# x = (3/(16 pi^2)) beta^2 (G M / (ell V^2)) = (3/(16 pi^2)) beta^2 [Vc(ell)/V]^2
coef = 3.0 / (16 * np.pi**2)
print(f"  (a) BOTH smoothed:  x(ell) = {coef:.5f} * beta^2 * [Vc(ell)/V]^2   -- ell cancels")
for b in (1.0, 1.1, 4.5):
    print(f"      ceiling at Vc<=V, beta={b:4.2f}:  x <= {coef*b**2:.5f}   C(gamma=2) = {C(coef*b**2):.4f}")
print("      (independently derived here; matches sibling explorer 2026-08-05)")

# --- (b) ell smooths rho ONLY; rho_crit fixed by GALAXY size (S53's R_half) ---
# NGC 3198, SPARC Lelli+2016 Table 1: Rdisk = 3.14 kpc, Vflat = 150.1 km/s
Rd, Vflat, h = 3.14, 150.1, 0.300      # kpc, km/s, kpc
M_disk = 47.0 * Vflat**4               # M_sun, the archive's BTFR normalization
Sigma0 = M_disk / (2 * np.pi * (Rd*1e3)**2)   # M_sun/pc^2
rho0 = Sigma0 / (2 * h*1e3)                   # M_sun/pc^3, midplane

R_half = 1.678 * Rd                    # kpc, exponential disk
beta_J = 1.1                           # S53's measured calibration
rho_crit_gal = Vflat**2 / (G * (beta_J * R_half*1e3)**2)   # M_sun/pc^3

print(f"\n  (b) rho smoothed at ell, rho_crit fixed by galaxy size (S53's own law):")
print(f"      NGC 3198 (SPARC Table 1): Rdisk={Rd} kpc, Vflat={Vflat} km/s")
print(f"      Sigma0 = {Sigma0:.2f} M_sun/pc^2   rho(0) = {rho0:.4f} M_sun/pc^3")
print(f"      R_half = {R_half:.2f} kpc   rho_crit = {rho_crit_gal:.5g} M_sun/pc^3")
print(f"      {'ell (kpc)':>10} {'rho_ell':>12} {'x':>12} {'C(g=2)':>10}")
for ell in (0.1, 0.3, 1.0, 3.0, 8.0, 20.0):
    # top-hat sphere radius ell at disk centre through an exponential disk of scale h
    ell_pc, h_pc, Rd_pc = ell*1e3, h*1e3, Rd*1e3
    r = np.linspace(1e-3, ell_pc, 4000)
    # cylindrical shell mass inside sphere: sech^2 vertical profile approximated by
    # exponential-in-|z| with scale h; integrate rho(R,z) over the sphere
    z = np.linspace(-ell_pc, ell_pc, 2001)
    RR, ZZ = np.meshgrid(r, z, indexing='ij')
    inside = (RR**2 + ZZ**2) <= ell_pc**2
    rho_grid = (Sigma0 / (2*h_pc)) * np.exp(-RR/Rd_pc) * np.exp(-np.abs(ZZ)/h_pc)
    dV = (2*np.pi*RR) * (r[1]-r[0]) * (z[1]-z[0])
    M_in = np.sum(rho_grid * dV * inside)
    rho_ell = M_in / (4/3*np.pi*ell_pc**3)
    xv = rho_ell / rho_crit_gal
    print(f"      {ell:>10.2f} {rho_ell:>12.5g} {xv:>12.5g} {C(xv):>10.5f}")
print("      => x DECREASES with ell (dilution), the opposite sign from 6.4 point 3,")
print("         and never approaches the knee. Neither consistent convention gives x prop ell^2.")

print()
print("=" * 78)
print("Q4 -- origin of the unexplained 1.515x in rho(0) for NGC 3198")
print("=" * 78)
for Rdv in (2.6, 3.14, 3.2):
    S0 = M_disk / (2*np.pi*(Rdv*1e3)**2)
    print(f"      R_d = {Rdv:4.2f} kpc -> Sigma0 = {S0:7.2f}  rho(0) = {S0/(2*h*1e3):.4f} M_sun/pc^3")
print(f"      (3.2/2.6)^2  = {(3.2/2.6)**2:.4f}")
print(f"      (3.14/2.6)^2 = {(3.14/2.6)**2:.4f}")
print("      => the 1.515x is a DISK SCALE LENGTH disagreement (2.6 vs 3.2 kpc), not a")
print("         density-formula error and not h=198 pc. SPARC Table 1 gives R_d = 3.14 kpc,")
print("         so 3.2 is right to 2% and 2.6 is 21% low. The 0.934 is the wrong number.")

print()
print("=" * 78)
print("Q5 -- the three conventions side by side at a common ell = 100 pc (NGC 3198)")
print("=" * 78)
ell = 0.1
ell_pc, h_pc, Rd_pc = ell*1e3, h*1e3, Rd*1e3
r = np.linspace(1e-3, ell_pc, 3000); z = np.linspace(-ell_pc, ell_pc, 1501)
RR, ZZ = np.meshgrid(r, z, indexing='ij')
rho_g = (Sigma0/(2*h_pc)) * np.exp(-RR/Rd_pc) * np.exp(-np.abs(ZZ)/h_pc)
dV = (2*np.pi*RR) * (r[1]-r[0]) * (z[1]-z[0])
M_in = np.sum(rho_g * dV * ((RR**2 + ZZ**2) <= ell_pc**2))
Vc = np.sqrt(G * M_in / ell_pc)
x_a = coef * beta_J**2 * (Vc/Vflat)**2                      # both smoothed
x_b = (M_in / (4/3*np.pi*ell_pc**3)) / rho_crit_gal         # rho only
x_c = x_from_C(0.0003)                                      # rho_crit only (whitepaper 08-06)
print(f"      M(<100 pc) = {M_in:.4g} M_sun, Vc = {Vc:.3f} km/s, Vc/Vflat = {Vc/Vflat:.4f}")
print(f"      (a) ell in BOTH        x = {x_a:.4g}")
print(f"      (b) ell in rho ONLY    x = {x_b:.4g}")
print(f"      (c) ell in rho_crit ONLY x = {x_c:.4g}   <- the 2026-08-06 whitepaper claim")
print(f"      spread (b)/(c) = {x_b/x_c:.4g}   ((a) and (c) coincide to ~30% by numerical accident)")
