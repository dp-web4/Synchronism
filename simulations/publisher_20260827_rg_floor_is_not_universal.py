#!/usr/bin/env python3
"""
Publisher 2026-08-27 — the RG floor eps_0 read from the SOURCE, not the abstract.

Closes the coverage caveat stated in the 2026-08-26 pass ("abstracts only; the
equation forms came through the routing note and were re-derived, not read in the
papers").  arXiv:1603.04943 was fetched and text-extracted locally; Eq. (2.6),
Eq. (4.1) and every fitted parameter below are read from the paper's own text.

Two results:

  (1) IDENTITY, against the VERBATIM published Eq. (4.1) rather than a re-derived
      form.  C_Omega === M&D's gravitational permittivity, exactly, in both log
      bases.  (The 08-26 check reached the same conclusion from a relayed form.)

  (2) THE FLOOR IS NOT A UNIVERSAL CONSTANT IN THE PRIOR ART.  The 08-26 pass
      recorded "RG fits eps_0 = 0.20-0.25" and opened a discriminator said to gate
      on execution of Cesare+2020.  That range is the two SPIRALS only.  The same
      paper fits eps_0 = 0.045 and 0.065 for two CLUSTERS and 1/6 for its model
      galaxy -- a factor 5.6 spread, every value BELOW the derived 0.315.

Run: python3 publisher_20260827_rg_floor_is_not_universal.py
"""
import math
import sympy as sp

FAIL = []


def check(label, cond, detail=""):
    print(f"  [{'OK ' if cond else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    if not cond:
        FAIL.append(label)


# ---------------------------------------------------------------- (1) identity
print("(1) C_Omega  vs  M&D 2016 Eq. (4.1), verbatim from arXiv:1603.04943 p.16")
print("    eps(rho) = eps_0 + (1-eps_0) * (1/2) * { tanh[ log (rho/rho_c)^q ] + 1 }")

u, q, Om = sp.symbols('u q Om', positive=True)
for name, ln_base in (('natural log', sp.Integer(1)), ('log base 10', sp.log(10))):
    z = q * sp.log(u) / ln_base                                  # log_base(u**q)
    eps = Om + (1 - Om) * sp.Rational(1, 2) * (sp.tanh(z) + 1)   # published form
    p = 2 * q / ln_base
    C = Om + (1 - Om) * (u ** p) / (1 + u ** p)                  # repo C_Omega
    d = sp.simplify(sp.expand(sp.powsimp((eps - C).rewrite(sp.exp), force=True)))
    check(f"{name:12s} p = 2q/ln(base) -> eps - C_Omega == 0", d == 0, f"got {d}")

# independent numeric confirmation over 8 decades
worst = 0.0
for base10 in (False, True):
    lb = math.log(10) if base10 else 1.0
    L = math.log10 if base10 else math.log
    for k in range(-40, 41):
        r = 10 ** (k / 10)                       # rho/rho_c
        a = 0.315 + (1 - 0.315) * 0.5 * (math.tanh(L(r ** 0.75)) + 1)
        X = r ** (2 * 0.75 / lb)
        b = 0.315 + (1 - 0.315) * X / (1 + X)
        worst = max(worst, abs(a - b))
check("numeric max|eps - C_Omega| over 8 decades, both bases", worst < 1e-12,
      f"{worst:.3e}")

# ------------------------------------------------------------------ (2) floors
OM = 0.315                       # framework: DERIVED, universal, = Omega_m
FITS = [                         # M&D 2016, read from the paper's own text
    ("NGC 6946   spiral   q=3/4", 0.25,  "sec 4.2, p.20"),
    ("NGC 1560   spiral   q=3/4", 0.20,  "sec 4.2, p.20"),
    ("model G0   disc     q=1/2", 1 / 6, "sec 4.1, fig.12 caption"),
    ("A1991      cluster  q=2",   0.045, "sec 4.3, p.24"),
    ("A1795      cluster  q=2",   0.065, "sec 4.3, p.24"),
]

print("\n(2) the floor: one DERIVED universal value vs five FITTED per-system values")
print(f"    framework derives a single universal eps_0 = Omega_m = {OM}")
for label, v, where in FITS:
    print(f"      {label:28s} eps_0 = {v:.4f}   0.315/eps_0 = {OM / v:5.2f}x"
          f"   ceiling 1/eps_0 = {1 / v:5.2f}   [{where}]")

lo = min(v for _, v, _ in FITS)
hi = max(v for _, v, _ in FITS)
print(f"\n    RG's fitted eps_0 spans {lo} -> {hi} = factor {hi / lo:.2f}")
check("every fitted eps_0 lies BELOW the derived floor", all(v < OM for _, v, _ in FITS))
check("the '0.20-0.25' range covers only the two spirals",
      sorted(v for _, v, _ in FITS)[-2:] == [0.20, 0.25])
check("max discrepancy is ~7x, not ~30%", abs(OM / lo - 7.0) < 0.05, f"{OM / lo:.2f}x")

# the ceiling inequality, which is TEST-09/10's inequality on a new dataset
B_FRAMEWORK = 1 / OM
print(f"\n    framework's boost ceiling  1/Omega_m = {B_FRAMEWORK:.2f}")
for label, v, _ in FITS:
    need = 1 / v
    if need > B_FRAMEWORK:
        print(f"      {label:28s} RG needs {need:5.2f}x  > {B_FRAMEWORK:.2f}  -> ceiling binds")
check("RG's cluster fits demand boosts the framework's ceiling forbids",
      1 / 0.045 > B_FRAMEWORK and 1 / 0.065 > B_FRAMEWORK)

# ------------------------------------------------- exponent: NOT base-invariant
phi = (1 + 5 ** 0.5) / 2
print("\n(3) the exponent comparison is NOT safe -- it flips with the log base")
print(f"    1/phi = 2q/ln(base)  =>  natural log: q = {1 / (2 * phi):.4f}"
      f"   |  log10: q = {math.log(10) / (2 * phi):.4f}")
print( "    M&D fit q = 0.75 (spirals), 2 (clusters); Cesare+2020 fit Q ~ 0.47.")
print( "    M&D write 'ln' elsewhere (sec 3.4) and 'log' in Eq. (4.1); Cesare+2020's")
print( "    reported form uses 'ln'.  The base is not stated in either relay, and it")
print( "    moves q by ln(10) = 2.30x.  Only the FLOOR is base-independent, so only")
print( "    the floor is quoted as a comparison.")

print("\n" + ("ALL CHECKS PASSED" if not FAIL else f"FAILURES: {FAIL}"))
raise SystemExit(1 if FAIL else 0)
