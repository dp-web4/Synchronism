#!/usr/bin/env python3
"""
Publisher 2026-08-19 — independent verification of the DE-sector locality-fork identities
(site explorer finding 2026-08-18, routed via Research/proposals/de_locality_fork_
perturbation_channel_factor_two_20260818.md).

Model (Sessions #100/#101, as the site states it):
    C(x)  = tanh(gamma * ln(1+x)),   x = rho_m / rho_crit
    rho_DE = rho_m (1 - C)/C

Checks (all symbolic, exact, for all x > 0 and gamma > 0):
  1. C closed form                 = ((1+x)^{2g} - 1)/((1+x)^{2g} + 1)
  2. rho_DE/rho_crit closed form   = 2x/((1+x)^{2g} - 1)
  3. gamma = 1/2  =>  rho_DE/rho_crit == 2 EXACTLY, and d(rho_DE)/dx == 0
     => rho_DE is constant in DENSITY, hence in SPACE: LCDM non-perturbatively.
  4. delta_DE/delta_m = dln rho_DE / dln x = 1 + w_DE, with leading term
     eps*(ln(1+x)/x - 1), eps = 2*gamma - 1  -- linear in eps, NO order-eps^0 term,
     and identically ZERO at gamma = 1/2.

Also recorded (arithmetic, not symbolic): the "locality improves the required precision
by 7%" claim in the source finding is a CENTRAL-VALUE artifact, not a gain from locality.
sigma_gamma = (1/2 - gamma_hat)/3 is horn-independent by construction.
"""
import sympy as sp

x, g, eps = sp.symbols('x gamma epsilon', positive=True)
u = (1 + x) ** (2 * g)

C = sp.tanh(g * sp.log(1 + x)).rewrite(sp.exp).simplify()
assert sp.simplify(sp.powsimp(C - (u - 1) / (u + 1), force=True)) == 0
print("1. C closed form                      OK")

rho = sp.simplify(sp.powsimp((x * (1 - C) / C).rewrite(sp.exp), force=True))
assert sp.simplify(sp.powsimp(rho - 2 * x / (u - 1), force=True)) == 0
print("2. rho_DE/rho_crit closed form        OK   ->", rho)

r_half = sp.simplify(sp.powsimp(rho.subs(g, sp.Rational(1, 2)), force=True))
assert r_half == 2 and sp.simplify(sp.diff(r_half, x)) == 0
print("3. rho_DE(gamma=1/2) = 2, drho/dx = 0 OK   (constant in density AND space)")

q = sp.simplify(sp.powsimp(x * sp.diff(rho, x) / rho, force=True))   # = 1 + w_DE
assert sp.simplify(q.subs(g, sp.Rational(1, 2))) == 0
lin = sp.simplify(sp.series(sp.simplify(sp.powsimp(q.subs(g, (1 + eps) / 2), force=True)),
                            eps, 0, 2).removeO())
assert sp.simplify(lin - eps * (sp.log(1 + x) / x - 1)) == 0
print("4. delta_DE/delta_m = 1+w             OK   -> eps*(ln(1+x)/x - 1) + O(eps^2), 0 at gamma=1/2")

# --- the 7% claim, arithmetic ---
for name, ghat, quoted in (("08-12 background horn", 0.487, 0.004),
                           ("08-18 locality horn",   0.489, 0.0037)):
    print(f"   sigma_gamma[{name}] = (0.5-{ghat})/3 = {(0.5-ghat)/3:.5f}  (source quotes {quoted})")
print("   => identical formula, no horn dependence: locality buys 0% in sigma_gamma.")
print("   The real x2 is in sigma(fsigma8): source's own Horn N/L peaks 0.223% / 0.444%")
print("   => required precision 0.0743% -> 0.148%; DESI DR2 delivers ~1-3%/bin (~10x short).")
print("\nALL CHECKS PASS")
