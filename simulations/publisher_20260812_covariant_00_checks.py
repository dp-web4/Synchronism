#!/usr/bin/env python3
"""Publisher 2026-08-12: independent verification of the covariant 00-component finding
(synchronism-site/explorer/findings/covariant-00-component-sign-lock-dies-desi-nogo-hardens.md,
2026-08-11) BEFORE propagating it into whitepaper §5.15 / Appendix D / Session #100 erratum.

Checks (each assert-fatal):
  1. Class identity: for any dark energy slaved to matter, rho_DE = rho_m * F(x) with
     x = x0 * a^-3, continuity gives w_DE = dlnF/dlnx exactly (sympy, symbolic).
  2. Completion A (Appendix D as written): Bianchi forces rho/C ~ a^-3, so H^2 ~ a^-3 (EdS);
     vacuum floor rho/C -> rho_crit/gamma as rho->0 (sympy limit), rho/C monotone => floor;
     a_end = (gamma*x0/Omega_m)^(1/3) = 1.0372 under Session #100's calibration
     (gamma=2, C0=Omega_m=0.3 => x0=0.16738).
  3. Completion B closure: substituting Cdot/C = -3*eps*H into the Brans-Dicke 00-component
     H^2 = 8piG rho_m/(3C) - H(Cdot/C) + (omega/6)(Cdot/C)^2 closes to
     H^2 = 8piG rho_m/(3 C B) with B = 1 - 3*eps - (3*omega/2)*eps^2 (sympy, symbolic);
     and eps = dlnC/dlnx = gamma*x*(1-C^2)/(C*(1+x)) (sympy, symbolic).
  4. Completion B numbers (omega=0): calibration C_eff(x0) = Omega_m = 0.3 forced; reproduce
     the finding's table -- x0, C0, w(0), w(1), a_rip -- for gamma in {0.2, 0.489, 0.5, 2.0};
     B -> -2 as x -> 0 (attractor destroyed); DESI-quadrant spot check: w0 > -1 members have
     wa > 0 (slope of w against (1-a) positive near a=1), i.e. the quadrant
     (w0 > -1 AND wa < 0) is not populated at the table points.

NOT re-run here: the 192-gamma x 4-omega dense scan, the CPL fits, and the BAO-shape rms
(those remain the site script's; this file reproduces the algebra and the table anchors).
"""

import sympy as sp
import numpy as np
from scipy.optimize import brentq

print("=" * 72)
print("CHECK 1 -- class identity: w_DE = dlnF/dlnx for rho_DE = rho_m*F(x), x ~ a^-3")
a, x0, gam, om = sp.symbols("a x0 gamma omega", positive=True)
F = sp.Function("F", positive=True)
x = x0 * a ** -3
rho_m = sp.Symbol("rho_m0", positive=True) * a ** -3
rho_DE = rho_m * F(x)
w = -1 - sp.Rational(1, 3) * sp.diff(sp.log(rho_DE), a) * a
# dlnF/dlnx evaluated at x = x0*a^-3
dlnF_dlnx = (sp.diff(F(x), a) / F(x)) / (sp.diff(sp.log(x), a))
assert sp.simplify(w - dlnF_dlnx) == 0
print("  PASS: w_DE = dlnF/dlnx, exactly, any F (symbolic)")

print("=" * 72)
print("CHECK 2 -- completion A: EdS + vacuum floor + a_end")
xs = sp.Symbol("x", positive=True)
C = sp.tanh(gam * sp.log(1 + xs))
# vacuum floor: rho/C = rho_crit * x / C(x) -> rho_crit/gamma as x -> 0
floor = sp.limit(xs / C, xs, 0)
assert sp.simplify(floor - 1 / gam) == 0
print("  PASS: rho/C -> rho_crit/gamma as rho -> 0 (symbolic limit)")
# monotone: d(x/C)/dx > 0 over wide range (floor is the infimum) -- symbolic derivative,
# evaluated numerically (a raw np.diff on the flat small-x tail is cancellation noise)
g2 = 2.0
dfdx = sp.lambdify(xs, sp.diff(xs / C, xs).subs(gam, g2), "numpy")
xv = np.logspace(-6, 4, 2000)
assert np.all(dfdx(xv) > 0)
print("  PASS: d(x/C)/dx > 0 (gamma=2 grid, 1e-6..1e4) => rho_crit/gamma is a floor")
# EdS: if rho/C = K a^-3 then H^2 = (8piG/3) K a^-3 -- trivially a^-3; deceleration q = +1/2
t = sp.Symbol("t", positive=True)
af = sp.Function("a", positive=True)
K8 = sp.Symbol("K", positive=True)
# H^2 = K a^-3  =>  adot = sqrt(K) a^(-1/2); q = -a*addot/adot^2
adot = sp.sqrt(K8) * af(t) ** sp.Rational(-1, 2)
addot = sp.diff(adot, t).subs(sp.Derivative(af(t), t), adot)
q = sp.simplify(-af(t) * addot / adot ** 2)
assert sp.simplify(q - sp.Rational(1, 2)) == 0
print("  PASS: H^2 ~ a^-3 gives q = +1/2 (Einstein-de Sitter, no acceleration)")
# a_end under Session #100 calibration: C(x0)=0.3 at gamma=2; a_end^3 = gamma*x0/Omega_m
Om_m = 0.3
Cnum = sp.lambdify(xs, C.subs(gam, g2), "numpy")
x0_100 = brentq(lambda u: Cnum(u) - Om_m, 1e-6, 10)
assert abs(x0_100 - 0.16738) < 5e-4, x0_100
a_end = (g2 * x0_100 / Om_m) ** (1.0 / 3.0)
assert abs(a_end - 1.0372) < 5e-4, a_end
print(f"  PASS: x0 = {x0_100:.5f} (S100's 0.16738); a_end = {a_end:.4f} (finding: 1.0372)")

print("=" * 72)
print("CHECK 3 -- completion B closure: B = 1 - 3*eps - (3*omega/2)*eps^2; eps formula")
H, rm, Csym, eps = sp.symbols("H rho_m C epsilon", positive=True)
omega_bd = sp.Symbol("omega_BD")
Cdot_over_C = -3 * eps * H
lhs = H ** 2
rhs = rm / (3 * Csym) - H * Cdot_over_C + (omega_bd / 6) * Cdot_over_C ** 2  # 8piG := 1
# solve lhs = rhs for H^2:  H^2 (1 - 3 eps - (3 omega/2) eps^2) = rho_m/(3C)
B_expr = sp.simplify((rhs - (rm / (3 * Csym))) / H ** 2)
B_closed = sp.simplify(1 - B_expr)
assert sp.simplify(B_closed - (1 - 3 * eps - sp.Rational(3, 2) * omega_bd * eps ** 2)) == 0
print("  PASS: 00-component closes to H^2 = 8piG rho_m/(3*C*B), B = 1 - 3eps - (3w/2)eps^2")
eps_derived = sp.simplify(sp.diff(sp.log(C), xs) * xs)
eps_claimed = gam * xs * (1 - C ** 2) / (C * (1 + xs))
assert sp.simplify(eps_derived - eps_claimed) == 0
print("  PASS: eps = dlnC/dlnx = gamma*x*(1-C^2)/(C*(1+x)) (symbolic)")
# B -> -2 as x -> 0 at omega=0 (eps -> 1): attractor destroyed
eps0 = sp.limit(eps_derived, xs, 0)
assert sp.simplify(eps0 - 1) == 0
print("  PASS: eps -> 1 as x -> 0, so B -> -2 (omega=0): w = -1 attractor destroyed")

print("=" * 72)
print("CHECK 4 -- completion B table anchors (omega=0), finding section 3")
def Cf(u, g):
    return np.tanh(g * np.log1p(u))

def epsf(u, g):
    c = Cf(u, g)
    return g * u * (1 - c * c) / (c * (1 + u))

def Ceff(u, g):
    return Cf(u, g) * (1 - 3 * epsf(u, g))

def w_of_x(u, g):
    # F(x) = (1 - Ceff)/Ceff ; w = dlnF/dlnx (check 1), numerically via log-derivative
    h = 1e-6
    lF = lambda v: np.log((1 - Ceff(v, g)) / Ceff(v, g))
    return (lF(u * (1 + h)) - lF(u * (1 - h))) / (2 * h)

table = {  # gamma: (x0, C0, a_rip, w0, w1)
    0.2:   (44.7, 0.64, 1.68, -0.92, -0.54),
    0.489: (7.94, 0.79, 1.25, -2.12, -1.06),
    0.5:   (7.65, 0.79, 1.24, -2.16, -1.08),
    2.0:   (1.08, 0.90, 1.10, -4.91, -3.50),
}
desi_quadrant_hit = False
for g, (x0_t, C0_t, arip_t, w0_t, w1_t) in table.items():
    # calibration: C_eff(x0) = Omega_m = 0.3 (C0*B0 = Om forced)
    x0n = brentq(lambda u: Ceff(u, g) - Om_m, 1e-4, 1e4)
    C0n = Cf(x0n, g)
    # rip: Ceff crosses zero going to low x (future); x_rip < x0
    x_rip = brentq(lambda u: Ceff(u, g), 1e-8, x0n)
    a_rip = (x0n / x_rip) ** (1.0 / 3.0)
    w0n = w_of_x(x0n, g)
    w1n = w_of_x(x0n * 8.0, g)  # z=1: x = x0*(1+z)^3
    ok = (abs(x0n - x0_t) / x0_t < 0.01 and abs(C0n - C0_t) < 0.01
          and abs(a_rip - arip_t) < 0.01 and abs(w0n - w0_t) < 0.02
          and abs(w1n - w1_t) < 0.02)
    # DESI quadrant needs w0 > -1 AND wa < 0. CPL: w(a) = w0 + wa*(1-a), so for a < 1
    # (the past) wa = (w_past - w0)/(1-a):  wa > 0 iff w was HIGHER in the past.
    # Spot check at a = 0.9 (z ~ 0.111).
    w_past = w_of_x(x0n * (1 / 0.9) ** 3, g)
    wa_sign_positive = (w_past - w0n) > 0
    if (w0n > -1) and not wa_sign_positive:
        desi_quadrant_hit = True
    print(f"  gamma={g:<6} x0={x0n:8.2f} C0={C0n:.2f} a_rip={a_rip:.2f} "
          f"w(0)={w0n:+.2f} w(1)={w1n:+.2f}  table-match={'PASS' if ok else 'FAIL'} "
          f"wa_sign={'+' if wa_sign_positive else '-'}")
    assert ok, f"table mismatch at gamma={g}"
assert not desi_quadrant_hit
print("  PASS: all four table rows reproduced; no table member reaches the DESI quadrant")
print("        (every w0 > -1 member has wa > 0; every wa-negative slot needs w0 < -1)")

print("=" * 72)
print("ALL CHECKS PASSED -- algebra + table anchors reproduced independently.")
print("Not re-run: 192-gamma dense scan, CPL fits, BAO rms (site script's numbers).")
