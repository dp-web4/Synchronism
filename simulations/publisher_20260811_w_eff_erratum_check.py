#!/usr/bin/env python3
"""Publisher 2026-08-11: independent verification of the 2026-08-10 dark-energy erratum.

Checks, against Research/proposals/dark_energy_sector_exists_and_forbids_desi_quadrant_20260810.md
and Research/Session100_Modified_Friedmann.md:

  1. Continuity-equation sign: w = -1 - (1/3) dln(rho_DE)/dln(a).
     Session #100/#101 published w_eff = -1 + (1/3) dln(rho_DE)/dln(a).
     Unit test both on matter (w=0) and radiation (w=1/3).
  2. The model: C(x) = tanh(gamma*ln(1+x)), x = x0*(1+z)^3, rho_DE = rho_m*(1-C)/C.
     Verify x0 from C(x0)=Omega_m at gamma=2, the corrected w(z) table
     (w(0) ~ -1.24, w(0.1) ~ -1.3186, ...), and that Session #100's published
     column equals the stated formula with the leading -1 dropped.
  3. gamma=1/2 identity: tanh((1/2)ln(1+x)) == x/(x+2), and with C0=Omega_m the
     coherence C(z) == Omega_m(z) (flat LCDM) exactly, hence w(z) == -1.
  4. Limits: w -> -2*gamma as z -> inf; w -> -1 as a -> inf (any gamma, x0).
  5. Quadrant no-go: sign(w0+1) == sign(wa) over gamma in [0.05, 20]
     (CPL wa = -dw/da at a=1); the DESI-preferred quadrant (w0 > -1 AND wa < 0)
     is unreachable.

All symbolic work in sympy (per 2026-08-06 rule: no np.gradient on log grids).
"""

import sympy as sp

a, x, x0, g, z = sp.symbols('a x x0 gamma z', positive=True)
Om = sp.Symbol('Omega_m', positive=True)

print("=" * 70)
print("CHECK 1: continuity-equation sign")
print("=" * 70)
# rho ~ a^(-3(1+w))  =>  dln rho/dln a = -3(1+w)  =>  w = -1 - (1/3) dln rho/dln a
w_sym = sp.Symbol('w')
rho = a ** (-3 * (1 + w_sym))
dlnrho_dlna = sp.simplify(a * sp.diff(sp.log(rho), a))
w_correct = sp.simplify(-1 - sp.Rational(1, 3) * dlnrho_dlna)
w_published = sp.simplify(-1 + sp.Rational(1, 3) * dlnrho_dlna)
print(f"  correct formula recovers w:   {w_correct}  (must be w)")
assert w_correct == w_sym
for name, wval in [("matter", 0), ("radiation", sp.Rational(1, 3))]:
    pub = w_published.subs(w_sym, wval)
    print(f"  published formula on {name:9s} (true w={wval}): returns {pub}")
assert w_published.subs(w_sym, 0) == -2          # proposal: -2 for matter
assert w_published.subs(w_sym, sp.Rational(1, 3)) == sp.Rational(-7, 3)  # -7/3 radiation
print("  PASS: published form is wrong-signed exactly as the erratum states")

print()
print("=" * 70)
print("CHECK 2: corrected w(z) at gamma=2, and the published table's diagnosis")
print("=" * 70)
C = sp.tanh(g * sp.log(1 + x))
xa = x0 * a ** (-3)                      # x(a), with a = 1/(1+z)
Ca = C.subs(x, xa)
rho_DE = a ** (-3) * (1 - Ca) / Ca       # rho_m0 factors out of dln
w_of_a = sp.simplify(-1 - sp.Rational(1, 3) * a * sp.diff(sp.log(rho_DE), a))
# closed form from the proposal: w = -gamma*(1+C)*x / (C*(1+x))
w_closed = (-g * (1 + C) * x / (C * (1 + x))).subs(x, xa)
diff_closed = sp.simplify(w_of_a - w_closed)
print(f"  w(a) - closed_form simplifies to: {diff_closed}")
assert diff_closed == 0
# x0 at gamma=2 from C(x0) = Omega_m = 0.30 (Session #100's own input, line 74)
x0_num = sp.nsolve(sp.tanh(2 * sp.log(1 + x)) - sp.Rational(30, 100), x, 0.2)
print(f"  x0 (gamma=2, C0=0.30): {float(x0_num):.5f}   (proposal: 0.16738)")
assert abs(float(x0_num) - 0.16738) < 5e-5
wfun = sp.lambdify(a, w_of_a.subs({g: 2, x0: x0_num}))
table = {0.0: -1.24, 0.1: -1.3186, 0.5: -1.7329, 1.0: -2.3690, 2.0: -3.2788}
for zz, expected in table.items():
    got = float(wfun(1 / (1 + zz)))
    print(f"  z={zz:3.1f}: corrected w = {got:+.4f}   (proposal: {expected:+.4f})")
    assert abs(got - expected) < 2e-2 if zz == 0 else abs(got - expected) < 2e-4
# Session #100's published column = stated formula with leading -1 dropped,
# i.e. +(1/3) dln rho_DE/dln a  ==  w_correct + 1 ... no: dropped -1 from the
# *stated* (wrong-signed) formula: +(1/3)dlnrho/dlna = -(w+1) under the correct sign.
w_nolead = sp.simplify(w_of_a + 1)                    # = -(1/3) dln rho_DE/dln a
pub_col = sp.lambdify(a, (-w_nolead).subs({g: 2, x0: x0_num}))
published = {0.1: 0.32, 0.5: 0.73, 1.0: 1.37, 2.0: 2.28}
for zz, pub in published.items():
    got = float(pub_col(1 / (1 + zz)))
    print(f"  z={zz:3.1f}: '+(1/3)dln rho/dln a' = {got:+.3f}   (Session #100 table: {pub:+.2f})")
    assert abs(got - pub) < 6e-3
print("  PASS: published column is the stated formula with the leading -1 dropped")

print()
print("=" * 70)
print("CHECK 3: gamma=1/2 branch is identically the LCDM background")
print("=" * 70)
mobius = sp.simplify((sp.tanh(sp.Rational(1, 2) * sp.log(1 + x)) - x / (x + 2)).rewrite(sp.exp))
print(f"  tanh(0.5*ln(1+x)) - x/(x+2) = {mobius}")
assert mobius == 0
# C0 = Omega_m  =>  x0 = 2*Om/(1-Om); then C(z) should equal Omega_m(z)
x0_half = 2 * Om / (1 - Om)
C_half = (x / (x + 2)).subs(x, x0_half * (1 + z) ** 3)
Om_z = Om * (1 + z) ** 3 / (Om * (1 + z) ** 3 + (1 - Om))
print(f"  C(z; gamma=1/2, C0=Om) - Omega_m(z) = {sp.simplify(C_half - Om_z)}")
assert sp.simplify(C_half - Om_z) == 0
w_half = sp.simplify(w_of_a.subs({g: sp.Rational(1, 2), x0: x0_half}).rewrite(sp.exp))
print(f"  w(a; gamma=1/2) = {w_half}")
assert sp.simplify((w_half + 1).rewrite(sp.exp)) == 0
print("  PASS: gamma=1/2 with C0=Om is Omega_m(z) exactly; w == -1 (it IS LCDM)")

print()
print("=" * 70)
print("CHECK 4: limits")
print("=" * 70)
lim_early = sp.limit(w_of_a, a, 0)       # z -> inf
lim_late = sp.limit(w_of_a, a, sp.oo)    # a -> inf
print(f"  w(a->0)  = {lim_early}   (proposal: -2*gamma)")
print(f"  w(a->oo) = {lim_late}   (proposal: -1)")
assert sp.simplify(lim_early + 2 * g) == 0
assert lim_late == -1
print("  PASS")

print()
print("=" * 70)
print("CHECK 5: quadrant no-go — sign(w0+1) == sign(wa), no phantom crossing")
print("=" * 70)
wa_expr = -sp.diff(w_of_a, a)            # CPL: w(a) ~ w0 + wa*(1-a); wa = -dw/da|a=1
bad = 0
gammas = [0.05, 0.1, 0.2, 0.3, 0.489, 0.5, 0.6, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0]
for gv in gammas:
    # C0 = Omega_m = 0.315 (Planck); the sign lock is calibration-independent.
    # Closed form: tanh(g*ln(1+x0)) = Om  =>  x0 = exp(atanh(Om)/g) - 1
    x0v = sp.exp(sp.atanh(sp.Rational(315, 1000)) / gv) - 1
    w0v = float(w_of_a.subs({g: gv, x0: x0v, a: 1}))
    wav = float(wa_expr.subs({g: gv, x0: x0v, a: 1}))
    same_sign = (w0v + 1) * wav >= -1e-12
    in_desi_quadrant = (w0v > -1) and (wav < -1e-12)
    flag = "" if same_sign and not in_desi_quadrant else "  <-- VIOLATION"
    print(f"  gamma={gv:6.3f}: w0={w0v:+.4f}  wa={wav:+.4f}  sign-locked={same_sign}{flag}")
    if not same_sign or in_desi_quadrant:
        bad += 1
print(f"  DESI-preferred quadrant (w0 > -1 AND wa < 0) reached: {bad}/{len(gammas)}")
assert bad == 0
print("  PASS: 0 of 16 gammas reach the DESI quadrant; w=-1 crossing impossible")

print()
print("ALL CHECKS PASS — the 2026-08-10 erratum and no-go verify independently.")
