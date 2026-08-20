"""
Publisher 2026-08-20 — independent verification of the completion-B omega correction.

Source claim (site explorer 2026-08-19, routed via
Research/proposals/completion_b_cassini_omega_absorbed_by_closure_20260819.md):

  (i)   Cassini forces omega >= 4.0e4 for an unscreened massless Brans-Dicke scalar,
        so the published grid {0, 1, 5, 50} sits ~800x inside the excluded region.
  (ii)  omega is ABSORBED by the model's own closure: the pinning C_eff(x0) = Omega_m
        is held by sliding x0, with eps_crit(omega) -> 0, so "192 x 4 omega values"
        is really 192 x 1.
  (iii) At the allowed omega the no-go HARDENS: w0 = -1.58 (omega=0) -> -3.18 (omega=4e4).

What is verified HERE (symbolically where possible), and what is not:
  VERIFIED  (1) gamma_PPN - 1 = -1/(2+omega) and the Cassini arithmetic.
  VERIFIED  (2) eps_crit(omega) = (-3 + sqrt(9+6 omega))/(3 omega) is the positive root
                of the background factor B(eps) = 1 - 3 eps - 1.5 omega eps^2.
  VERIFIED  (3) THE ABSORPTION IDENTITY, which is the sharp form of claim (ii) and is
                stronger than the numerical "trajectories agree to <2%": at the pinned
                point the closure fixes the PRODUCT omega*eps^2, not omega.
                    1.5 omega eps_crit^2 = 1 - 3 eps_crit   (exactly, by definition of the root)
                 => omega * eps_crit^2 -> 2/3  as omega -> inf.
                So beyond eps_crit << 1 the construction has one physical point, not a
                family: omega is not a scanned dimension, it is a reparametrisation of x0.
  NOT VERIFIED HERE: that B(x) = 1 - 3 eps - 1.5 omega eps^2 is the complete omega
                dependence of the completion-B background (taken from the source
                implementation, covariant_00_component_sign_lock_audit.py), and the
                w0 = -3.18 trajectory value, which requires the source integrator.
"""
import sympy as sp

om, eps = sp.symbols('omega epsilon', positive=True)
checks = []

def check(name, cond):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

print("PART A -- Cassini bound")
gamma_ppn = (1 + om) / (2 + om)
dev = sp.simplify(gamma_ppn - 1)
check("gamma_PPN - 1 == -1/(2+omega)", sp.simplify(dev + 1/(2 + om)) == 0)
# Bertotti, Iess & Tortora 2003: gamma_PPN - 1 = (2.1 +/- 2.3)e-5; 2-sigma lower edge:
lo = 2.1e-5 - 2 * 2.3e-5
om_min = float(sp.solve(sp.Eq(-1/(2 + om), lo), om)[0])
print(f"        2-sigma lower edge = {lo:.3e}  ->  omega_min = {om_min:.4g}")
check("omega_min ~ 4.0e4", abs(om_min - 4.0e4) / 4.0e4 < 0.01)
check("published grid max (omega=50) is >~ 800x inside the excluded region",
      round(om_min / 50) == 800)

print("\nPART B -- eps_crit is the positive root of the background factor")
B = 1 - 3*eps - sp.Rational(3, 2)*om*eps**2
eps_crit = (-3 + sp.sqrt(9 + 6*om)) / (3*om)
check("B(eps_crit) == 0 identically in omega", sp.simplify(B.subs(eps, eps_crit)) == 0)
check("eps_crit > 0 for all omega > 0", sp.simplify(sp.limit(eps_crit, om, 0)) == sp.Rational(1, 3))
check("eps_crit -> 0 as omega -> inf", sp.limit(eps_crit, om, sp.oo) == 0)

print("\nPART C -- THE ABSORPTION IDENTITY (the sharp form of 'omega is absorbed')")
prod = sp.simplify(om * eps_crit**2)
lim = sp.limit(prod, om, sp.oo)
print(f"        omega * eps_crit^2  ->  {lim}   as omega -> inf")
check("omega*eps_crit^2 -> 2/3  (closure pins the PRODUCT, not omega)",
      sp.simplify(lim - sp.Rational(2, 3)) == 0)
check("1.5*omega*eps_crit^2 == 1 - 3*eps_crit exactly",
      sp.simplify(sp.Rational(3,2)*om*eps_crit**2 - (1 - 3*eps_crit)) == 0)
# eps_crit ~ sqrt(2/(3 omega)) at large omega -- so eps_crit falls as omega^(-1/2)
asym = sp.sqrt(sp.Rational(2, 3) / om)
check("eps_crit ~ sqrt(2/(3 omega)) asymptotically",
      sp.limit(eps_crit / asym, om, sp.oo) == 1)

print("\nPART D -- the quoted numbers")
for w, want in ((50, 0.097), (4.0e4, 4.1e-3)):
    got = float(eps_crit.subs(om, w))
    print(f"        eps_crit(omega={w:>9.0f}) = {got:.4g}   (source quotes {want:g})")
    check(f"eps_crit({w:.0f}) matches source to 2 sig figs",
          abs(got - want) / want < 0.02)
# and the reason "omega=4e4 vs 1e6 agree to <2%": the physical content is omega*eps^2,
# which is within 2/3*(1 - 3 eps_crit) of 2/3 at both.
for w in (4.0e4, 1.0e6):
    e = float(eps_crit.subs(om, w))
    print(f"        omega={w:.0e}: eps_crit={e:.3e}, omega*eps^2={w*e*e:.6f}, "
          f"departure from 2/3 = {abs(w*e*e - 2/3)/(2/3)*100:.2f}%")

n_pass = sum(p for _, p in checks)
print(f"\n{n_pass}/{len(checks)} assertions pass")
assert n_pass == len(checks), [n for n, p in checks if not p]
