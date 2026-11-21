"""Publisher 2026-11-13: is "eta == Abrikosov-Gor'kov" decided by "AG lowers T_c"?

The paper (superconductivity.md, S616 audit) prints T_c = Delta / (1.76 k_B eta), 607 K for YBCO, and
"eta == AG pair-breaking efficiency ... identical to the standard AG framework".
PREDICTIONS.md 2026-11-12 (iii) calls the identity a label: "AG pair breaking enters through the digamma
equation and lowers T_c". That compares against T_c0. The printed formula compares T_c against the GAP:
eta = Delta / (1.764 k T_c). So the question is what AG does to Delta/(1.764 k T_c), and that depends on
which Delta is meant:
  (a) the AG order parameter Delta(0), or
  (b) the spectral gap Omega_g = Delta (1 - (alpha/Delta)^(2/3))^(3/2), what tunnelling/ARPES measures.
Only the gapped branch (alpha <= Delta, alpha/alpha_c <= 0.912) is used; units Tc0 = k = 1.

Run: python3 -I eta_vs_abrikosov_gorkov_20261113.py
"""
import numpy as np
from scipy.special import digamma
from scipy.optimize import brentq

kB = 8.617333262e-5  # eV/K
print("Printed formula, site inputs (Delta = 35 meV):")
for eta in (0.38, 1.0):
    print(f"  eta = {eta:4.2f}: T_c = {0.035 / (1.76 * kB * eta):6.1f} K")

D00 = np.pi / np.exp(np.euler_gamma)          # 1.764, BCS Delta(0)/k Tc0
a_c = D00 / 2                                  # AG critical pair-breaking energy at T = 0
a_gapless = D00 * np.exp(-np.pi / 4)           # alpha = Delta at the gapless onset


def Tc(a):
    if a == 0:
        return 1.0
    return brentq(lambda t: np.log(1 / t) - (digamma(0.5 + a / (2 * np.pi * t)) - digamma(0.5)), 1e-12, 1.0)


def Delta0(a):  # gapped branch: ln(D00/D) = pi a / (4 D)
    if a == 0:
        return D00
    return brentq(lambda D: np.log(D00 / D) - np.pi * a / (4 * D), a, D00)


print(f"\nAG, gapped branch to alpha/alpha_c = {a_gapless / a_c:.3f}")
print(" a/a_c   Tc/Tc0   eta(a) order-param   eta(b) spectral")
rows = []
for frac in np.linspace(0, 0.91, 92):
    a = frac * a_c
    t, d = Tc(a), Delta0(a)
    om = d * (1 - (a / d) ** (2 / 3)) ** 1.5
    rows.append((frac, t, d / (D00 * t), om / (D00 * t)))
for r in rows[::10]:
    print(f" {r[0]:5.2f}   {r[1]:.4f}   {r[2]:.4f}               {r[3]:.4f}")
ea = np.array([r[2] for r in rows]); eb = np.array([r[3] for r in rows])
print(f"\n(a) order parameter: min eta = {ea.min():.4f} -> eta < 1 unreachable on this branch")
cross = [r for r in rows if r[3] <= 0.38]
print(f"(b) spectral gap:    eta falls monotonically: {bool(np.all(np.diff(eb) < 0))}; "
      f"eta = 0.38 first reached at alpha/alpha_c ~ {cross[0][0]:.2f} (Tc/Tc0 = {cross[0][1]:.3f})")
print("\nVerdict: 'AG lowers T_c' does not decide the identity. Under (a) eta = 0.38 is outside AG;"
      " under (b) AG reaches it. Which Delta the formula uses is not stated in the paper.")
