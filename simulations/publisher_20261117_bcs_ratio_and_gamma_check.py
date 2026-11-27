"""Publisher check 2026-11-17: the paper's "BCS ratio derived: 2*sqrt(pi) = 3.54, observed 3.52 (<1% error)".

Printed at superconductivity.md (Derived Gap Ratio), conclusion.md (BCS ratio bullet, with a check mark)
and appendix_b_chemistry.md L88. Source: Research/Chemistry/Session1_Superconductivity_Coherence.md §3.3.

What this script does (no data, textbook BCS only):
  1. Solves the BCS gap equation (weak coupling, cutoff hbar*omega_D) for Delta(0) and T_c at several
     couplings lambda, and prints 2*Delta(0)/(k_B*T_c). Identity control: the closed form
     2*pi*exp(-gamma_Euler) = 3.5278.
  2. Compares that with 2*sqrt(pi) = 3.5449.
  3. Prints the archive's own "Synchronism" ratio formula from simulations/chemistry/superconductor_coherence.py
     (2*pi / tanh(gamma*ln 2), and the plotting line's 2*pi / tanh(gamma*ln 3)) at gamma = 2.

Units: k_B = 1, hbar*omega_D = 1. Run: python3 -I <this file>
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

EULER = 0.5772156649015329


def gap_eq_T0(delta, lam):
    # 1 = lam * int_0^1 de / sqrt(e^2 + d^2)  ->  asinh(1/d)
    return lam * np.arcsinh(1.0 / delta) - 1.0


def gap_eq_Tc(T, lam):
    # 1 = lam * int_0^1 tanh(e / 2T) / e de
    val, _ = quad(lambda e: np.tanh(e / (2 * T)) / e if e > 0 else 1 / (2 * T), 0, 1, limit=400)
    return lam * val - 1.0


closed = 2 * np.pi * np.exp(-EULER)
print(f"closed form 2*pi*exp(-gamma_E)          = {closed:.4f}")
print(f"2*sqrt(pi)                              = {2*np.sqrt(np.pi):.4f}")
print(f"difference                              = {100*(2*np.sqrt(np.pi)/closed-1):.2f}%")
print()
print("numerical BCS gap equation (weak coupling):")
for lam in [0.15, 0.2, 0.25, 0.3]:
    d0 = brentq(gap_eq_T0, 1e-12, 10, args=(lam,))
    tc = brentq(gap_eq_Tc, 1e-8, 5, args=(lam,), xtol=1e-14)
    print(f"  lambda={lam:.2f}: Delta0={d0:.4e}  Tc={tc:.4e}  2*Delta0/Tc={2*d0/tc:.4f}")
print()
print("archive's own 'Synchronism' ratio formulas at gamma = 2:")
for g in [2.0]:
    print(f"  2*pi/tanh(g*ln2) = {2*np.pi/np.tanh(g*np.log(2)):.3f}   (QUANTITATIVE TEST block)")
    print(f"  2*pi/tanh(g*ln3) = {2*np.pi/np.tanh(g*np.log(3)):.3f}   (plot block, assigned but not drawn)")
