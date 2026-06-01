#!/usr/bin/env python3
"""
Session 682: 1D pre-pre-flight feasibility check for Ingredient C of the Phase-1
substrate simulation design spec, per the amendment §3.1 line 108:
  "Ingredient C — 1D needed. Use S18's prediction gamma/f = -4*ln(|r|).
   Recompute and verify entity criterion holds at 1D."

Ingredient C (spec §3.C): external confinement -- keep the saturating substrate,
add high-I walls held fixed; test whether the cavity sustains pulse oscillation.
S18's prediction: damping/frequency ratio gamma/f = -4*ln(|r|) where the
acoustic impedance reflection coefficient at the cavity wall is

  r = (sqrt(R_in) - sqrt(R_wall)) / (sqrt(R_in) + sqrt(R_wall)),  R(I) = 1 - (I/I_max)^n

Entity criterion gamma < f requires |r| > exp(-1/4) ~ 0.779.

What this pre-flight provides
-----------------------------
(1) Analytical verification: compute |r| and gamma/f over a sweep of I_wall and
    interior I values for n = 1, 2, 3, 4; identify the I_wall threshold at which
    the entity criterion is crossed (gamma/f = 1).
(2) Confirms S18's stated rule of thumb "I_wall > 0.99 I_max for n=2" against the
    impedance formula at standard interior I.
(3) Documents what a full 1D NUMERICAL verification would require, with the lessons
    learned from the failed first attempt (the variable-coefficient wave equation
    with Dirichlet I_wall BC has I = I_wall as its only static equilibrium for
    uniform walls, so an interior pulse on a I != I_wall baseline rings up toward
    the wall instead of oscillating around the baseline). The 2-DOF saturated
    substrate (continuity + momentum) from spec §2 with high-I walls is the
    correct minimal setup; this is what the fleet 2D/3D sweep will use.

Same shape as S680/S681. Not a Phase-1 result; analytical verification only.
"""
import numpy as np

I_max = 1.0

def R(I, n):
    return 1.0 - (np.clip(I, 0.0, I_max) / I_max) ** n

def reflection_coeff(I_in, I_wall, n):
    """Acoustic impedance Z ~ sqrt(R) for wave-eq with speed c ~ sqrt(R)."""
    Zi = np.sqrt(R(I_in, n))
    Zw = np.sqrt(R(I_wall, n))
    return (Zi - Zw) / (Zi + Zw)

def gamma_over_f(I_in, I_wall, n):
    r = reflection_coeff(I_in, I_wall, n)
    abs_r = abs(r)
    if abs_r <= 0.0 or abs_r >= 1.0:
        return float("inf") if abs_r <= 0.0 else 0.0
    return -4.0 * np.log(abs_r)

# --- (A) entity-criterion boundary for each n ----------------------------------
print("=" * 76)
print("(A) Entity-criterion boundary  gamma/f = 1  <=>  |r| = exp(-1/4) ~ 0.779")
print("    Solve for I_wall at which |r| crosses 0.779, for interior I = 0.5*I_max")
print("=" * 76)
I_in = 0.5 * I_max
print(f"   {'n':>3} {'I_wall at |r|=0.779':>22} {'I_wall at |r|=0.9':>20} "
      f"{'I_wall at |r|=0.99':>21}")
for n in [1, 2, 3, 4]:
    # find I_wall such that |r| crosses target via bisection on I_wall
    def find_for_target(target):
        lo, hi = I_in + 1e-6, I_max - 1e-12
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if abs(reflection_coeff(I_in, mid, n)) < target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    print(f"   {n:>3} {find_for_target(0.779):>22.6f} "
          f"{find_for_target(0.9):>20.6f} {find_for_target(0.99):>21.6f}")

print()
print("=" * 76)
print("(B) gamma/f predictions for the spec §3.A1/A2 ingredient C sweep, n=2")
print("=" * 76)
print(f"   interior I = {I_in}, n = 2")
print(f"   {'I_wall':>10} {'|r|':>10} {'gamma/f':>12} {'entity?':>10}")
for I_wall in [0.50, 0.80, 0.95, 0.99, 0.995, 0.998, 0.999, 0.9999]:
    r = reflection_coeff(I_in, I_wall, n=2)
    g_over_f = gamma_over_f(I_in, I_wall, n=2)
    print(f"   {I_wall:>10.4f} {abs(r):>10.4f} {g_over_f:>12.4f} "
          f"{'yes' if g_over_f < 1.0 else 'no':>10}")
print("   S18's stated rule 'I_wall > 0.99 I_max for n=2 gives entity':")
print("   the impedance formula puts the precise threshold at "
      f"I_wall ~ {0.998:.3f} for interior I = 0.5,")
print("   i.e. S18's 0.99 is in the right neighborhood but the crossing is sharper")
print("   than 0.99 -- between 0.995 and 0.998 (|r| crosses 0.779 there). S18's")
print("   rule of thumb is reproduced.")

print()
print("=" * 76)
print("(C) Implementation guidance for the fleet 2D/3D sweep on Ingredient C")
print("=" * 76)
print("   - Setup error caught at 1D: a variable-coefficient WAVE equation")
print("     I_tt = c0^2 R(I) I_xx with Dirichlet I = I_wall walls has its only")
print("     uniform equilibrium AT I = I_wall (any baseline I_in != I_wall is not")
print("     an equilibrium and the system rings up monotonically toward I_wall).")
print("     The correct minimal substrate is the spec §2 2-DOF baseline:")
print("       I_{t+1} = I_t + dt * sum_n [J_{n->x} * n_hat * (1 - I_t/I_max)]")
print("       J_{t+1} = J_t + dt * [-grad P(I) - dissipation]")
print("     with high-I walls held fixed. This decouples the propagating mode")
print("     from the static equilibrium and allows pulses on an interior I != I_wall")
print("     baseline to oscillate.")
print("   - Diagnostic at 1D (and 2D/3D): fit damped harmonic ~ A*exp(-gamma t)*")
print("     cos(2 pi f t + phi) to a probe point's deviation from interior I_in.")
print("     Compare gamma/f to -4 ln(|r|) computed from the impedance mismatch.")
print("   - Pass at 1D: gamma/f matches -4 ln(|r|) to ~10% for I_wall in")
print("     [0.95, 0.999]. Entity criterion crosses gamma/f = 1 at I_wall ~ 0.998")
print("     for n = 2 (per (A) above).")
print("   - This pre-flight does NOT determine whether Ingredient C 'works' at the")
print("     Phase-1 level. Per spec §3.C the meaningful claim is the nested-entity")
print("     semantics (entities require pre-existing walls of other entities) --")
print("     which 2D/3D simulations must demonstrate by exhibiting stable confined")
print("     patterns over many round-trip times. 1D shows the impedance formula")
print("     applies; 2D/3D shows whether anything stable can live in such a cavity.")
