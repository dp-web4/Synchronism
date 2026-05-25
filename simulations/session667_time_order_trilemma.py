#!/usr/bin/env python3
"""
Session 667: The time-order trilemma — why the substrate findings keep recurring.

S665: the transfer rule's induced velocity is irrotational -> no vortices (spatial).
S666: the transfer rule is dissipative -> no oscillation; Schrödinger needs an
      inserted i (temporal).
This session asks the organizing question: WHY do these keep happening? The common
root is the ORDER-IN-TIME of the evolution equation. The transfer rule is
first-order-in-time and parabolic. That single fact forces THREE consequences, and
the three things the framework needs map to three INCOMPATIBLE PDE classes:

  Demand (FUNDAMENTALS)                  | Requires                | PDE class
  ---------------------------------------|-------------------------|-------------
  Saturation / arrow of time (Found. 3)  | dissipative term        | parabolic (1st order, real)
  Stable de Broglie entities (Core Def.) | UNDAMPED oscillation    | unitary / wave
  c emerges, light cone (substrate)      | finite propagation speed| hyperbolic (2nd order)

This script demonstrates the NEW leg (causality) and ties the three together:

 (A) The DISCRETE transfer rule respects a strict light cone: a disturbance reaches
     cell d no earlier than tick d  (speed = 1 cell/tick = ell_P/t_P = c). FINITE.
 (B) The CONTINUUM limit the framework actually uses (the diffusion PDE) is ACAUSAL:
     the heat kernel is nonzero everywhere instantly -> infinite propagation speed.
     So "c = one cell per tick" holds for the discrete rule and is DESTROYED by the
     continuum limit used for Schrödinger / N-S / C(rho).
 (C) The only single equation that is causal AND dissipative AND oscillatory is the
     TELEGRAPHER (damped wave) equation -- and its oscillation is DAMPED, so it gives
     no stable entities (matches the framework's own 2-DOF result, Sessions 17-18:
     "damped oscillation, amplitude decays"). Undamped-oscillation + dissipation is a
     contradiction in one linear equation.
"""
import numpy as np

I_MAX = 1.0
N_EXP = 2
def R(I):
    return 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** N_EXP

# ---------------------------------------------------------------------------
print("=" * 72)
print("(A) DISCRETE transfer rule: strict light cone (finite speed = 1 cell/tick)")
print("=" * 72)
n = 401
c0 = n // 2
bg = 0.10
k = 0.20                      # coupling
def discrete_step(I):
    # nearest-neighbor saturation-limited transfer: dI = k * sum_nb (I_nb - I) * R(I)
    lap = (np.roll(I, 1) + np.roll(I, -1) - 2 * I)
    return I + k * lap * R(I)
# Causality = nothing propagates beyond the cone: after t ticks every cell beyond
# c0+t must be EXACTLY background (nearest-neighbor coupling moves influence at most
# 1 cell/tick). (The edge cell's own amplitude ~k^t underflows into bg in float64
# for large t -- 0.10 + 1e-36 == 0.10 -- so we test the causal bound, not the edge.)
print(f"   {'tick t':>8} {'edge reached (cell c0+t)':>26} {'ALL cells > c0+t are bg?':>28}")
ok = True
for t in [5, 20, 50, 100, 150]:
    I = np.full(n, bg); I[c0] = 0.90
    for _ in range(t):
        I = discrete_step(I)
    edge_reached = abs(I[c0 + t] - bg) > 0.0            # resolvable only at small t
    beyond_all_bg = np.all(I[c0 + t + 1:] == bg)        # the causality bound
    ok = ok and beyond_all_bg
    edge_str = "yes" if edge_reached else "(<eps, underflow)"
    print(f"   {t:>8} {edge_str:>26} {str(bool(beyond_all_bg)):>28}")
print(f"   causal bound (no disturbance beyond c0+t) holds for all ticks: {ok}")
print("   => influence reaches cell d no sooner than tick d, never beyond: speed = 1")
print("      cell/tick = ell_P/t_P = c. The DISCRETE rule is causal (has a light cone).")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(B) CONTINUUM limit (heat kernel): infinite propagation speed (ACAUSAL)")
print("=" * 72)
# True dt->0 solution of dI/dt = D d2I/dx2 from a point source is the Gaussian
# kernel G(x,t)=exp(-x^2/4Dt)/sqrt(4 pi D t) > 0 for ALL x at any t>0.
# Compare the signal at a point OUTSIDE the discrete light cone (x = 2 c t).
D = 1.0
for tphys in [1.0, 5.0, 25.0]:
    # discrete cone edge sits at x_cone = c * t = 1*t (cells); probe at x = 2*t (outside)
    x_out = 2.0 * tphys
    G = np.exp(-x_out**2 / (4 * D * tphys)) / np.sqrt(4 * np.pi * D * tphys)
    print(f"   t={tphys:>5}:  discrete cone edge at x={tphys:>5.0f};  "
          f"heat-kernel value at x={x_out:>5.0f} (OUTSIDE cone) = {G:.3e}  (> 0!)")
print("   => the continuum solution is strictly positive beyond the light cone for")
print("      every t>0: exponentially small but NONZERO. Parabolic = infinite speed.")
print("      The framework computes Schrödinger / N-S / C(rho) on THIS acausal form.")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(C) TELEGRAPHER (damped wave): the only causal+dissipative single equation")
print("    -- but its oscillation is DAMPED (no stable entities)")
print("=" * 72)
# u_tt + (1/tau) u_t = v^2 u_xx  : finite speed v, dissipation 1/tau, damped waves.
# Explicit leapfrog. Show (1) a sharp finite-speed front (zero beyond x=v t), and
# (2) a probed oscillation that DECAYS.
m = 1201
mid = m // 2
dx = 1.0
v = 1.0
tau = 40.0
dt = dx / v                        # CFL = 1: numerical cone coincides with physical v
u = np.zeros(m); u_prev = np.zeros(m)
u[mid] = 1.0                       # initial localized kick
u_prev[mid] = 1.0
def lap1d(f):
    return (np.roll(f, 1) + np.roll(f, -1) - 2 * f) / dx**2
probe = mid + 25
series = []
front_zero_beyond = True
steps = 600
for t in range(steps):
    # damped-wave leapfrog: u_tt + (1/tau)u_t = v^2 u_xx
    u_next = (dt**2 * v**2 * lap1d(u) + 2*u - (1 - dt/(2*tau))*u_prev) / (1 + dt/(2*tau))
    u_prev, u = u, u_next
    series.append(u[probe])
    # causality: at CFL=1 the exact front is at x = (t+1) cells; beyond must stay 0
    edge = (t + 1) + 1
    if mid + edge < m and np.any(u[mid + edge:] != 0.0):
        front_zero_beyond = False
series = np.array(series)
# count sign changes (oscillation) and check amplitude decay (damping)
sign_changes = int(np.sum(np.diff(np.sign(series[series != 0])) != 0))
first_peak = np.max(np.abs(series[:len(series)//3])) if len(series) > 3 else 0
last_peak = np.max(np.abs(series[2*len(series)//3:])) if len(series) > 3 else 0
print(f"   finite-speed front (zero beyond x = v t): {front_zero_beyond}")
print(f"   probe oscillation: {sign_changes} sign changes (OSCILLATES)")
print(f"   amplitude early {first_peak:.4e} -> late {last_peak:.4e}  "
      f"(ratio {last_peak/max(first_peak,1e-30):.3f}) => DAMPED" )
print("   => causal + dissipative + oscillatory, but the oscillation DECAYS.")
print("      Stable entities need UNDAMPED oscillation; damping = dissipation of it.")
print("      This is exactly the framework's own 2-DOF result (Sessions 17-18).")

print()
print("=" * 72)
print("TRILEMMA (no single linear evolution equation satisfies all three):")
print("=" * 72)
print("   equation            | saturation/arrow | undamped oscillation | finite-speed cone")
print("   transfer rule (par.)|        YES       |         no           |        no")
print("   Schrödinger (unit.) |        no        |        YES           |        no (inf. speed)")
print("   wave (hyperbolic)   |        no        |        YES           |       YES")
print("   telegrapher (damped)|        YES       |     no (damped)      |       YES")
print("   -> the two the framework most needs (saturation + UNDAMPED oscillation) are")
print("      never both satisfied: an undamped oscillator conserves energy; a")
print("      dissipative one cannot oscillate forever. Causality is a third miss.")
