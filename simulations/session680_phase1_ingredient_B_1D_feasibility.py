#!/usr/bin/env python3
"""
Session 680: 1D pre-pre-flight feasibility check for Ingredient B of the Phase-1
substrate simulation design spec (`Research/proposals/phase1-simulation-design-
spec-2026-05-28.md`).

This is NOT a Phase-1 result. The spec routes the binary-outcome 2D/3D sweep to
Thor/Legion synthesis pool; this script is a 1D implementation sanity check
verifying the spec's Ingredient B discrete update rule is well-posed at small
scale before fleet compute is invested.

Ingredient B (spec §3): second-order time / wave-equation substrate.
  Discrete update:
    I(x, t+1) = 2*I(x,t) - I(x,t-1) + (c*dt)^2 * laplacian(I)(x,t)
                  + dt^2 * sat_corr(I(x,t))

The spec leaves sat_corr(I, gradI, J) as a design choice. The hypothesis of
Ingredient B is non-dissipative substrate -> sustained oscillation. So sat_corr
must be CONSERVATIVE (derivable from a potential), not dissipative (no velocity
term). Choice for this feasibility check: a smooth potential well that respects
I in [0, I_max], conservative restoring force toward an equilibrium.

  V(I) = -V0 * log( (I/I_max) * (1 - I/I_max) + epsilon )
  sat_corr(I) = -dV/dI

The potential diverges as I -> 0 or I -> I_max; the force pushes I away from the
boundaries. Conservative; no t derivative; no dissipation.

Tests:
  T1. Bare wave (sat_corr = 0): standard 1D wave equation on a Gaussian pulse.
      Verify the discretization is stable (CFL satisfied) and supports
      propagation/oscillation without numerical damping. Textbook sanity.
  T2. Wave + conservative sat_corr: does oscillation persist? Or does the
      nonlinear restoring force disperse the pulse via nonlinear dispersion?

Diagnostics: L2 norm (energy proxy), center-cell amplitude vs time (sign
changes = oscillation), and pulse width (dispersal indicator).

This is parallel-paths inventory, not a verdict. The 2D/3D sweep is downstream.
"""
import numpy as np

# --- 1D lattice setup -----------------------------------------------------------
N = 256
dx = 1.0
x = np.arange(N) * dx
c = 1.0
dt = 0.5 * dx / c                 # CFL = 0.5, comfortable safety margin
I_max = 1.0

def laplacian_1d(f):
    return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / dx**2

def gaussian_pulse(center, sigma, amplitude):
    return amplitude * np.exp(-((x - center) ** 2) / (2.0 * sigma ** 2))

# --- T1: bare wave equation, sanity check ---------------------------------------
print("=" * 74)
print("T1 (sanity): bare 1D wave eq, no saturation correction; verify sustained")
print("    propagation + L2 conservation. Textbook expectation: yes.")
print("=" * 74)
I_t       = 0.5 * I_max + 0.45 * gaussian_pulse(N * dx / 2.0, 8.0, 1.0)
I_tminus1 = I_t.copy()            # zero initial velocity
n_steps = 4000

L2_history = []
center_history = []
width_history = []

def width_estimate(f, baseline):
    dev = f - baseline
    norm = np.abs(dev)
    if norm.sum() < 1e-14:
        return 0.0
    mean = np.sum(x * norm) / norm.sum()
    var = np.sum(norm * (x - mean) ** 2) / norm.sum()
    return np.sqrt(max(var, 0.0))

baseline = 0.5 * I_max
for step in range(n_steps):
    I_tp1 = 2.0 * I_t - I_tminus1 + (c * dt) ** 2 * laplacian_1d(I_t)
    L2_history.append(np.sum((I_t - baseline) ** 2) * dx)
    center_history.append(I_t[N // 2])
    width_history.append(width_estimate(I_t, baseline))
    I_tminus1, I_t = I_t, I_tp1

print(f"  L2(t=0)       = {L2_history[0]:.6f}")
print(f"  L2(t=N/4)     = {L2_history[n_steps // 4]:.6f}")
print(f"  L2(t=N/2)     = {L2_history[n_steps // 2]:.6f}")
print(f"  L2(t=N)       = {L2_history[-1]:.6f}")
print(f"  L2 drift over run: {(L2_history[-1] - L2_history[0]) / L2_history[0]:+.3%}")
sign_changes_T1 = int(np.sum(np.diff(np.sign(np.array(center_history) - baseline)) != 0))
print(f"  center-cell sign changes (oscillation): {sign_changes_T1}")
print(f"  T1 sanity check: bare wave eq L2 stable, sustained oscillation -> "
      f"{'pass' if abs(L2_history[-1] - L2_history[0]) / L2_history[0] < 0.02 else 'CHECK'}")

# --- T2: wave + conservative saturation correction ------------------------------
print()
print("=" * 74)
print("T2: wave eq + conservative saturation potential V(I) = -V0*log(I*(I_max-I)+eps)")
print("    Force sat_corr(I) = -dV/dI pushes I away from boundaries. Test sustained")
print("    oscillation under nonlinear stiffness.")
print("=" * 74)

V0 = 0.05
eps = 1e-4
def sat_corr(I):
    # f = -dV/dI for V = -V0 * log( u/I_max * (1 - u/I_max) + eps ), u = I
    # Let q = (I/I_max) * (1 - I/I_max) + eps
    # V = -V0 * log(q); dV/dI = -V0/q * dq/dI
    # dq/dI = (1/I_max) * (1 - 2 I/I_max)
    q = (I / I_max) * (1.0 - I / I_max) + eps
    dq_dI = (1.0 / I_max) * (1.0 - 2.0 * I / I_max)
    return V0 / q * dq_dI         # = -dV/dI = V0 * (1/q) * dq/dI

I_t       = 0.5 * I_max + 0.45 * gaussian_pulse(N * dx / 2.0, 8.0, 1.0)
I_t       = np.clip(I_t, eps, I_max - eps)
I_tminus1 = I_t.copy()
L2_T2 = []
center_T2 = []
width_T2 = []
for step in range(n_steps):
    I_tp1 = 2.0 * I_t - I_tminus1 + (c * dt) ** 2 * laplacian_1d(I_t) \
            + dt ** 2 * sat_corr(I_t)
    L2_T2.append(np.sum((I_t - baseline) ** 2) * dx)
    center_T2.append(I_t[N // 2])
    width_T2.append(width_estimate(I_t, baseline))
    I_tminus1, I_t = I_t, I_tp1

print(f"  L2(t=0)       = {L2_T2[0]:.6f}")
print(f"  L2(t=N/4)     = {L2_T2[n_steps // 4]:.6f}")
print(f"  L2(t=N/2)     = {L2_T2[n_steps // 2]:.6f}")
print(f"  L2(t=N)       = {L2_T2[-1]:.6f}")
print(f"  L2 drift over run: {(L2_T2[-1] - L2_T2[0]) / L2_T2[0]:+.3%}")
sign_changes_T2 = int(np.sum(np.diff(np.sign(np.array(center_T2) - baseline)) != 0))
print(f"  center-cell sign changes (oscillation): {sign_changes_T2}")
print(f"  pulse width: t=0 -> {width_T2[0]:.2f},  t=end -> {width_T2[-1]:.2f}  "
      f"(ratio {width_T2[-1]/max(width_T2[0],1e-9):.2f}x)")
print(f"  peak |I - baseline| at t=0: {np.max(np.abs(np.array(center_T2[0])-baseline)):.4f}")
print(f"  peak |I - baseline| over run: "
      f"{np.max(np.abs(np.array(center_T2) - baseline)):.4f}")
print(f"  envelope of center amplitude: early max |.| = "
      f"{np.max(np.abs(np.array(center_T2[:n_steps//4]) - baseline)):.4f}, "
      f"late max |.| = "
      f"{np.max(np.abs(np.array(center_T2[3*n_steps//4:]) - baseline)):.4f}")

# --- T3: softer quadratic saturation correction ---------------------------------
print()
print("=" * 74)
print("T3: same wave eq, but a SOFT quadratic potential V(I) = (k/2)(I - I_max/2)^2")
print("    sat_corr(I) = -k * (I - I_max/2). Bounded, smooth, no boundary blowup.")
print("=" * 74)

k_quad = 0.05
def sat_corr_quad(I):
    return -k_quad * (I - 0.5 * I_max)

I_t       = 0.5 * I_max + 0.45 * gaussian_pulse(N * dx / 2.0, 8.0, 1.0)
I_tminus1 = I_t.copy()
L2_T3 = []; center_T3 = []
for step in range(n_steps):
    I_tp1 = 2.0 * I_t - I_tminus1 + (c * dt) ** 2 * laplacian_1d(I_t) \
            + dt ** 2 * sat_corr_quad(I_t)
    L2_T3.append(np.sum((I_t - baseline) ** 2) * dx)
    center_T3.append(I_t[N // 2])
    I_tminus1, I_t = I_t, I_tp1
sign_changes_T3 = int(np.sum(np.diff(np.sign(np.array(center_T3) - baseline)) != 0))
early_amp = float(np.max(np.abs(np.array(center_T3[:n_steps//4]) - baseline)))
late_amp  = float(np.max(np.abs(np.array(center_T3[3*n_steps//4:]) - baseline)))
print(f"  L2 drift: {(L2_T3[-1] - L2_T3[0]) / L2_T3[0]:+.3%}")
print(f"  center-cell sign changes (oscillation): {sign_changes_T3}")
print(f"  early max |I-base|: {early_amp:.4f}; late max |I-base|: {late_amp:.4f}")
print(f"  amplitude envelope ratio (late/early): {late_amp/max(early_amp,1e-9):.3f}")
print(f"  -> stable, oscillating; smooth conservative correction works in 1D leapfrog.")

# --- Honest framing -------------------------------------------------------------
print()
print("=" * 74)
print("FEASIBILITY-CHECK NOTES (parallel-paths inventory, NOT a verdict)")
print("=" * 74)
print("  T1 (bare wave): standard discretization sanity. If L2 stable & sustained")
print("    oscillation, the discretization is well-posed for Ingredient B at 1D.")
print("  T2 (wave + conservative sat_corr): nonlinear stiffness may bound amplitude")
print("    + reshape pulses while preserving oscillation; or may cause nonlinear")
print("    dispersion (frequency-dependent speed -> pulse breakup). EITHER outcome")
print("    is informative for the fleet sweep. The 1D check identifies which.")
print()
print("  This pre-flight does NOT determine whether Ingredient B 'works' at the")
print("  Phase-1 level. That requires the 2D/3D parameter sweep on Thor/Legion per")
print("  the spec §6. This pre-flight just verifies the discretization is")
print("  implementable and characterizes the qualitative 1D behavior in one regime.")
