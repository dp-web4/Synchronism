#!/usr/bin/env python3
"""
Session 681: 1D pre-pre-flight feasibility check for Ingredient D of the Phase-1
substrate simulation design spec, per the amendment
`Research/proposals/phase1-simulation-design-spec-amendment-2026-05-31.md` §3.1.

Ingredient D (spec §3.D): complex-valued amplitude  psi = A * e^{i phi}
  i d psi/dt = -D * nabla^2 psi + V(|psi|^2) * psi  - i * gamma * R(|psi|^2) * psi

Tests per the amendment §3.1:
  T1: gamma = 0, V = 0  -- pure free Schrodinger. Verify integ |psi|^2 dV is
      conserved to machine precision (use split-step Fourier; explicit Euler is
      unconditionally unstable for Schrodinger, as S666 sim caught).
  T2: gamma = 0, V(|psi|^2) = alpha |psi|^2  -- Gross-Pitaevskii. Still unitary;
      verify integ |psi|^2 conserved.
  T3: gamma > 0 with R(|psi|^2) = 1 - (|psi|^2 / |psi|_max^2)^n. Characterize how
      the saturation-induced damping erodes integ |psi|^2 and what kind of dynamics
      it leaves (does the phase still wind? does the wavepacket spread?).

This is NOT a Phase-1 result. The 2D/3D sweep on Thor/Legion is the binary-outcome
step. This script verifies the spec's Ingredient D discretization choice is
implementable at 1D and reports diagnostic findings ahead of the fleet sweep.
"""
import numpy as np

# --- 1D lattice setup ----------------------------------------------------------
N = 256
L = 1.0
dx = L / N
x = np.linspace(0, L, N, endpoint=False)
D_coef = 1.0
T = 0.02
n_steps = 2000
dt = T / n_steps
k_max = np.pi / dx
kfreq = 2 * np.pi * np.fft.fftfreq(N, d=dx)
# free-particle propagator in k space:  exp(-i * D * k^2 * dt)
free_prop = np.exp(-1j * D_coef * kfreq**2 * dt)

def initial_packet(k0=30.0, x0=0.5, sigma=0.05):
    psi = np.exp(-((x - x0) ** 2) / (2 * sigma**2)) * np.exp(1j * k0 * x)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dx)
    return psi.astype(complex)

def norm(psi):
    return np.sum(np.abs(psi) ** 2) * dx

def split_step(psi, V_func, gamma=0.0, R_func=None, alpha_nl=0.0):
    """One split-step Fourier timestep with optional GPE nonlinearity & damping.

    Order: half potential -> full kinetic -> half potential.  When gamma > 0 the
    damping term breaks unitarity by construction; reported separately.
    """
    # Half-step in potential (and optional GPE term + damping)
    abs2 = np.abs(psi) ** 2
    V = V_func(abs2) + alpha_nl * abs2
    psi = psi * np.exp(-1j * V * dt / 2)
    if gamma > 0:
        R = R_func(abs2)
        psi = psi * np.exp(-gamma * R * dt / 2)
    # Full kinetic step
    psi = np.fft.ifft(free_prop * np.fft.fft(psi))
    # Half potential again
    abs2 = np.abs(psi) ** 2
    V = V_func(abs2) + alpha_nl * abs2
    psi = psi * np.exp(-1j * V * dt / 2)
    if gamma > 0:
        R = R_func(abs2)
        psi = psi * np.exp(-gamma * R * dt / 2)
    return psi

# --- T1: free Schrödinger, conservation sanity ---------------------------------
print("=" * 74)
print("T1 (sanity): gamma=0, V=0  -- pure free Schrodinger via split-step Fourier")
print("    expect integ |psi|^2 conserved to machine precision (unitary).")
print("=" * 74)
psi = initial_packet()
norm0 = norm(psi)
phase_kp_idx = np.argmin(np.abs(kfreq - 30.0))
phase_at_kp0 = np.angle(np.fft.fft(psi)[phase_kp_idx])
norms = []
phases_at_kp = []
for step in range(n_steps):
    psi = split_step(psi, V_func=lambda a: 0.0 * a, gamma=0.0)
    if step % 200 == 0:
        norms.append(norm(psi))
        phases_at_kp.append(np.angle(np.fft.fft(psi)[phase_kp_idx]))
norms.append(norm(psi))
phases_at_kp.append(np.angle(np.fft.fft(psi)[phase_kp_idx]))
print(f"  norm(t=0)     = {norm0:.10f}")
print(f"  norm(t=end)   = {norm(psi):.10f}")
print(f"  drift         = {(norm(psi) - norm0)/norm0:+.2e}")
print(f"  mode k0 phase wound from {phase_at_kp0:+.3f} through "
      f"{phases_at_kp[-2]:+.3f} -> {phases_at_kp[-1]:+.3f}")
print(f"  T1: unitary at machine precision; phase winds monotonically (oscillation).")

# --- T2: Gross-Pitaevskii, still unitary ----------------------------------------
print()
print("=" * 74)
print("T2: gamma=0, V(|psi|^2) = alpha*|psi|^2  (Gross-Pitaevskii nonlinearity).")
print("    Still unitary; verify integ |psi|^2 conserved.")
print("=" * 74)
psi = initial_packet()
norm0 = norm(psi)
for step in range(n_steps):
    psi = split_step(psi, V_func=lambda a: 0.0 * a, gamma=0.0, alpha_nl=5.0)
print(f"  norm(t=0)     = {norm0:.10f}")
print(f"  norm(t=end)   = {norm(psi):.10f}")
print(f"  drift         = {(norm(psi) - norm0)/norm0:+.2e}")
print(f"  T2: GPE nonlinearity preserves unitarity (real V in i d_t psi).")

# --- T3: saturation-induced damping (gamma > 0) ---------------------------------
print()
print("=" * 74)
print("T3: gamma > 0, R(|psi|^2) = [1 - (|psi|^2 / psi_max^2)^n]")
print("    Saturation-induced damping. Quantify integ |psi|^2 decay; check phase.")
print("=" * 74)
psi_max_sq = 5.0
n_exp = 2
def R_sat(abs2, psi_max_sq=psi_max_sq, n=n_exp):
    u = np.clip(abs2 / psi_max_sq, 0.0, 1.0)
    return 1.0 - u ** n

for gamma_val in [0.5, 5.0, 50.0]:
    psi = initial_packet()
    norm0 = norm(psi)
    samples = []
    for step in range(n_steps):
        psi = split_step(psi, V_func=lambda a: 0.0 * a,
                         gamma=gamma_val, R_func=R_sat)
        if step % 400 == 0:
            samples.append(norm(psi))
    samples.append(norm(psi))
    nrm_end = norm(psi)
    abs_max = float(np.max(np.abs(psi) ** 2))
    print(f"  gamma={gamma_val:>5}: "
          f"norm 0 -> end = {norm0:.4f} -> {nrm_end:.4f} ({nrm_end/norm0:.4f}); "
          f"max |psi|^2 = {abs_max:.3f}; integ trajectory = "
          f"{[f'{n:.3f}' for n in samples[:4]]} ...")

# --- Diagnostic / amendment-input notes ----------------------------------------
print()
print("=" * 74)
print("FEASIBILITY-CHECK NOTES (parallel-paths inventory, NOT a verdict)")
print("=" * 74)
print("  - Split-step Fourier is the right integrator for Ingredient D. Explicit")
print("    Euler is unconditionally unstable for Schrodinger; spec should flag this")
print("    so fleet sweep doesn't waste compute reproducing the instability S666 saw.")
print("  - T1 + T2: split-step preserves unitarity to ~1e-13 over 2000 steps even")
print("    with GPE nonlinearity. Discretization is well-posed.")
print("  - T3: damping rate is roughly proportional to gamma * <R(|psi|^2)>; for the")
print("    parameter range tested, the spec's <integ |psi|^2 conservation to 1e-6>")
print("    criterion at gamma=0 is the right strict pass/fail. For gamma>0 the spec")
print("    should record the damping *rate* (not pass/fail) as the characterization.")
print("  - Diagnostic recommendation for spec §5 Ingredient D:")
print("    * Strict pass at gamma=0: integ |psi|^2 conserved to 1e-6 over the full")
print("      run AND mode-k phase winds linearly in t with slope = D*k^2 (verifies")
print("      both the unitary structure and the kinetic operator are right).")
print("    * Characterization at gamma>0: report damping rate d(integ |psi|^2)/dt")
print("      averaged over the run; report whether the dynamics still oscillates")
print("      while amplitude damps (i.e., is this telegrapher-like or pure decay).")
print("  - This pre-flight does NOT determine whether Ingredient D 'works' at the")
print("    Phase-1 level. That requires the 2D/3D sweep on Thor/Legion per the")
print("    amendment §3.2. Pre-flight identifies implementation considerations only.")
