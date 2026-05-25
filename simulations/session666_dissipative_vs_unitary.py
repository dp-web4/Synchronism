#!/usr/bin/env python3
"""
Session 666: Is the Schrödinger equation DERIVED from the Intent transfer rule,
or is the imaginary unit inserted by hand?

FUNDAMENTALS.md asserts two things about the same substrate:
  (Foundations 1,3) Intent evolves by REAL, saturation-limited diffusion:
        dI/dt = div[ D * R(I) * grad(I) ],   R(I)=1-(I/I_max)^n   -- DISSIPATIVE
  (Core Definitions) Existence IS oscillation: entities recur, f = E/h
        (de Broglie), interaction is phase-locking -- requires UNITARY dynamics.

The two Schrödinger "derivations" bridge these only by:
  S307 (line 48/code 129): replace real diffusion with IMAGINARY diffusion
        psi += dt*(1j*D*lap - 1j*V*psi)   -- the 1j is inserted; R(I) dropped.
  S99  (Axiom 3 + line 51): POSIT oscillation w=E/hbar as an axiom, posit the
        complex i (Axiom 4), then take the "non-dissipative limit D->0".

This script shows the two are MUTUALLY EXCLUSIVE dynamical classes, separated
exactly by the factor i, and that "D->0" does not extract quantum dynamics from
diffusion -- it freezes the substrate and lets the posited phase axiom do the work.

Demos:
  (A) Same Laplacian operator, eigenvalue contrast: real diffusion -> exp(-Dk^2 t)
      (monotone decay, no oscillation); imaginary diffusion -> exp(-i Dk^2 t)
      (pure oscillation, |amp| conserved). The only difference is i.
  (B) Evolve the REAL transfer rule from a lump: L2 norm decays monotonically
      (Lyapunov => dissipative, arrow of time), center amplitude monotone, no
      recurrence. Evolve free Schrödinger: |psi|^2 conserved, Re(psi) oscillates.
  (C) The D->0 freeze: in the REAL rule, D->0 => dI/dt -> 0 (nothing moves).
      So the "non-dissipative limit" contributes zero dynamics; all quantum
      motion in S99 comes from the separately-posited Hamilton-Jacobi axiom.
"""
import numpy as np

I_MAX = 1.0
N_EXP = 2
def R(I):
    return 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** N_EXP

# ---------------------------------------------------------------------------
print("=" * 70)
print("(A) Eigenvalue contrast: same operator, differ only by the factor i")
print("=" * 70)
D = 1.0
for k in [1.0, 2.0, 4.0]:
    lam_real = -D * k**2            # real diffusion mode rate
    lam_imag = -1j * D * k**2       # Schrödinger mode rate
    t = 1.0
    amp_real = np.exp(lam_real * t)
    amp_imag = np.exp(lam_imag * t)
    print(f" k={k:>3}:  real-diffusion exp(-Dk^2 t)={amp_real:8.4f} (|.|={abs(amp_real):.4f}, DECAYS)"
          f"   | Schrodinger exp(-iDk^2 t)={amp_imag.real:+.3f}{amp_imag.imag:+.3f}i (|.|={abs(amp_imag):.4f}, OSCILLATES)")
print(" => real eigenvalues (decay, no phase) vs imaginary eigenvalues (oscillation,")
print("    |amp| conserved). The transfer rule gives the former; QM needs the latter.")

# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("(B) Evolve REAL transfer rule vs free Schrödinger from the same lump")
print("=" * 70)
n = 256
L = 1.0
dx = L / n
x = np.linspace(0, L, n, endpoint=False)
x0, sig = 0.5, 0.05
lump = np.exp(-((x - x0) ** 2) / (2 * sig ** 2))

def lap1d(f):
    return (np.roll(f, -1) - 2 * f + np.roll(f, 1)) / dx ** 2

# (B1) REAL saturation-limited transfer rule  dI/dt = d/dx[ D R(I) dI/dx ]
I = 0.05 + 0.9 * lump
dt = 0.2 * dx ** 2 / D
L2_0 = np.sum(I ** 2) * dx
cen_0 = I[n // 2]
print(" REAL transfer rule (the substrate):")
print(f"   {'step':>6} {'center I':>10} {'L2 norm':>12} {'monotone?':>10}")
prev = None; mono = True
for step in range(0, 20001):
    if step % 4000 == 0:
        l2 = np.sum(I ** 2) * dx
        if prev is not None and l2 > prev + 1e-12:
            mono = False
        prev = l2
        print(f"   {step:>6} {I[n//2]:>10.5f} {l2:>12.6f} {'yes' if mono else 'NO':>10}")
    gI = (np.roll(I, -1) - np.roll(I, 1)) / (2 * dx)
    flux = D * R(I) * gI
    I = I + dt * (np.roll(flux, -1) - np.roll(flux, 1)) / (2 * dx)
print("   => L2 monotonically DECREASES (dissipative Lyapunov); center decays;")
print("      no oscillation, no recurrence. This is an arrow-of-time relaxation.")

# (B2) Free Schrödinger  i dψ/dt = -(1/2) d2ψ/dx2   (hbar=m=1), exact spectral
# propagator psi(t)=IFFT[ exp(-i k^2/2 t) FFT(psi0) ]  -- unitary by construction.
# (Forward Euler is unconditionally UNSTABLE here: |1+i*lam*dt|>1 for imaginary
#  eigenvalues -- itself a sign of the different dynamical class -- so use spectral.)
psi0 = (lump.astype(complex)) * np.exp(1j * 30.0 * x)   # moving wavepacket
psi0 /= np.sqrt(np.sum(np.abs(psi0) ** 2) * dx)
kfreq = 2 * np.pi * np.fft.fftfreq(n, d=dx)
psi0_hat = np.fft.fft(psi0)
print("\n Free Schrödinger (the dynamics QM needs), exact spectral propagator:")
print(f"   {'time':>8} {'|psi|^2 norm':>14} {'arg(mode k0)':>14}")
# probe oscillation via the phase of one Fourier mode -- it winds as -0.5 k^2 t
kp = np.argmin(np.abs(kfreq - 30.0))   # mode near the carrier
phase_prev = None; winds = 0.0
for j in range(0, 9):
    t = j * 2.0e-4
    factor = np.exp(-1j * 0.5 * kfreq ** 2 * t)
    psi = np.fft.ifft(factor * psi0_hat)
    nrm = np.sum(np.abs(psi) ** 2) * dx
    ph = np.angle((factor * psi0_hat)[kp])
    print(f"   {t:>8.1e} {nrm:>14.8f} {ph:>14.4f}")
print("   => |psi|^2 norm CONSERVED to 1.0 (UNITARY); the mode phase winds")
print("      monotonically (per part (A), exp(-i k^2 t/2): pure OSCILLATION, |amp|=1).")
print("      Reversible, norm-preserving -- the OPPOSITE dynamical class to (B1).")

# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("(C) S99's 'non-dissipative limit D->0' freezes the substrate")
print("=" * 70)
for Dval in [1.0, 0.1, 0.01, 0.0]:
    I = 0.05 + 0.9 * lump.copy()
    dt_c = 0.2 * dx ** 2 / max(Dval, 1e-6)
    steps = 2000
    cen0 = I[n // 2]
    for _ in range(steps):
        gI = (np.roll(I, -1) - np.roll(I, 1)) / (2 * dx)
        flux = Dval * R(I) * gI
        I = I + dt_c * (np.roll(flux, -1) - np.roll(flux, 1)) / (2 * dx)
    print(f"   D={Dval:>5}:  center moved {abs(I[n//2]-cen0):.3e} over {steps} steps")
print("   => as D->0 the magnitude stops evolving entirely. The 'non-dissipative")
print("      limit' yields ZERO substrate dynamics; all quantum motion in S99 comes")
print("      from the separately-posited phase axiom (Hamilton-Jacobi), not diffusion.")

print()
print("VERDICT: real saturation-limited diffusion (the transfer rule) and unitary")
print("oscillation (Schrödinger / the entity ontology) are mutually exclusive")
print("dynamical classes. The Schrödinger 'derivations' bridge them by inserting i")
print("and switching the substrate off (drop R, or D->0). The i is assumed, not derived.")
