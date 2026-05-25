#!/usr/bin/env python3
"""
Session 665: Does the Intent transfer rule actually give Navier-Stokes?

The CFD reframing (2026-03-08) claims: "The Navier-Stokes equations are not an
analogy applied to Intent dynamics. They ARE the Intent dynamics, in continuum
form." It identifies the velocity field as v = J/I = -D*R(I)*grad(I)/I.

This script stress-tests that claim on the one property that separates genuine
N-S from scalar diffusion: VORTICITY.

Claim under test (analytic): any velocity field of the form v = -g(I)*grad(I),
where g is any scalar function of the scalar field I, is irrotational
(curl v = 0) everywhere, for all time. If true, the "Intent fluid" can never
support vortices, vortex stretching, or a turbulent cascade -- which the CFD
document's downstream phenomenology (qualia as vortex modes, consciousness as
self-similar turbulence, dark matter as viscous vortex drag) all require.

Two numerical demonstrations:
  (A) Take a generic scalar field I(x,y), build v = -g(I)*grad(I), measure curl.
      Compare against a genuinely rotational field. Expect curl(v) ~ machine eps.
  (B) Evolve the ACTUAL scalar transfer rule  dI/dt = div[D*R(I)*grad(I)]  from a
      vortex-like (ring) initial condition. Track the induced velocity's vorticity
      and the structure's amplitude. Expect: vorticity stays ~0, ring disperses.
      (Independently consistent with grid Sessions 19-22: R(I) is defocusing.)
"""
import numpy as np

np.random.seed(0)

I_MAX = 1.0
N_EXP = 2          # R(I) = 1 - (I/I_max)^n
D = 1.0

def R(I):
    return 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** N_EXP

def grad(F, dx):
    gy, gx = np.gradient(F, dx, dx)   # np.gradient returns d/d(axis0), d/d(axis1)
    return gx, gy

def curl_z(vx, vy, dx):
    # 2D scalar vorticity  w = dvy/dx - dvx/dy
    dvy_dx = np.gradient(vy, dx, axis=1)
    dvx_dy = np.gradient(vx, dx, axis=0)
    return dvy_dx - dvx_dy

def divergence(vx, vy, dx):
    dvx_dx = np.gradient(vx, dx, axis=1)
    dvy_dy = np.gradient(vy, dx, axis=0)
    return dvx_dx + dvy_dy

def intent_velocity(I, dx):
    """v = J/I = -D*R(I)*grad(I)/I  (the CFD doc's identification, line 146)."""
    gx, gy = grad(I, dx)
    g = D * R(I) / np.maximum(I, 1e-9)
    return -g * gx, -g * gy

# ---------------------------------------------------------------------------
# (A) Irrotationality of the induced velocity for a generic scalar field
# ---------------------------------------------------------------------------
print("=" * 70)
print("(A) Curl of the Intent velocity v = -g(I) grad(I) for a generic field")
print("=" * 70)

n = 128
L = 1.0
dx = L / n
xs = np.linspace(0, L, n, endpoint=False)
X, Y = np.meshgrid(xs, xs)

# A lumpy, non-trivial scalar field (several Gaussian blobs), in (0, I_max)
I_field = 0.15 * np.ones((n, n))
for _ in range(6):
    cx, cy = np.random.rand(2) * L
    w = 0.05 + 0.1 * np.random.rand()
    I_field += 0.6 * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * w ** 2))
I_field = np.clip(I_field, 1e-3, 0.98 * I_MAX)

vx, vy = intent_velocity(I_field, dx)
w_intent = curl_z(vx, vy, dx)
speed = np.sqrt(vx ** 2 + vy ** 2)

# Reference: a genuinely rotational field of comparable magnitude
vx_rot = -(Y - L / 2)
vy_rot = (X - L / 2)
scale = speed.mean() / (np.sqrt(vx_rot ** 2 + vy_rot ** 2).mean() + 1e-12)
vx_rot *= scale; vy_rot *= scale
w_rot = curl_z(vx_rot, vy_rot, dx)

# dimensionless: |curl| * (length / typical speed)
typ = speed.mean() + 1e-12
print(f"Intent velocity:   mean|v| = {speed.mean():.4e}")
print(f"  mean|curl v| (interior)        = {np.abs(w_intent)[2:-2,2:-2].mean():.4e}")
print(f"  dimensionless |curl|*L/|v|      = {np.abs(w_intent)[2:-2,2:-2].mean()*L/typ:.4e}")
print(f"Reference solid-body rotation of same |v|:")
print(f"  mean|curl v_rot|                = {np.abs(w_rot)[2:-2,2:-2].mean():.4e}")
print(f"  dimensionless |curl|*L/|v|      = {np.abs(w_rot)[2:-2,2:-2].mean()*L/typ:.4e}")
ratio = np.abs(w_intent)[2:-2,2:-2].mean() / (np.abs(w_rot)[2:-2,2:-2].mean() + 1e-30)
print(f"\nVorticity ratio  intent/rotational = {ratio:.3e}")
print("=> induced curl is at discretization-noise level; analytically curl=0.")

# Also confirm divergence is NOT zero (so it is NOT incompressible N-S)
dv = divergence(vx, vy, dx)
print(f"\nIncompressibility check: mean|div v|*L/|v| = "
      f"{np.abs(dv)[2:-2,2:-2].mean()*L/typ:.4e}  (>>0 => compressible, NOT div-free)")

# ---------------------------------------------------------------------------
# (B) Evolve the actual scalar transfer rule from a vortex-ring IC
# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("(B) Evolve dI/dt = div[D R(I) grad(I)] from a ring; track vorticity & ring")
print("=" * 70)

def laplacian_flux_step(I, dx, dt):
    """One explicit step of dI/dt = div[ D R(I) grad(I) ] with periodic BCs."""
    gx, gy = grad(I, dx)
    Rf = R(I)
    Jx = D * Rf * gx
    Jy = D * Rf * gy
    # divergence of J via periodic-aware gradient
    dJx = np.gradient(Jx, dx, axis=1)
    dJy = np.gradient(Jy, dx, axis=0)
    return I + dt * (dJx + dJy)

n = 128
L = 1.0
dx = L / n
xs = np.linspace(0, L, n, endpoint=False)
X, Y = np.meshgrid(xs, xs)
r = np.sqrt((X - L/2) ** 2 + (Y - L/2) ** 2)

# Ring-shaped (annular) high-Intent structure -- the natural "vortex" candidate
r0, wr = 0.22, 0.04
I = 0.1 + 0.8 * np.exp(-((r - r0) ** 2) / (2 * wr ** 2))
I = np.clip(I, 1e-3, 0.98 * I_MAX)

dt = 0.2 * dx ** 2 / D     # parabolic stability
I0_peak = I.max()
ring_mass0 = I[(r > r0 - 2*wr) & (r < r0 + 2*wr)].sum()

print(f"{'step':>6} {'peak I':>10} {'ring frac':>10} {'mean|curl|*L/|v|':>18}")
for step in range(0, 4001):
    if step % 800 == 0:
        vx, vy = intent_velocity(I, dx)
        sp = np.sqrt(vx**2 + vy**2).mean() + 1e-12
        w = np.abs(curl_z(vx, vy, dx))[2:-2,2:-2].mean()
        ring_now = I[(r > r0 - 2*wr) & (r < r0 + 2*wr)].sum()
        print(f"{step:>6} {I.max():>10.4f} {ring_now/ring_mass0:>10.4f} {w*L/sp:>18.3e}")
    I = laplacian_flux_step(I, dx, dt)

print("\nInterpretation:")
print(" - peak I decays and the ring fraction drops => structure DISPERSES.")
print(" - induced vorticity stays at noise level throughout => never rotational.")
print(" - matches grid Sessions 19-22 (R(I) defocusing) -- analytically forced,")
print("   not a resolution/initialization artifact to be fixed on a bigger grid.")
