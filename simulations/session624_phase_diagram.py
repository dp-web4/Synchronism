"""
Session 624 Part 2: Phase Diagram and 2-DOF Test

Finding from Part 1: Non-monotonic R produces Class 3 (chaotic) CA —
the first non-trivial dynamics in 624 sessions. But self-confinement
and gravity+waves still fail.

This part:
1. Map the full (A, k) phase diagram — where exactly is chaos? Is there Class 4?
2. Test 2-DOF (I + velocity) with non-monotonic R — does combining the
   two fixes (momentum + non-monotonicity) produce self-confinement?
3. Test whether the chaotic regime has STRUCTURE (gliders, memory) or is
   just random noise (which would be Class 3, not Class 4).
"""

import numpy as np
from collections import Counter

I_MAX = 1.0


def R_nonmonotonic(I, n=2, A=0.5):
    x = np.clip(I, 0, I_MAX) / I_MAX
    base = 1.0 - x**n
    bump = 1.0 + A * np.sin(np.pi * x)
    return np.maximum(0, base * bump)


def step_1d(I, k, A=0.5, n=2):
    I_left = np.roll(I, 1)
    I_right = np.roll(I, -1)
    R_left = R_nonmonotonic(I_left, n, A)
    R_right = R_nonmonotonic(I_right, n, A)
    dI = k * ((I_left - I) * R_left + (I_right - I) * R_right)
    return np.clip(I + dI, 0.0, I_MAX)


def spatial_entropy(I, n_bins=20):
    hist, _ = np.histogram(I, bins=n_bins, range=(0, I_MAX))
    p = hist / hist.sum()
    p = p[p > 0]
    return -np.sum(p * np.log2(p))


def lyapunov_1d(I0, k, A, steps=300, eps=1e-8):
    I = I0.copy()
    I2 = I0.copy()
    I2[len(I)//2] += eps
    for _ in range(steps):
        I = step_1d(I, k, A)
        I2 = step_1d(I2, k, A)
    diff = np.max(np.abs(I - I2))
    return np.log10(diff / eps) if diff > 0 else -20


# ============================================================
# Phase Diagram: A vs k
# ============================================================
print("=" * 70)
print("PHASE DIAGRAM: Non-monotonic amplitude A vs coupling k")
print("=" * 70)

N = 256
np.random.seed(42)
I0 = np.random.uniform(0.1, 0.9, N)

A_vals = np.arange(0.0, 2.51, 0.1)
k_vals = np.arange(0.1, 1.01, 0.05)

# Store results
phase_map = {}
class_counts = Counter()

print(f"\n{'A':>5s} {'k':>5s} {'S_final':>8s} {'S_std':>8s} {'Lyap':>8s} {'Patterns':>8s} {'Class':>10s}")
print("-" * 60)

for A in A_vals:
    for k in k_vals:
        I = I0.copy()
        entropies = []
        for t in range(500):
            entropies.append(spatial_entropy(I))
            I = step_1d(I, k, A)

        S_final = np.mean(entropies[-50:])
        S_std = np.std(entropies[-50:])

        # Discretize and count patterns
        levels = np.clip((I / I_MAX * 10).astype(int), 0, 9)
        blocks = set()
        for i in range(len(levels) - 2):
            blocks.add(tuple(levels[i:i+3]))
        n_patterns = len(blocks)

        lyap = lyapunov_1d(I0.copy(), k, A)

        # Classify
        if S_final < 0.5 and S_std < 0.01:
            wclass = 1
        elif lyap > 3:
            wclass = 3
        elif S_std > 0.05 and n_patterns > 30 and lyap < 3 and lyap > -2:
            wclass = 4  # Edge of chaos: some complexity, bounded divergence
        elif n_patterns < 30:
            wclass = 2
        else:
            wclass = 2  # Default to periodic if unclear

        phase_map[(A, k)] = wclass
        class_counts[wclass] += 1

# Print only transitions (where class changes)
for A in [0.0, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0]:
    prev_class = None
    for k in k_vals:
        A_key = min(A_vals, key=lambda x: abs(x - A))
        k_key = min(k_vals, key=lambda x: abs(x - k))
        c = phase_map.get((A_key, k_key), 0)
        if c != prev_class:
            I = I0.copy()
            for t in range(500):
                I = step_1d(I, k_key, A_key)
            S = spatial_entropy(I)
            levels = np.clip((I / I_MAX * 10).astype(int), 0, 9)
            blocks = set()
            for i in range(len(levels) - 2):
                blocks.add(tuple(levels[i:i+3]))
            n_pat = len(blocks)
            lyap = lyapunov_1d(I0.copy(), k_key, A_key)
            print(f"  A={A_key:.1f} k={k_key:.2f}: Class {c} (S={S:.2f}, pat={n_pat}, lyap={lyap:+.1f})")
            prev_class = c

print(f"\nClass distribution: {dict(class_counts)}")

# ============================================================
# ASCII Phase Map
# ============================================================
print("\n" + "=" * 70)
print("ASCII PHASE MAP (A rows, k columns)")
print("=" * 70)
symbols = {1: '.', 2: 'o', 3: '#', 4: '*'}
print(f"  Legend: .=Class1 o=Class2 #=Class3 *=Class4")
print(f"  {'A\\k':>5s}", end="")
for k in k_vals[::2]:
    print(f" {k:.2f}", end="")
print()

for A in A_vals:
    print(f"  {A:5.1f} ", end="")
    for k in k_vals[::2]:
        c = phase_map.get((A, k), 0)
        print(f"  {symbols.get(c, '?')} ", end="")
    print()


# ============================================================
# 2-DOF with Non-monotonic R: Does momentum + chaos = confinement?
# ============================================================
print("\n" + "=" * 70)
print("2-DOF TEST: Momentum + Non-monotonic R")
print("=" * 70)
print("Adding velocity field (momentum conservation) to the non-monotonic CA.")
print("This combines S17's fix (2-DOF) with the new fix (non-monotonic R).")

def step_2dof(I, v, k, A=1.0, n=2, damping=0.0):
    """2-DOF update: I (density) + v (velocity).

    Euler equations with non-monotonic viscosity:
    ∂I/∂t + ∂(Iv)/∂x = 0  (continuity)
    ∂v/∂t + v·∂v/∂x = -(1/I)·∂P/∂x + viscous  (momentum)

    P = integral of R_nm(I) — so dP/dI = R_nm(I) which is non-monotonic.
    """
    dx = 1.0

    # Pressure from non-monotonic R: P = ∫R dI
    # dP/dI = R_nm(I)
    R_vals = R_nonmonotonic(I, n, A)

    # Pressure gradient
    P_right = np.cumsum(R_nonmonotonic(np.roll(I, -1), n, A)) * dx
    P_left = np.cumsum(R_nonmonotonic(np.roll(I, 1), n, A)) * dx

    # Better: compute P directly and take gradient
    # P(I) = integral from 0 to I of R(I') dI'
    # Approximate: P ≈ R(I) * I (crude but captures the sign structure)
    # Actually, for the Euler equations, we just need dP/dx:

    I_safe = np.maximum(I, 1e-10)

    # Sound speed: c² = dP/dρ = R_nm(ρ)
    # For non-monotonic R, c² can be large at intermediate ρ
    cs2 = R_vals

    # Pressure gradient (centered finite difference)
    I_left = np.roll(I, 1)
    I_right = np.roll(I, -1)
    dPdx = (R_nonmonotonic(I_right, n, A) * I_right - R_nonmonotonic(I_left, n, A) * I_left) / (2 * dx)

    # Advection (upwind)
    dIdx = np.where(v > 0, I - I_left, I_right - I) / dx
    dvdx = np.where(v > 0, v - np.roll(v, 1), np.roll(v, -1) - v) / dx

    # Continuity: dI/dt = -d(Iv)/dx
    v_right = np.roll(v, -1)
    v_left = np.roll(v, 1)
    flux = (I_right * v_right - I_left * v_left) / (2 * dx)

    # Momentum: dv/dt = -v·dv/dx - (1/I)·dP/dx - damping·v
    dv = -v * dvdx - dPdx / I_safe - damping * v

    # Simple Euler step with CFL-limited dt
    cs_max = np.sqrt(np.maximum(np.max(np.abs(cs2)), 1e-10))
    v_max = np.max(np.abs(v)) + cs_max
    dt = min(0.3 * dx / max(v_max, 1e-10), 0.1)

    I_new = I - dt * flux
    v_new = v + dt * dv

    I_new = np.clip(I_new, 1e-8, I_MAX)
    v_new = np.clip(v_new, -10, 10)

    return I_new, v_new, dt


# Test: Gaussian pulse in 2-DOF non-monotonic system
N = 256
x = np.arange(N, dtype=float)

for A in [0.5, 1.0, 1.5]:
    for bg in [0.1, 0.3]:
        I = np.ones(N) * bg
        I += 0.6 * np.exp(-((x - N//2)**2) / (2 * 8**2))
        I = np.clip(I, 1e-8, I_MAX)
        v = np.zeros(N)

        init_width = np.sum(I > bg + 0.1)

        stable = True
        widths = []
        max_vals = []

        for t in range(2000):
            try:
                I, v, dt = step_2dof(I, v, k=0.5, A=A)
            except:
                stable = False
                break

            if np.any(np.isnan(I)) or np.any(np.isnan(v)):
                stable = False
                break

            w = np.sum(I > bg + 0.05)
            widths.append(w)
            max_vals.append(np.max(I))

        if stable and len(widths) > 100:
            final_width = widths[-1]
            final_max = max_vals[-1]

            # Check for oscillation in width (sign of self-organization)
            w_arr = np.array(widths[100:])
            w_std = np.std(w_arr)
            w_mean = np.mean(w_arr)

            if final_width <= init_width * 1.5 and final_width > 0:
                status = f"CONFINED (w: {init_width}→{final_width}, max: {max_vals[0]:.3f}→{final_max:.3f})"
            elif w_std > 5 and w_mean < N * 0.8:
                status = f"OSCILLATING (w: {init_width}→{final_width}±{w_std:.0f}, max: {final_max:.3f})"
            else:
                status = f"Dispersed (w: {init_width}→{final_width}, max: {final_max:.3f})"
        elif not stable:
            status = f"UNSTABLE at t={t}"
        else:
            status = "Incomplete"

        print(f"  A={A:.1f} bg={bg:.1f}: {status}")


# ============================================================
# Test 7: The Real Question — Does chaos + momentum = structure?
# ============================================================
print("\n" + "=" * 70)
print("TEST: RANDOM IC IN 2-DOF NON-MONOTONIC — SPONTANEOUS STRUCTURE?")
print("=" * 70)
print("Does starting from random IC produce spontaneous localized structures?")

np.random.seed(42)

for A in [0.8, 1.0, 1.5]:
    I = np.random.uniform(0.2, 0.8, N)
    v = np.random.uniform(-0.1, 0.1, N)

    stable = True
    for t in range(3000):
        try:
            I, v, dt = step_2dof(I, v, k=0.5, A=A, damping=0.01)
        except:
            stable = False
            break
        if np.any(np.isnan(I)):
            stable = False
            break

    if stable:
        # Look for peaks (localized structures)
        mean_I = np.mean(I)
        std_I = np.std(I)
        peaks = np.sum(I > mean_I + 3 * std_I)

        # Pattern diversity
        levels = np.clip((I / I_MAX * 20).astype(int), 0, 19)
        blocks = set()
        for i in range(len(levels) - 2):
            blocks.add(tuple(levels[i:i+3]))

        S = spatial_entropy(I, 30)
        print(f"  A={A:.1f}: S={S:.3f} std={std_I:.4f} peaks={peaks} "
              f"patterns={len(blocks)} range=[{np.min(I):.3f}, {np.max(I):.3f}]")
    else:
        print(f"  A={A:.1f}: Unstable at t={t}")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
