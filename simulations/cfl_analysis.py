"""
CFL Analysis of Discrete Intent Dynamics
=========================================
Session 13: Actually computing the cellular automaton dynamics
instead of critiquing the vocabulary.

The transfer rule: ΔI(x→y) = k · (I_x - I_y) · R(I_y)
where R(I) = [1 - (I/I_max)^n]

Question: Under what conditions does the discrete dynamics
produce oscillating modes (CFL violation)?

Linear stability analysis (1D, uniform background I_0):
  Amplification factor G(q) = 1 - 4k·R(I_0)·(1-cos(q))
  CFL condition: k·R(I_0) ≤ 1/(4d) where d = dimension

If CFL is violated, short-wavelength modes oscillate and grow.
Nonlinear saturation (R → 0 at I_max) caps the growth.
"""

import numpy as np
import json
from pathlib import Path


def R(I, I_max, n):
    """Saturation resistance function."""
    return 1.0 - (I / I_max) ** n


def simulate_1d(N, k, I_max, n, I_init, steps):
    """
    Simulate 1D discrete Intent dynamics with periodic boundaries.

    Returns: array of shape (steps+1, N) with Intent values.
    """
    history = np.zeros((steps + 1, N))
    I = I_init.copy()
    history[0] = I.copy()

    for t in range(steps):
        I_new = I.copy()
        for x in range(N):
            xp = (x + 1) % N
            xm = (x - 1) % N

            # Net change for cell x:
            # Received from x+1: k * (I[xp] - I[x]) * R(I[x])
            # Sent to x+1:       k * (I[x] - I[xp]) * R(I[xp])
            # Net from x+1:      k * (I[xp] - I[x]) * (R(I[x]) + R(I[xp]))

            dI_right = k * (I[xp] - I[x]) * (R(I[x], I_max, n) + R(I[xp], I_max, n))
            dI_left  = k * (I[xm] - I[x]) * (R(I[x], I_max, n) + R(I[xm], I_max, n))

            I_new[x] = I[x] + dI_right + dI_left

            # Enforce bounds
            I_new[x] = np.clip(I_new[x], 0, I_max)

        I = I_new
        history[t + 1] = I.copy()

    return history


def linear_stability_1d(k, I_0, I_max, n, d=1):
    """
    Compute amplification factor for q=pi mode (worst case).

    CFL condition: |G(pi)| <= 1
    G(pi) = 1 - 8*d*k*R(I_0)  (1D: d=1, factor 8; but our symmetric
    formulation has factor 4k*(R_x + R_y) = 8k*R for uniform state)

    Wait, let me recompute. For uniform I_0:
    G(q) = 1 - 4k*R(I_0)*(1-cos(q))  [per dimension]
    Total in d dimensions: 1 - 4k*R(I_0)*sum(1-cos(q_i))

    Worst case q_i = pi for all i:
    G = 1 - 4k*R*2d = 1 - 8dk*R

    CFL: 8dk*R <= 2, i.e., k*R <= 1/(4d)
    """
    R_0 = R(I_0, I_max, n)
    G_worst = 1 - 8 * d * k * R_0
    cfl_threshold = 1.0 / (4 * d)
    kR = k * R_0

    return {
        'R_0': R_0,
        'kR': kR,
        'G_worst': G_worst,
        'cfl_threshold': cfl_threshold,
        'cfl_violated': kR > cfl_threshold,
        'oscillating': G_worst < -1,
        'growing': abs(G_worst) > 1,
    }


def scan_stability(k_values, I_0_values, I_max, n, d=1):
    """Scan parameter space for CFL violation regions."""
    results = []
    for k in k_values:
        for I_0 in I_0_values:
            info = linear_stability_1d(k, I_0, I_max, n, d)
            results.append({
                'k': k,
                'I_0': I_0,
                'I_0_normalized': I_0 / I_max,
                **info
            })
    return results


def find_oscillation_modes(k, I_max, n, d=1):
    """
    For given k, find range of I_0 where CFL is violated.

    CFL violated when k * R(I_0) > 1/(4d)
    i.e., k * [1 - (I_0/I_max)^n] > 1/(4d)
    i.e., (I_0/I_max)^n < 1 - 1/(4dk)
    i.e., I_0/I_max < [1 - 1/(4dk)]^(1/n)

    For this to have solutions: 1 - 1/(4dk) > 0, i.e., k > 1/(4d)
    """
    if k <= 1.0 / (4 * d):
        return {
            'has_cfl_violation': False,
            'msg': f'k={k:.4f} <= 1/(4d)={1/(4*d):.4f}: no CFL violation possible',
            'I_threshold': None,
        }

    rho_threshold = (1.0 - 1.0 / (4 * d * k)) ** (1.0 / n)
    I_threshold = rho_threshold * I_max

    return {
        'has_cfl_violation': True,
        'I_threshold': I_threshold,
        'rho_threshold': rho_threshold,
        'msg': (f'CFL violated for I_0 < {I_threshold:.4f} '
                f'(ρ < {rho_threshold:.4f}). '
                f'Fraction of I-space that oscillates: {rho_threshold:.1%}'),
    }


def run_oscillation_test(k, I_max, n, N=64, steps=200):
    """
    Test whether CFL violation actually produces oscillations
    in the nonlinear dynamics.
    """
    # Initialize with small perturbation around low-I background
    I_bg = 0.1 * I_max  # Low I → R ≈ 1 → possible CFL violation
    I_init = I_bg * np.ones(N) + 0.01 * I_max * np.sin(2 * np.pi * np.arange(N) / 2)  # q=pi mode

    history = simulate_1d(N, k, I_max, n, I_init, steps)

    # Check for oscillation: does cell 0's value oscillate?
    cell_0 = history[:, 0]
    # Compute sign changes in derivative
    dI = np.diff(cell_0)
    sign_changes = np.sum(np.diff(np.sign(dI)) != 0)

    # Compute amplitude growth
    amplitude = np.abs(cell_0 - np.mean(cell_0))
    early_amp = np.mean(amplitude[1:11]) if len(amplitude) > 10 else amplitude[1]
    late_amp = np.mean(amplitude[-10:]) if len(amplitude) > 10 else amplitude[-1]

    return {
        'oscillates': sign_changes > steps * 0.3,  # Many sign changes = oscillation
        'sign_changes': int(sign_changes),
        'early_amplitude': float(early_amp),
        'late_amplitude': float(late_amp),
        'amplitude_ratio': float(late_amp / early_amp) if early_amp > 1e-15 else float('inf'),
        'final_range': float(np.max(history[-1]) - np.min(history[-1])),
        'cell_0_trace': cell_0[:20].tolist(),  # First 20 steps
    }


def main():
    print("=" * 70)
    print("CFL ANALYSIS OF DISCRETE INTENT DYNAMICS")
    print("=" * 70)

    I_max = 1.0  # Normalized

    print("\n--- Linear Stability Analysis ---")
    print(f"Transfer rule: ΔI(x→y) = k·(I_x - I_y)·R(I_y)")
    print(f"R(I) = [1 - (I/I_max)^n]")
    print(f"CFL condition (1D): k·R(I_0) ≤ 1/4")
    print(f"CFL condition (3D): k·R(I_0) ≤ 1/12")

    # Test different n values
    for n in [1, 2, 4]:
        print(f"\n{'='*50}")
        print(f"n = {n}")
        print(f"{'='*50}")

        for d in [1, 3]:
            print(f"\n  d = {d} dimensions:")

            # Find oscillation threshold
            for k in [0.1, 0.25, 0.3, 0.5, 1.0]:
                modes = find_oscillation_modes(k, I_max, n, d)
                marker = "✓ OSCILLATES" if modes['has_cfl_violation'] else "  stable"
                print(f"    k={k:.2f}: {marker} — {modes['msg']}")

    # Nonlinear simulation
    print("\n" + "=" * 70)
    print("NONLINEAR SIMULATION TESTS")
    print("=" * 70)

    n = 2  # Quadratic saturation

    for k in [0.1, 0.26, 0.5, 1.0]:
        print(f"\n--- k = {k} (CFL threshold at 0.25 for 1D) ---")

        # Linear prediction
        lin = linear_stability_1d(k, 0.1 * I_max, I_max, n, d=1)
        print(f"  Linear: G(π) = {lin['G_worst']:.4f}, "
              f"CFL violated: {lin['cfl_violated']}, "
              f"oscillating: {lin['oscillating']}")

        # Nonlinear simulation
        result = run_oscillation_test(k, I_max, n, N=64, steps=200)
        print(f"  Nonlinear: oscillates={result['oscillates']}, "
              f"sign_changes={result['sign_changes']}, "
              f"amp_ratio={result['amplitude_ratio']:.4f}")
        print(f"  Cell 0 first 10 steps: {[f'{v:.4f}' for v in result['cell_0_trace'][:10]]}")

    # KEY QUESTION: What happens at the CFL boundary with saturation?
    print("\n" + "=" * 70)
    print("KEY TEST: CFL + SATURATION = STABLE OSCILLATION?")
    print("=" * 70)
    print("If k > 1/4 AND saturation caps growth → standing wave?")

    k = 0.5
    n = 2
    N = 32

    for I_bg_frac in [0.05, 0.1, 0.3, 0.5, 0.8]:
        I_bg = I_bg_frac * I_max
        I_init = I_bg * np.ones(N)
        # Add q=pi perturbation (alternating high/low)
        I_init += 0.02 * I_max * np.array([(-1)**i for i in range(N)])
        I_init = np.clip(I_init, 0, I_max)

        history = simulate_1d(N, k, I_max, n, I_init, 500)

        cell_0 = history[:, 0]
        dI = np.diff(cell_0)
        sign_changes = np.sum(np.diff(np.sign(dI)) != 0)

        # Check if it reached steady oscillation (last 100 steps)
        late = cell_0[-100:]
        late_range = np.max(late) - np.min(late)
        late_sign_changes = np.sum(np.diff(np.sign(np.diff(late))) != 0)

        R_bg = R(I_bg, I_max, n)
        G = 1 - 8 * k * R_bg

        status = "SUSTAINED" if late_sign_changes > 60 and late_range > 0.001 else \
                 "DAMPED" if late_range < 0.001 else "GROWING/SATURATED"

        print(f"  I_bg/I_max={I_bg_frac:.2f}, R={R_bg:.4f}, G(π)={G:.4f}: "
              f"{status}, late_range={late_range:.6f}, "
              f"late_osc={late_sign_changes}")

    # Entity test: Can a localized high-I region oscillate?
    print("\n" + "=" * 70)
    print("ENTITY TEST: LOCALIZED HIGH-I PATTERN")
    print("=" * 70)

    k = 0.5
    n = 2
    N = 64

    # Create a "particle": gaussian peak of Intent
    x = np.arange(N)
    for peak_height in [0.3, 0.5, 0.7, 0.9]:
        I_init = 0.05 * I_max * np.ones(N)
        I_init += peak_height * I_max * np.exp(-((x - N//2) ** 2) / (2 * 3**2))
        I_init = np.clip(I_init, 0, I_max)

        history = simulate_1d(N, k, I_max, n, I_init, 1000)

        # Track the peak cell
        peak_cell = history[:, N//2]
        dI = np.diff(peak_cell)
        sign_changes = np.sum(np.diff(np.sign(dI)) != 0)

        # Check late behavior
        late = peak_cell[-200:]
        late_range = np.max(late) - np.min(late)
        late_sign = np.sum(np.diff(np.sign(np.diff(late))) != 0)

        # Does the peak survive?
        final_peak = np.max(history[-1]) - np.min(history[-1])

        print(f"  Peak height={peak_height:.1f}: "
              f"final_contrast={final_peak:.6f}, "
              f"peak_osc={late_sign}, peak_range={late_range:.6f}")
        if final_peak < 0.001:
            print(f"    → DIFFUSED AWAY (no entity)")
        elif late_sign > 100:
            print(f"    → OSCILLATING ENTITY")
        else:
            print(f"    → STATIC REMNANT (not oscillating)")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
Key results:
1. CFL threshold (1D): k·R(I_0) > 1/4 → short-wavelength oscillation
2. CFL threshold (3D): k·R(I_0) > 1/12 → more restrictive
3. R(I) decreases with I → CFL violated at LOW I, satisfied at HIGH I
4. This means BACKGROUND oscillates, entities (high I) are STABLE
5. Localized peaks: do they oscillate or diffuse?
   → This determines whether the cellular automaton produces entities

The CFL analysis reveals an INVERSION:
- The framework needs high-I patterns to oscillate (entities at f = E/h)
- The mathematics shows high-I regions are CFL-STABLE (R→0, no oscillation)
- Low-I regions (background) are CFL-UNSTABLE (R≈1, oscillation possible)
- Entities should be BACKGROUNDS, not peaks — or the model needs revision
""")


if __name__ == '__main__':
    main()
