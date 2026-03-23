"""
Elastic Boundary Test — Conservation Bug Hypothesis
=====================================================
Tests whether adding momentum reflection at saturation boundaries
produces oscillating entities in the discrete Intent dynamics.

The original rule: ΔI = k·Σ(I_n - I)·R(I_n)
Problem: when R → 0, momentum vanishes (implicit sink)
Fix: track velocity field v, reflect at saturation boundaries

If this produces sustained oscillation, the 810 failed configurations
from Sessions #18-27 were testing a broken conservation law.
"""

import numpy as np


def R(I, I_max, n):
    """Saturation resistance function."""
    return 1.0 - (I / I_max) ** n


def simulate_elastic_1d(N, k, I_max, n, I_init, v_init, steps,
                        reflect_threshold=0.01):
    """
    1D Intent dynamics with momentum reflection at saturation boundaries.

    Two fields:
      I[x] = Intent density (scalar, conserved)
      v[x] = velocity (directed, reflected at saturation)

    Update:
      1. Compute advective flux: J = I * v
      2. Compute diffusive flux: J_diff = -D * R(I) * grad(I)
      3. Update I from total flux (conserving total I)
      4. When transfer blocked by R → 0, reflect velocity
    """
    I_history = np.zeros((steps + 1, N))
    v_history = np.zeros((steps + 1, N))
    energy_history = np.zeros(steps + 1)

    I = I_init.copy()
    v = v_init.copy()
    I_history[0] = I.copy()
    v_history[0] = v.copy()
    energy_history[0] = np.sum(I) + 0.5 * np.sum(I * v**2)

    dt = 1.0  # one tick

    for t in range(steps):
        I_new = I.copy()
        v_new = v.copy()

        for x in range(N):
            xp = (x + 1) % N
            xm = (x - 1) % N

            # --- Advective transport (momentum-carrying) ---
            # Flow from x in direction of v[x]
            if v[x] > 0:
                target = xp
                R_target = R(I[target], I_max, n)
                if R_target > reflect_threshold:
                    # Transfer allowed: move intent in flow direction
                    transfer = min(k * I[x] * abs(v[x]) * R_target, I[x] * 0.25)
                    I_new[x] -= transfer
                    I_new[target] += transfer
                else:
                    # ELASTIC REFLECTION: boundary blocks, reverse velocity
                    v_new[x] = -abs(v[x])
            elif v[x] < 0:
                target = xm
                R_target = R(I[target], I_max, n)
                if R_target > reflect_threshold:
                    transfer = min(k * I[x] * abs(v[x]) * R_target, I[x] * 0.25)
                    I_new[x] -= transfer
                    I_new[target] += transfer
                else:
                    # ELASTIC REFLECTION
                    v_new[x] = abs(v[x])

            # --- Diffusive transport (non-directional) ---
            for neighbor in [xp, xm]:
                R_n = R(I[neighbor], I_max, n)
                R_x = R(I[x], I_max, n)
                diff_flux = k * 0.1 * (I[neighbor] - I[x]) * (R_x + R_n) * 0.5
                I_new[x] += diff_flux

        # Enforce bounds
        I_new = np.clip(I_new, 0, I_max)

        # Velocity damping from viscosity (small, proportional to R)
        # High I → low R → low damping (entities are low-viscosity)
        for x in range(N):
            visc_damp = 1.0 - 0.01 * R(I_new[x], I_max, n)
            v_new[x] *= visc_damp

        I = I_new
        v = v_new
        I_history[t + 1] = I.copy()
        v_history[t + 1] = v.copy()
        energy_history[t + 1] = np.sum(I) + 0.5 * np.sum(I * v**2)

    return I_history, v_history, energy_history


def check_oscillation(trace, last_frac=0.5):
    """Check if a trace oscillates in its later portion."""
    n = len(trace)
    start = int(n * (1 - last_frac))
    late = trace[start:]
    if len(late) < 10:
        return False, 0, 0

    dI = np.diff(late)
    sign_changes = np.sum(np.diff(np.sign(dI)) != 0)
    amplitude = np.max(late) - np.min(late)

    # Oscillation: many sign changes AND non-trivial amplitude
    oscillates = sign_changes > len(late) * 0.2 and amplitude > 0.001
    return oscillates, sign_changes, amplitude


def main():
    print("=" * 70)
    print("ELASTIC BOUNDARY TEST — Conservation Bug Hypothesis")
    print("=" * 70)

    I_max = 1.0
    n = 2
    N = 32

    # ================================================================
    # TEST 1: Gaussian peak with initial velocity → does it bounce?
    # ================================================================
    print("\n--- TEST 1: Gaussian peak with rightward velocity ---")

    x = np.arange(N)
    I_init = 0.05 * np.ones(N)
    I_init += 0.5 * np.exp(-((x - N // 4) ** 2) / (2 * 2 ** 2))
    I_init = np.clip(I_init, 0, I_max)

    v_init = np.zeros(N)
    v_init[N // 4 - 3:N // 4 + 4] = 0.5  # rightward velocity at peak

    for k_val in [0.1, 0.3, 0.5]:
        I_hist, v_hist, E_hist = simulate_elastic_1d(
            N, k_val, I_max, n, I_init, v_init, 500)

        peak_trace = I_hist[:, N // 4]
        osc, sc, amp = check_oscillation(peak_trace)

        # Check where the peak is at different times
        peak_pos = [np.argmax(I_hist[t]) for t in [0, 50, 100, 200, 500]]

        E_conserved = abs(E_hist[-1] - E_hist[0]) / max(E_hist[0], 1e-10) < 0.1

        print(f"  k={k_val:.1f}: peak_positions={peak_pos}, "
              f"osc={osc}, sign_changes={sc}, amp={amp:.4f}, "
              f"E_conserved={E_conserved}")

    # ================================================================
    # TEST 2: Confined cavity — high-I walls with energy inside
    # ================================================================
    print("\n--- TEST 2: Confined cavity with reflecting walls ---")

    for cavity_size in [8, 12, 16]:
        I_init = np.ones(N) * 0.95  # high-I background (walls)
        center = N // 2
        half = cavity_size // 2
        I_init[center - half:center + half] = 0.3  # low-I cavity
        I_init[center] = 0.6  # energy inside cavity

        v_init = np.zeros(N)
        v_init[center] = 0.5  # initial velocity

        I_hist, v_hist, E_hist = simulate_elastic_1d(
            N, 0.3, I_max, n, I_init, v_init, 2000)

        center_trace = I_hist[:, center]
        osc, sc, amp = check_oscillation(center_trace)

        # Check multiple cells in cavity
        any_osc = False
        for cell in range(center - half, center + half):
            cell_osc, _, _ = check_oscillation(I_hist[:, cell])
            if cell_osc:
                any_osc = True
                break

        print(f"  cavity_size={cavity_size}: center_osc={osc}, any_cell_osc={any_osc}, "
              f"center_sign_changes={sc}, center_amp={amp:.6f}")
        print(f"    center trace [0,100,500,1000,2000]: "
              f"{[f'{I_hist[t, center]:.4f}' for t in [0, 100, 500, 1000, min(1999, 2000)]]}")

    # ================================================================
    # TEST 3: Systematic sweep — k, cavity size, initial velocity
    # ================================================================
    print("\n--- TEST 3: Systematic sweep ---")
    print(f"{'k':>6} {'cavity':>6} {'v0':>6} {'osc':>5} {'sign_ch':>8} {'amp':>10} {'status':>12}")

    oscillation_count = 0
    total_count = 0

    for k_val in [0.1, 0.2, 0.3, 0.5]:
        for cavity_size in [6, 10, 16]:
            for v0 in [0.2, 0.5, 1.0]:
                I_init = np.ones(N) * 0.95
                center = N // 2
                half = cavity_size // 2
                I_init[center - half:center + half] = 0.2
                I_init[center] = 0.5

                v_init = np.zeros(N)
                v_init[center] = v0

                I_hist, v_hist, E_hist = simulate_elastic_1d(
                    N, k_val, I_max, n, I_init, v_init, 1000)

                # Check ALL cavity cells for oscillation
                best_osc = False
                best_sc = 0
                best_amp = 0
                for cell in range(center - half, center + half):
                    o, s, a = check_oscillation(I_hist[:, cell])
                    if s > best_sc:
                        best_sc = s
                        best_amp = a
                        best_osc = o

                total_count += 1
                if best_osc:
                    oscillation_count += 1

                status = "OSCILLATES" if best_osc else "static"
                print(f"{k_val:>6.1f} {cavity_size:>6d} {v0:>6.1f} "
                      f"{str(best_osc):>5} {best_sc:>8d} {best_amp:>10.6f} "
                      f"{status:>12}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTotal configurations tested: {total_count}")
    print(f"Configurations with oscillation: {oscillation_count}")
    print(f"Oscillation rate: {oscillation_count}/{total_count} "
          f"({100*oscillation_count/max(total_count,1):.1f}%)")

    if oscillation_count > 0:
        print("\n*** ELASTIC BOUNDARIES PRODUCE OSCILLATION ***")
        print("The conservation bug hypothesis is SUPPORTED.")
        print("The 810 prior failures tested a broken conservation law.")
    else:
        print("\n*** NO OSCILLATION EVEN WITH ELASTIC BOUNDARIES ***")
        print("The conservation bug hypothesis is NOT SUPPORTED.")
        print("The problem is deeper than momentum conservation.")


if __name__ == '__main__':
    main()
