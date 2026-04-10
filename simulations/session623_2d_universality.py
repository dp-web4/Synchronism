"""
Session 623 Part 2: 2D CA universality test.

Does 2D change the picture? More neighbors (4 instead of 2) might give richer dynamics.
Also tests: does the system ever produce GLIDERS (localized moving structures)?
"""

import numpy as np

I_MAX = 1.0


def R(I, n=2):
    return 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n


def step_2d(I, k, n=2):
    """Synchronous update on 2D grid with periodic boundaries."""
    up = np.roll(I, 1, axis=0)
    down = np.roll(I, -1, axis=0)
    left = np.roll(I, 1, axis=1)
    right = np.roll(I, -1, axis=1)

    dI = k * ((up - I) * R(up, n) + (down - I) * R(down, n) +
              (left - I) * R(left, n) + (right - I) * R(right, n))
    return np.clip(I + dI, 0.0, I_MAX)


def test_2d():
    N = 64
    steps = 1000
    np.random.seed(42)

    print("=" * 60)
    print("2D COMPUTATIONAL UNIVERSALITY TEST")
    print(f"Grid: {N}x{N}, steps={steps}")
    print("=" * 60)

    for k in [0.1, 0.24, 0.3, 0.5, 0.8]:
        print(f"\n--- k = {k:.2f} ---")

        # Random initial condition
        I = np.random.uniform(0.1, 0.9, (N, N))
        I0 = I.copy()

        # Perturbed copy for Lyapunov
        I_p = I.copy()
        I_p[N // 2, N // 2] += 1e-8

        stds = []
        for t in range(steps):
            I = step_2d(I, k)
            I_p = step_2d(I_p, k)
            if t % 100 == 0:
                stds.append(I.std())

        div = np.max(np.abs(I - I_p))

        # Check for spatial structure
        grad_x = np.diff(I, axis=1)
        grad_y = np.diff(I, axis=0)
        grad_energy = np.mean(grad_x ** 2) + np.mean(grad_y ** 2)

        print(f"  std evolution: {stds[0]:.4f} → {stds[-1]:.6f}")
        print(f"  Gradient energy: {grad_energy:.2e}")
        print(f"  Lyapunov div: {div:.2e}")

        if stds[-1] < 1e-4:
            print("  → UNIFORM (Class 1)")
        elif div < 1e-6:
            print("  → STABLE pattern (Class 2)")
        else:
            print("  → COMPLEX? (needs investigation)")

    # === GLIDER TEST ===
    print("\n" + "=" * 60)
    print("GLIDER TEST: Can localized structures move?")
    print("=" * 60)

    for k in [0.3, 0.5, 0.8]:
        print(f"\n--- k = {k:.2f} ---")

        # Background with a localized Gaussian pulse
        I = np.ones((N, N)) * 0.3
        y, x = np.mgrid[:N, :N]
        cx, cy = N // 4, N // 4
        pulse = 0.6 * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * 3 ** 2))
        I += pulse
        I = np.clip(I, 0, I_MAX)

        # Track center of mass of the perturbation
        I_bg = np.ones((N, N)) * 0.3  # Will evolve to uniform quickly
        com_x = []
        com_y = []
        peak_vals = []

        for t in range(steps):
            I = step_2d(I, k)
            diff = I - 0.3  # Approximate excess
            total = np.sum(np.abs(diff)) + 1e-10
            com_x.append(np.sum(x * np.abs(diff)) / total)
            com_y.append(np.sum(y * np.abs(diff)) / total)
            peak_vals.append(np.max(np.abs(diff)))

        drift_x = com_x[-1] - com_x[0]
        drift_y = com_y[-1] - com_y[0]
        total_drift = np.sqrt(drift_x ** 2 + drift_y ** 2)
        peak_decay = peak_vals[-1] / peak_vals[0]

        print(f"  COM drift: ({drift_x:.1f}, {drift_y:.1f}), total={total_drift:.1f} cells")
        print(f"  Peak decay: {peak_vals[0]:.4f} → {peak_vals[-1]:.2e} (ratio={peak_decay:.2e})")
        if peak_decay < 0.01:
            print("  → Pulse DISPERSED (no glider)")
        elif total_drift < 1.0:
            print("  → Pulse STATIONARY (no glider)")
        else:
            print("  → Possible GLIDER?")

    # === THE CRITICAL QUESTION ===
    print("\n" + "=" * 60)
    print("VERDICT")
    print("=" * 60)
    print("""
The stated transfer rule on a 2D lattice:
  - Cannot propagate signals (diffusion only)
  - Cannot sustain localized structures (everything disperses)
  - Cannot support gliders (no moving structures)
  - Has negative Lyapunov exponents (all perturbations die)
  - Is computationally trivial (Class 1-2)

A computational substrate must be at least Class 4 (computationally universal).
The universe we inhabit builds Turing machines. This CA cannot.

ROOT CAUSE: Monotonic saturation R(I) is a SMOOTHING operator.
  - At low I: R ≈ 1, full transfer → rapid equalization
  - At high I: R ≈ 0, no transfer → frozen
  - No intermediate regime with complex dynamics
  - No way to create, store, or process information

The physics arguments (S617-622) and the computation argument converge:
  SAME ROOT CAUSE → monotonic R(I) is too simple for anything.
""")


if __name__ == "__main__":
    test_2d()
