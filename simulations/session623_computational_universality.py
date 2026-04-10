"""
Session 623: Computational Universality Test

The stated transfer rule defines a cellular automaton:
  I_new(x) = I_old(x) + k * sum_n [I_old(n) - I_old(x)] * R(I_old(n))
  R(I) = [1 - (I/I_max)^n]

Question: What Wolfram class is this CA? Can it support complex computation?

A universe substrate must be at least Class 4 (computationally universal).
If the CA is Class 1-2, it's provably too simple to be the universe.

Tests:
  1. Entropy evolution from random initial conditions (does complexity grow, decay, or stabilize?)
  2. Mutual information between distant cells (does information propagate non-trivially?)
  3. Lyapunov-like sensitivity (does a 1-bit perturbation grow, shrink, or stay bounded?)
  4. Pattern diversity (does the system generate diverse long-lived structures?)
"""

import numpy as np
from collections import Counter

I_MAX = 1.0


def R(I, n=2):
    """Saturation resistance function."""
    return 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n


def step_1d(I, k, n=2):
    """One synchronous update step in 1D with periodic boundaries."""
    I_left = np.roll(I, 1)
    I_right = np.roll(I, -1)
    dI = k * ((I_left - I) * R(I_left, n) + (I_right - I) * R(I_right, n))
    I_new = np.clip(I + dI, 0.0, I_MAX)
    return I_new


def spatial_entropy(I, n_bins=20):
    """Shannon entropy of the spatial distribution."""
    hist, _ = np.histogram(I, bins=n_bins, range=(0, I_MAX))
    p = hist / hist.sum()
    p = p[p > 0]
    return -np.sum(p * np.log2(p))


def block_entropy(I, block_size=3, n_levels=10):
    """Block entropy: entropy of consecutive block patterns."""
    # Discretize
    levels = np.clip((I / I_MAX * n_levels).astype(int), 0, n_levels - 1)
    blocks = []
    for i in range(len(levels) - block_size + 1):
        blocks.append(tuple(levels[i:i + block_size]))
    counts = Counter(blocks)
    total = sum(counts.values())
    p = np.array([c / total for c in counts.values()])
    return -np.sum(p * np.log2(p))


def mutual_information_distance(I, separation, n_bins=10):
    """Mutual information between cells separated by `separation`."""
    x = I[:-separation] if separation > 0 else I
    y = np.roll(I, -separation)[:-separation] if separation > 0 else I

    # Joint histogram
    hist_xy, _, _ = np.histogram2d(x, y, bins=n_bins, range=[[0, I_MAX], [0, I_MAX]])
    pxy = hist_xy / hist_xy.sum()
    px = pxy.sum(axis=1)
    py = pxy.sum(axis=0)

    # MI = sum p(x,y) log(p(x,y) / p(x)p(y))
    mi = 0.0
    for i in range(n_bins):
        for j in range(n_bins):
            if pxy[i, j] > 0 and px[i] > 0 and py[j] > 0:
                mi += pxy[i, j] * np.log2(pxy[i, j] / (px[i] * py[j]))
    return mi


def lyapunov_test(I0, k, n, steps, perturbation=1e-8):
    """Track divergence of two initially close states."""
    I_a = I0.copy()
    I_b = I0.copy()
    # Perturb one cell
    mid = len(I0) // 2
    I_b[mid] += perturbation

    divergences = []
    for t in range(steps):
        I_a = step_1d(I_a, k, n)
        I_b = step_1d(I_b, k, n)
        div = np.max(np.abs(I_a - I_b))
        divergences.append(div)
    return divergences


def count_distinct_patterns(I, threshold=0.01, window=5):
    """Count distinct local patterns above a threshold."""
    levels = (I / I_MAX * 10).astype(int)
    patterns = set()
    for i in range(len(levels) - window + 1):
        patterns.add(tuple(levels[i:i + window]))
    return len(patterns)


def run_test(grid_size=256, k_values=None, n=2, steps=2000):
    """Run full computational universality test suite."""
    if k_values is None:
        k_values = [0.1, 0.3, 0.49, 0.55, 0.8, 1.0]

    np.random.seed(42)

    print("=" * 70)
    print("COMPUTATIONAL UNIVERSALITY TEST")
    print(f"Grid: {grid_size}, n={n}, steps={steps}")
    print("=" * 70)

    for k in k_values:
        print(f"\n--- k = {k:.2f} ---")

        # Random initial condition
        I0 = np.random.uniform(0.1, 0.9, grid_size)

        # Evolve and track metrics
        I = I0.copy()
        entropies = []
        block_entropies = []
        pattern_counts = []
        mi_near = []
        mi_far = []

        for t in range(steps):
            I = step_1d(I, k, n)
            if t % 100 == 0:
                entropies.append(spatial_entropy(I))
                block_entropies.append(block_entropy(I))
                pattern_counts.append(count_distinct_patterns(I))
                mi_near.append(mutual_information_distance(I, 1))
                mi_far.append(mutual_information_distance(I, grid_size // 4))

        # Final state analysis
        I_range = I.max() - I.min()
        I_std = I.std()

        print(f"  Final state: mean={I.mean():.4f}, std={I_std:.6f}, range={I_range:.6f}")
        print(f"  Entropy: initial={entropies[0]:.3f}, final={entropies[-1]:.3f}")
        print(f"  Block entropy: initial={block_entropies[0]:.3f}, final={block_entropies[-1]:.3f}")
        print(f"  Pattern diversity: initial={pattern_counts[0]}, final={pattern_counts[-1]}")
        print(f"  MI(1): initial={mi_near[0]:.4f}, final={mi_near[-1]:.4f}")
        print(f"  MI(N/4): initial={mi_far[0]:.4f}, final={mi_far[-1]:.4f}")

        # Wolfram classification heuristic
        if I_std < 1e-4:
            wolfram_class = "Class 1 (uniform)"
        elif I_range < 0.01 and entropies[-1] < 0.5:
            wolfram_class = "Class 1 (uniform)"
        elif pattern_counts[-1] < 10:
            wolfram_class = "Class 2 (periodic/simple)"
        elif entropies[-1] > entropies[0] * 0.8 and pattern_counts[-1] > 50:
            wolfram_class = "Class 3 (chaotic) or Class 4 (complex)"
        else:
            wolfram_class = "Class 2 (periodic/simple)"
        print(f"  Classification: {wolfram_class}")

        # Lyapunov test
        divs = lyapunov_test(I0.copy(), k, n, steps=500)
        max_div = max(divs)
        final_div = divs[-1]
        print(f"  Lyapunov: max_div={max_div:.2e}, final_div={final_div:.2e}", end="")
        if max_div < 1e-6:
            print(" → STABLE (perturbation dies)")
        elif final_div > 0.1:
            print(" → UNSTABLE (perturbation grows)")
        else:
            print(" → BOUNDED (perturbation persists)")

    # === CRITICAL TEST: Can the CA transmit a signal? ===
    print("\n" + "=" * 70)
    print("SIGNAL PROPAGATION TEST")
    print("=" * 70)

    for k in [0.3, 0.55, 1.0]:
        print(f"\n--- k = {k:.2f} ---")
        # Uniform background with a localized pulse
        I = np.ones(grid_size) * 0.3
        I[grid_size // 2] = 0.9

        # Track where the pulse "information" is
        I_bg = np.ones(grid_size) * 0.3  # Reference: no pulse
        for t in range(steps):
            I = step_1d(I, k, n)
            I_bg = step_1d(I_bg, k, n)

        diff = np.abs(I - I_bg)
        max_diff = diff.max()
        spread = np.sum(diff > 1e-6)

        print(f"  Max difference from no-pulse: {max_diff:.2e}")
        print(f"  Cells affected (>1e-6): {spread}/{grid_size}")

        if max_diff < 1e-6:
            print("  → Signal VANISHED (diffusion erased it)")
        elif spread > grid_size * 0.9:
            print("  → Signal SPREAD uniformly (diffusion)")
        else:
            print("  → Signal LOCALIZED (propagation)")

    # === GATE TEST: Can two signals interact? ===
    print("\n" + "=" * 70)
    print("SIGNAL INTERACTION (GATE) TEST")
    print("=" * 70)
    print("For computational universality, two signals must interact non-trivially.")
    print("Testing: does pulse A + pulse B produce something ≠ pulse A + pulse B separately?")

    for k in [0.3, 0.55, 1.0]:
        print(f"\n--- k = {k:.2f} ---")

        I_bg = np.ones(grid_size) * 0.3

        # Pulse A only
        I_a = I_bg.copy()
        I_a[grid_size // 3] = 0.9

        # Pulse B only
        I_b = I_bg.copy()
        I_b[2 * grid_size // 3] = 0.9

        # Both pulses
        I_ab = I_bg.copy()
        I_ab[grid_size // 3] = 0.9
        I_ab[2 * grid_size // 3] = 0.9

        for t in range(steps):
            I_a = step_1d(I_a, k, n)
            I_b = step_1d(I_b, k, n)
            I_ab = step_1d(I_ab, k, n)
            I_bg_t = step_1d(I_bg, k, n) if t == 0 else I_bg
            # (I_bg is uniform — stays uniform under the rule)

        # Superposition test: does I_ab = I_a + I_b - I_bg?
        I_superposition = I_a + I_b - 0.3  # Linear superposition prediction
        nonlinearity = np.max(np.abs(I_ab - I_superposition))

        print(f"  Max nonlinear interaction: {nonlinearity:.2e}")
        if nonlinearity < 1e-6:
            print("  → LINEAR: signals don't interact (no gates possible)")
        elif nonlinearity < 0.01:
            print("  → WEAK nonlinearity: minimal interaction")
        else:
            print(f"  → NONLINEAR interaction detected (strength: {nonlinearity:.4f})")
            # But does it produce structure or just mixing?
            diff_ab = I_ab - 0.3
            structure = np.std(diff_ab) / (np.mean(np.abs(diff_ab)) + 1e-10)
            print(f"     Structure ratio (std/mean|diff|): {structure:.3f}")


if __name__ == "__main__":
    run_test()
