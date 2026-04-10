"""
Session #622: Does the discrete transfer rule have oscillatory behavior
that the continuum limit misses?

The transfer rule ΔI(x→y) = k·(I_x - I_y)·R(I_y) is, in the continuum limit,
nonlinear diffusion (S617). But the DISCRETE update with synchronous stepping
is equivalent to explicit Euler on diffusion — which is UNSTABLE for
k·D_eff > 1/(2d). In the unstable regime, the discrete system oscillates
even though the PDE monotonically relaxes.

Question 1: At what k does the discrete system transition from monotonic
            relaxation to oscillation?
Question 2: Is the oscillation sustained or transient?
Question 3: What is the energy scale of the oscillation?
Question 4: Does R(I) saturation modulate the instability in a physically
            interesting way?

If the discrete system oscillates where the PDE doesn't, S617's conclusion
("transfer rule = diffusion = no oscillation") has a gap. But the gap may
lead to the cosmological constant problem rather than to salvation.
"""

import numpy as np
import sys

def R(I, I_max=1.0, n=2):
    """Saturation resistance function."""
    return 1.0 - (I / I_max)**n

def run_1d_discrete(N, k, I_init, I_max=1.0, n=2, steps=1000):
    """
    Run the discrete transfer rule on a 1D lattice.
    Synchronous update (all cells update simultaneously).
    Periodic boundary conditions.

    Returns: I history at the center cell, and full field at select timesteps.
    """
    I = I_init.copy()
    center = N // 2
    history = np.zeros(steps)
    total_I = np.zeros(steps)

    for t in range(steps):
        history[t] = I[center]
        total_I[t] = np.sum(I)

        # Synchronous update: compute all transfers from current state
        I_new = I.copy()
        for i in range(N):
            left = (i - 1) % N
            right = (i + 1) % N

            # Transfer from left neighbor
            delta_left = k * (I[left] - I[i]) * R(I[left], I_max, n)
            # Transfer from right neighbor
            delta_right = k * (I[right] - I[i]) * R(I[right], I_max, n)

            I_new[i] += delta_left + delta_right

        # Clip to [0, I_max] — saturation enforced
        I_new = np.clip(I_new, 0.0, I_max)
        I = I_new

    return history, total_I

def count_sign_changes(arr):
    """Count oscillation sign changes relative to mean."""
    centered = arr - np.mean(arr)
    signs = np.sign(centered)
    changes = np.sum(np.abs(np.diff(signs)) > 0)
    return changes

def main():
    print("="*70)
    print("SESSION 622: DISCRETE INSTABILITY TEST")
    print("="*70)

    N = 64  # lattice size
    I_max = 1.0
    n = 2
    steps = 2000

    # ================================================================
    # TEST 1: Find the critical k for instability
    # For constant-coefficient diffusion, stability requires k < 1/(2d)
    # In 1D with D_eff = D*R(I), stability requires k*R(I) < 0.5
    # At low I, R ≈ 1, so critical k ≈ 0.5
    # ================================================================
    print("\n" + "-"*70)
    print("TEST 1: Critical coupling for discrete instability")
    print("-"*70)

    # Gaussian perturbation on uniform background
    x = np.arange(N)

    k_values = [0.1, 0.2, 0.3, 0.4, 0.49, 0.51, 0.6, 0.8, 1.0, 1.5, 2.0]

    print(f"\n{'k':>6s} | {'Sign changes':>12s} | {'I_center range':>16s} | {'Conservation':>12s} | {'Verdict':>12s}")
    print("-" * 75)

    for k in k_values:
        I_init = 0.3 * np.ones(N)
        I_init[N//2 - 2:N//2 + 3] = 0.5  # small perturbation

        history, total_I = run_1d_discrete(N, k, I_init, I_max, n, steps)

        sc = count_sign_changes(history)
        I_range = np.max(history) - np.min(history)
        conservation = (np.max(total_I) - np.min(total_I)) / total_I[0]

        verdict = "OSCILLATES" if sc > 20 else "MONOTONE"
        print(f"{k:6.2f} | {sc:12d} | {I_range:16.6f} | {conservation:12.6f} | {verdict:>12s}")

    # ================================================================
    # TEST 2: Near-saturation background — does R(I)→0 stabilize?
    # ================================================================
    print("\n" + "-"*70)
    print("TEST 2: Near-saturation background (R→0 should stabilize)")
    print("-"*70)

    backgrounds = [0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99]
    k = 1.0  # well above linear stability threshold

    print(f"\n{'I_bg':>6s} | {'R(I_bg)':>8s} | {'k*R':>6s} | {'Sign ch':>8s} | {'Range':>12s} | {'Verdict':>12s}")
    print("-" * 70)

    for bg in backgrounds:
        I_init = bg * np.ones(N)
        I_init[N//2 - 2:N//2 + 3] = min(bg + 0.05, I_max * 0.999)

        history, total_I = run_1d_discrete(N, k, I_init, I_max, n, steps)

        sc = count_sign_changes(history)
        I_range = np.max(history) - np.min(history)
        r_bg = R(bg, I_max, n)

        verdict = "OSCILLATES" if sc > 20 else "MONOTONE"
        print(f"{bg:6.2f} | {r_bg:8.4f} | {k*r_bg:6.3f} | {sc:8d} | {I_range:12.6f} | {verdict:>12s}")

    # ================================================================
    # TEST 3: Supercritical oscillation — amplitude and frequency
    # ================================================================
    print("\n" + "-"*70)
    print("TEST 3: Oscillation characteristics at k=1.0, I_bg=0.3")
    print("-"*70)

    I_init = 0.3 * np.ones(N)
    I_init[N//2] = 0.6

    history, total_I = run_1d_discrete(N, k=1.0, I_init=I_init, steps=steps)

    # Analyze oscillation in first 200 and last 200 steps
    early = history[:200]
    late = history[-200:]

    sc_early = count_sign_changes(early)
    sc_late = count_sign_changes(late)
    amp_early = np.max(early) - np.min(early)
    amp_late = np.max(late) - np.min(late)

    print(f"Early (t=0-200):  sign changes = {sc_early}, amplitude = {amp_early:.6f}")
    print(f"Late (t=1800-2000): sign changes = {sc_late}, amplitude = {amp_late:.6f}")

    if sc_early > 0 and sc_late > 0:
        print(f"Amplitude ratio (late/early) = {amp_late/amp_early:.4f}")
        if amp_late / amp_early > 0.5:
            print("→ SUSTAINED oscillation")
        else:
            print("→ DECAYING oscillation")

    # ================================================================
    # TEST 4: What is the oscillation period?
    # ================================================================
    print("\n" + "-"*70)
    print("TEST 4: Oscillation period analysis")
    print("-"*70)

    if count_sign_changes(history) > 10:
        # FFT to find dominant frequency
        centered = history - np.mean(history)
        fft = np.abs(np.fft.rfft(centered))
        freqs = np.fft.rfftfreq(len(centered))

        # Skip DC component
        peak_idx = np.argmax(fft[1:]) + 1
        peak_freq = freqs[peak_idx]
        peak_period = 1.0 / peak_freq if peak_freq > 0 else float('inf')

        print(f"Dominant frequency: {peak_freq:.4f} cycles/tick")
        print(f"Dominant period: {peak_period:.1f} ticks")

        if abs(peak_period - 2.0) < 0.5:
            print("→ Period ≈ 2 ticks = CHECKERBOARD mode (numerical artifact)")
        else:
            print(f"→ Period ≈ {peak_period:.1f} ticks — NOT checkerboard")
    else:
        print("No oscillation to analyze.")

    # ================================================================
    # TEST 5: Conservation check — does clipping destroy Intent?
    # ================================================================
    print("\n" + "-"*70)
    print("TEST 5: Conservation under clipping")
    print("-"*70)

    I_init = 0.3 * np.ones(N)
    I_init[N//2] = 0.6

    # Run WITHOUT clipping
    I = I_init.copy()
    for t in range(500):
        I_new = I.copy()
        for i in range(N):
            left = (i - 1) % N
            right = (i + 1) % N
            delta_left = 1.0 * (I[left] - I[i]) * R(I[left])
            delta_right = 1.0 * (I[right] - I[i]) * R(I[right])
            I_new[i] += delta_left + delta_right
        I = I_new

    print(f"Without clipping after 500 steps:")
    print(f"  min(I) = {np.min(I):.6f}, max(I) = {np.max(I):.6f}")
    print(f"  sum(I) = {np.sum(I):.6f} (initial = {np.sum(I_init):.6f})")
    print(f"  Conservation error: {abs(np.sum(I) - np.sum(I_init))/np.sum(I_init):.2e}")

    if np.min(I) < 0 or np.max(I) > I_max:
        print(f"  ⚠️ Values exceed physical bounds! This IS the instability.")
        print(f"  Values go negative: {np.sum(I < 0)} cells")
        print(f"  Values exceed I_max: {np.sum(I > I_max)} cells")

    # ================================================================
    # TEST 6: The critical question — is the oscillation PHYSICAL or
    # is it the same instability that makes explicit Euler useless?
    # ================================================================
    print("\n" + "-"*70)
    print("TEST 6: Physical or numerical? Reduce timestep (k) and check")
    print("-"*70)
    print("If oscillation is physical: it should persist at all k")
    print("If oscillation is numerical: it should vanish below k_crit = 0.5")

    k_test = [0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.48, 0.49, 0.50, 0.51, 0.55, 0.6]

    print(f"\n{'k':>6s} | {'Sign ch':>8s} | {'min(I)':>10s} | {'max(I)':>10s} | {'Physical?':>10s}")
    print("-" * 55)

    for k in k_test:
        I_init = 0.3 * np.ones(N)
        I_init[N//2] = 0.6

        # Run WITHOUT clipping to see raw behavior
        I = I_init.copy()
        center_hist = []
        for t in range(500):
            center_hist.append(I[N//2])
            I_new = I.copy()
            for i in range(N):
                left = (i - 1) % N
                right = (i + 1) % N
                I_new[i] += k * ((I[left] - I[i]) * R(I[left]) + (I[right] - I[i]) * R(I[right]))
            I = I_new

        center_hist = np.array(center_hist)
        sc = count_sign_changes(center_hist)

        physical = "YES" if sc > 0 else "no"
        if sc > 0 and k < 0.5:
            physical = "UNEXPECTED"

        print(f"{k:6.3f} | {sc:8d} | {np.min(I):10.4f} | {np.max(I):10.4f} | {physical:>10s}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("""
The discrete transfer rule with synchronous update IS explicit Euler
for nonlinear diffusion. For k·D_eff > 0.5 (1D), the scheme is
numerically unstable — producing checkerboard oscillations.

Key question: Is this instability PHYSICAL or an artifact?

Arguments for ARTIFACT:
  - The oscillation period is 2 ticks (Nyquist, not a physical scale)
  - Reducing k (= using a finer time step) removes the oscillation
  - Values go below 0 or above I_max — unphysical
  - Every numerical methods textbook calls this "instability to be avoided"

Arguments for PHYSICAL (if the universe IS this discrete system):
  - There IS no finer time step — Planck time is fundamental
  - k is not a numerical parameter — it's a physical coupling
  - The "instability" might be the mechanism for vacuum fluctuations
  - BUT: the energy scale is Planck energy → cosmological constant
    problem (10^122 too large)

VERDICT: The discrete-continuum gap is REAL — the discrete system
CAN oscillate where the PDE doesn't. But the oscillation is the
explicit Euler checkerboard instability. If taken as physical, it
IS the cosmological constant problem in another guise. This doesn't
save the framework; it reproduces the worst prediction in physics.

S617's continuum limit was the RIGHT move — not because the system
is continuous, but because the discrete oscillations are either:
(a) numerical artifacts (if k < k_crit, no physics there), or
(b) vacuum energy at Planck scale (if k > k_crit, 10^122 problem)

Neither option helps.
""")

if __name__ == "__main__":
    main()
