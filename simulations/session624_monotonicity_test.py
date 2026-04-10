"""
Session 624: The Monotonicity Constraint

ROOT CAUSE HYPOTHESIS: Sessions 617-623 identified six independent structural failures
of the Synchronism transfer rule. All trace back to monotonic R(I). Monotonic saturation
means the CA is a contraction mapping (always smoothing). This session tests whether
non-monotonic R(I) rescues any of the failed properties.

Test 1: Wolfram class — does non-monotonic R give Class 3/4?
Test 2: Signal propagation — do perturbations propagate as waves (not diffuse)?
Test 3: Self-confinement — do localized structures persist?
Test 4: Gravity + waves — does non-monotonic R give a P(ρ) with both?

The non-monotonic R:
  R_nm(I) = [1 - (I/I_max)^n] * [1 + A * sin(π * I/I_max)]

For A > 0: R has a bump (increases then decreases). This means transfer
is ENHANCED at intermediate I — the opposite of pure saturation.
Physics: the medium becomes MORE transmissive at intermediate density
before saturating at high density. Like a phase transition in conductivity.
"""

import numpy as np
from collections import Counter
import sys

I_MAX = 1.0


def R_monotonic(I, n=2):
    """Original monotonic saturation."""
    return np.maximum(0, 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n)


def R_nonmonotonic(I, n=2, A=0.5):
    """Non-monotonic: bump at intermediate I."""
    x = np.clip(I, 0, I_MAX) / I_MAX
    base = 1.0 - x**n
    bump = 1.0 + A * np.sin(np.pi * x)
    return np.maximum(0, base * bump)


def step_1d(I, k, R_func, **R_kwargs):
    """One synchronous update step in 1D with periodic boundaries."""
    I_left = np.roll(I, 1)
    I_right = np.roll(I, -1)
    R_left = R_func(I_left, **R_kwargs)
    R_right = R_func(I_right, **R_kwargs)
    dI = k * ((I_left - I) * R_left + (I_right - I) * R_right)
    I_new = np.clip(I + dI, 0.0, I_MAX)
    return I_new


def spatial_entropy(I, n_bins=20):
    """Shannon entropy of the spatial distribution."""
    hist, _ = np.histogram(I, bins=n_bins, range=(0, I_MAX))
    p = hist / hist.sum()
    p = p[p > 0]
    return -np.sum(p * np.log2(p))


def lyapunov_test(I, k, R_func, steps=500, eps=1e-8, **R_kwargs):
    """Test sensitivity to perturbation."""
    I2 = I.copy()
    I2[len(I)//2] += eps

    for _ in range(steps):
        I = step_1d(I, k, R_func, **R_kwargs)
        I2 = step_1d(I2, k, R_func, **R_kwargs)

    diff = np.max(np.abs(I - I2))
    return np.log10(diff / eps) if diff > 0 else -20


def count_distinct_patterns(I, n_levels=10, block_size=3):
    """Count distinct patterns in discretized field."""
    levels = np.clip((I / I_MAX * n_levels).astype(int), 0, n_levels - 1)
    blocks = set()
    for i in range(len(levels) - block_size + 1):
        blocks.add(tuple(levels[i:i + block_size]))
    return len(blocks)


def signal_propagation_test(I_bg, k, R_func, steps=200, **R_kwargs):
    """Test whether a localized perturbation propagates as wave or diffuses."""
    N = len(I_bg)
    I = I_bg.copy()
    I[N // 2] = min(I[N // 2] + 0.3, I_MAX)

    # Track the width of the perturbation over time
    widths = []
    max_vals = []

    for t in range(steps):
        diff = np.abs(I - I_bg)
        if np.max(diff) > 1e-10:
            # Width = number of cells where perturbation > 10% of max
            threshold = 0.1 * np.max(diff)
            width = np.sum(diff > threshold)
            widths.append(width)
            max_vals.append(np.max(diff))
        else:
            widths.append(N)
            max_vals.append(0)
        I = step_1d(I, k, R_func, **R_kwargs)
        # Also evolve the background
        I_bg = step_1d(I_bg, k, R_func, **R_kwargs)

    return widths, max_vals


def confinement_test(k, R_func, N=256, steps=2000, **R_kwargs):
    """Test whether a localized high-I pulse stays localized."""
    I = np.ones(N) * 0.1
    # Gaussian pulse
    x = np.arange(N)
    I += 0.7 * np.exp(-((x - N//2)**2) / (2 * 5**2))
    I = np.clip(I, 0, I_MAX)

    initial_width = np.sum(I > 0.3)

    widths = []
    for t in range(steps):
        width = np.sum(I > 0.3)
        widths.append(width)
        I = step_1d(I, k, R_func, **R_kwargs)

    return initial_width, widths


def pressure_analysis(R_func, **R_kwargs):
    """Analyze the effective equation of state P(ρ) implied by R."""
    rho = np.linspace(0.01, 0.99, 1000)
    R_vals = R_func(rho, **R_kwargs)

    # P = integral of R(ρ) dρ (the only identification that gives waves, from S619)
    P = np.cumsum(R_vals) * (rho[1] - rho[0])

    # Sound speed squared = dP/dρ = R(ρ)
    cs2 = R_vals

    # For gravity: need dP/dρ < 0 somewhere (attractive)
    # For waves: need dP/dρ > 0 somewhere (propagating)
    # Non-monotonic R gives both if R changes sign or has different regimes

    # Check for non-monotonic P
    dP = np.diff(P)
    has_min = np.any(np.diff(np.sign(dP)) > 0)

    # Check where R < 0 (if ever)
    R_negative = np.any(R_vals < 0)

    # For gravity via different mechanism:
    # The key insight from S619 is that P = I_max - ρ gives gravity but no waves
    # and P = ∫R dρ gives waves but no gravity
    # Non-monotonic R might give a P with BOTH regimes

    return {
        'rho': rho,
        'R': R_vals,
        'P': P,
        'cs2': cs2,
        'has_P_minimum': has_min,
        'R_goes_negative': R_negative,
        'R_max': np.max(R_vals),
        'R_min': np.min(R_vals),
        'R_peak_location': rho[np.argmax(R_vals)]
    }


# ============================================================
# MAIN TESTS
# ============================================================

print("=" * 70)
print("SESSION 624: THE MONOTONICITY CONSTRAINT")
print("=" * 70)

N = 256
np.random.seed(42)
I0 = np.random.uniform(0.1, 0.9, N)

# ============================================================
# Test 1: Wolfram Class — Monotonic vs Non-monotonic
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: WOLFRAM CLASS COMPARISON")
print("=" * 70)

for A_val in [0.0, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
    if A_val == 0.0:
        R_func = R_monotonic
        R_kwargs = {'n': 2}
        label = "Monotonic"
    else:
        R_func = R_nonmonotonic
        R_kwargs = {'n': 2, 'A': A_val}
        label = f"Non-mono A={A_val}"

    for k in [0.3, 0.5, 0.7]:
        I = I0.copy()
        entropies = []
        for t in range(500):
            entropies.append(spatial_entropy(I))
            I = step_1d(I, k, R_func, **R_kwargs)

        # Classify based on entropy trajectory
        S_init = entropies[0]
        S_final = np.mean(entropies[-50:])
        S_std = np.std(entropies[-50:])

        n_patterns = count_distinct_patterns(I)
        lyap = lyapunov_test(I0.copy(), k, R_func, **R_kwargs)

        if S_final < 1.0 and S_std < 0.01:
            wclass = "Class 1 (fixed)"
        elif S_std < 0.1 and n_patterns < 20:
            wclass = "Class 2 (periodic)"
        elif lyap > 2:
            wclass = "Class 3 (chaotic)"
        elif S_std > 0.1 and n_patterns > 50:
            wclass = "Class 4? (complex)"
        else:
            wclass = f"Unclear (S={S_final:.2f}, std={S_std:.3f})"

        print(f"  {label:20s} k={k:.1f}: S_final={S_final:.3f} std={S_std:.4f} "
              f"patterns={n_patterns:4d} lyap={lyap:+.2f} → {wclass}")


# ============================================================
# Test 2: Signal Propagation
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: SIGNAL PROPAGATION (wave vs diffusion)")
print("=" * 70)

I_bg_uniform = np.ones(N) * 0.3

for A_val in [0.0, 0.5, 1.0, 2.0]:
    if A_val == 0.0:
        R_func = R_monotonic
        R_kwargs = {'n': 2}
        label = "Monotonic"
    else:
        R_func = R_nonmonotonic
        R_kwargs = {'n': 2, 'A': A_val}
        label = f"Non-mono A={A_val}"

    for k in [0.3, 0.5]:
        widths, maxes = signal_propagation_test(
            I_bg_uniform.copy(), k, R_func, steps=200, **R_kwargs
        )

        # Wave: width grows linearly, max stays roughly constant
        # Diffusion: width grows as sqrt(t), max decreases as 1/sqrt(t)
        if len(widths) > 100 and maxes[50] > 1e-8:
            width_ratio = widths[100] / max(widths[10], 1)
            max_ratio = maxes[100] / max(maxes[10], 1e-15)

            if max_ratio > 0.5 and width_ratio < 3:
                behavior = "WAVE-LIKE (amplitude preserved)"
            elif max_ratio < 0.1:
                behavior = "Diffusive (amplitude decays)"
            else:
                behavior = f"Mixed (amp_ratio={max_ratio:.3f}, w_ratio={width_ratio:.1f})"
        else:
            behavior = "Signal died"

        print(f"  {label:20s} k={k:.1f}: {behavior}")


# ============================================================
# Test 3: Self-Confinement
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: SELF-CONFINEMENT (pulse persistence)")
print("=" * 70)

for A_val in [0.0, 0.5, 1.0, 1.5, 2.0]:
    if A_val == 0.0:
        R_func = R_monotonic
        R_kwargs = {'n': 2}
        label = "Monotonic"
    else:
        R_func = R_nonmonotonic
        R_kwargs = {'n': 2, 'A': A_val}
        label = f"Non-mono A={A_val}"

    for k in [0.3, 0.5]:
        init_w, widths = confinement_test(k, R_func, **R_kwargs)
        final_w = widths[-1] if widths[-1] > 0 else 0

        # Check if the pulse dispersed or stayed
        if final_w <= init_w * 1.5 and final_w > 0:
            status = f"CONFINED (w: {init_w}→{final_w})"
        elif final_w == 0:
            status = f"Dissolved (w: {init_w}→0)"
        else:
            status = f"Dispersed (w: {init_w}→{final_w})"

        print(f"  {label:20s} k={k:.1f}: {status}")


# ============================================================
# Test 4: Equation of State Analysis
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: EQUATION OF STATE (gravity + waves compatibility)")
print("=" * 70)

for A_val in [0.0, 0.5, 1.0, 1.5, 2.0]:
    if A_val == 0.0:
        R_func = R_monotonic
        R_kwargs = {'n': 2}
        label = "Monotonic"
    else:
        R_func = R_nonmonotonic
        R_kwargs = {'n': 2, 'A': A_val}
        label = f"Non-mono A={A_val}"

    eos = pressure_analysis(R_func, **R_kwargs)

    # Key question: does R go negative anywhere? (needed for gravity in ∫R framework)
    # Or: does R have distinct regimes that could separate wave/gravity scales?
    print(f"  {label:20s}: R_min={eos['R_min']:.4f} R_max={eos['R_max']:.4f} "
          f"peak_at={eos['R_peak_location']:.3f} "
          f"P_has_minimum={eos['has_P_minimum']} "
          f"R_negative={eos['R_goes_negative']}")


# ============================================================
# Test 5: CRITICAL — Non-monotonic R with negative region
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: R WITH NEGATIVE REGION (anti-diffusion)")
print("=" * 70)
print("  If R goes negative, transfer REVERSES — clustering instead of spreading.")
print("  This could give both gravity (R<0 at high ρ) and waves (R>0 at low ρ).")

def R_phase_transition(I, n=2, I_crit=0.5):
    """R that changes sign: positive below I_crit, negative above.

    Physics: below critical density, Intent diffuses normally.
    Above critical density, Intent CLUSTERS (anti-diffusion).
    This is like a phase transition — and it naturally gives
    both wave propagation (dilute) and gravitational clustering (dense).
    """
    x = np.clip(I, 0, I_MAX) / I_MAX
    # Smooth transition from +1 at x=0 to -1 at x=1
    # Zero crossing at x = I_crit
    return 1.0 - 2.0 * (x / I_crit) ** n


for I_crit in [0.3, 0.5, 0.7]:
    R_func = R_phase_transition
    R_kwargs = {'n': 2, 'I_crit': I_crit}

    print(f"\n  I_crit = {I_crit}:")

    # Check R values
    rho = np.linspace(0, 1, 100)
    R_vals = R_func(rho, **R_kwargs)
    zero_cross = rho[np.argmin(np.abs(R_vals))]
    print(f"    R zero-crossing at ρ ≈ {zero_cross:.3f}")
    print(f"    R(0.1) = {R_func(np.array([0.1]), **R_kwargs)[0]:.3f} "
          f"R(0.5) = {R_func(np.array([0.5]), **R_kwargs)[0]:.3f} "
          f"R(0.9) = {R_func(np.array([0.9]), **R_kwargs)[0]:.3f}")

    # Test stability
    for k in [0.1, 0.3, 0.5]:
        I = I0.copy() * 0.5  # Start in the mixed regime
        stable = True
        for t in range(500):
            I_new = step_1d(I, k, R_func, **R_kwargs)
            if np.any(np.isnan(I_new)) or np.any(np.isinf(I_new)):
                print(f"    k={k}: UNSTABLE at t={t}")
                stable = False
                break
            I = I_new

        if stable:
            S = spatial_entropy(I)
            n_pat = count_distinct_patterns(I)
            # Check if clustering occurred
            I_max_val = np.max(I)
            I_min_val = np.min(I)
            contrast = I_max_val - I_min_val

            print(f"    k={k}: Stable. S={S:.3f} patterns={n_pat} "
                  f"contrast={contrast:.4f} range=[{I_min_val:.4f}, {I_max_val:.4f}]")


# ============================================================
# Test 6: THE FOUNDATIONAL QUESTION
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: FOUNDATIONAL QUESTION — WHAT DOES NON-MONOTONIC R CHANGE?")
print("=" * 70)

# Compare monotonic vs best non-monotonic across all metrics
print("\nSummary comparison (k=0.5, 500 steps, N=256):")
print(f"{'Metric':<30s} {'Monotonic':>12s} {'Non-mono A=1':>12s} {'Phase trans':>12s}")
print("-" * 70)

configs = [
    ("Monotonic", R_monotonic, {'n': 2}),
    ("Non-mono A=1", R_nonmonotonic, {'n': 2, 'A': 1.0}),
    ("Phase trans", R_phase_transition, {'n': 2, 'I_crit': 0.5}),
]

k = 0.5
for label, R_func, R_kwargs in configs:
    I = I0.copy()
    for t in range(500):
        I_new = step_1d(I, k, R_func, **R_kwargs)
        if np.any(np.isnan(I_new)):
            break
        I = I_new

    if not np.any(np.isnan(I)):
        S = spatial_entropy(I)
        n_pat = count_distinct_patterns(I)
        lyap = lyapunov_test(I0.copy(), k, R_func, steps=500, **R_kwargs)
        contrast = np.max(I) - np.min(I)
    else:
        S = n_pat = lyap = contrast = float('nan')

# Need to collect these properly
results = {}
for label, R_func, R_kwargs in configs:
    I = I0.copy()
    valid = True
    for t in range(500):
        I_new = step_1d(I, k, R_func, **R_kwargs)
        if np.any(np.isnan(I_new)):
            valid = False
            break
        I = I_new

    if valid:
        results[label] = {
            'S': spatial_entropy(I),
            'patterns': count_distinct_patterns(I),
            'lyap': lyapunov_test(I0.copy(), k, R_func, steps=500, **R_kwargs),
            'contrast': np.max(I) - np.min(I),
        }
    else:
        results[label] = {'S': 'NaN', 'patterns': 'NaN', 'lyap': 'NaN', 'contrast': 'NaN'}

for metric in ['S', 'patterns', 'lyap', 'contrast']:
    vals = []
    for label, _, _ in configs:
        v = results[label][metric]
        if isinstance(v, str):
            vals.append(v)
        elif isinstance(v, int):
            vals.append(f"{v}")
        else:
            vals.append(f"{v:.4f}")
    print(f"  {metric:<28s} {vals[0]:>12s} {vals[1]:>12s} {vals[2]:>12s}")

print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)
