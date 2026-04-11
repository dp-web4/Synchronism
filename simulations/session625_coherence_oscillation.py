"""
Session 625: Coherence vs Oscillation — Conflict or Harmony?

SESSION_FOCUS Open Question #2 treats oscillation basis and C(ρ) as
potentially complementary: "coherence C may be stability measure of
recurring pattern." But what if they CONFLICT?

Oscillation = temporal recurrence (the pattern returns after τ ticks)
Coherence = spatial correlation (nearby cells have similar values)

These are INDEPENDENT properties. A field can be:
  - Temporally periodic but spatially uniform (checkerboard) — no entity
  - Spatially structured but temporally static (wall) — no entity
  - Both periodic and structured — entity candidate

The question: Is there a TRADE-OFF between spatial coherence and
temporal stability? If so, it constrains what entities can exist.
If not, they're truly independent and the framework's attempt to
unify them through C(ρ) is wrong.

Test in the non-monotonic CA at edge-of-chaos (A=1.0, k=0.40),
which is the ONLY regime with non-trivial dynamics (S624).
"""

import numpy as np

I_MAX = 1.0


def R_nonmonotonic(I, n=2, A=0.5):
    x = np.clip(I, 0, I_MAX) / I_MAX
    base = 1.0 - x**n
    bump = 1.0 + A * np.sin(np.pi * x)
    return np.maximum(0, base * bump)


def R_monotonic(I, n=2):
    return np.maximum(0, 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n)


def step_1d(I, k, R_func, **R_kwargs):
    I_left = np.roll(I, 1)
    I_right = np.roll(I, -1)
    R_left = R_func(I_left, **R_kwargs)
    R_right = R_func(I_right, **R_kwargs)
    dI = k * ((I_left - I) * R_left + (I_right - I) * R_right)
    return np.clip(I + dI, 0.0, I_MAX)


def spatial_autocorrelation(I, max_lag=50):
    """Compute spatial autocorrelation function."""
    N = len(I)
    I_centered = I - np.mean(I)
    var = np.var(I)
    if var < 1e-15:
        return np.zeros(max_lag)

    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(I_centered * np.roll(I_centered, lag)) / var
    return acf


def correlation_length(I, max_lag=50):
    """Correlation length = integral of autocorrelation function."""
    acf = spatial_autocorrelation(I, max_lag)
    # Integrate until first zero crossing
    for i in range(1, len(acf)):
        if acf[i] < 0:
            return np.sum(acf[:i])
    return np.sum(np.maximum(acf, 0))


def temporal_spectrum(time_series):
    """FFT of a cell's time series. Returns dominant period and spectral width."""
    ts = time_series - np.mean(time_series)
    if np.std(ts) < 1e-15:
        return {'period': np.inf, 'width': 0, 'power': 0}

    fft = np.abs(np.fft.rfft(ts))
    freqs = np.fft.rfftfreq(len(ts))

    # Skip DC
    fft[0] = 0

    if np.max(fft) < 1e-15:
        return {'period': np.inf, 'width': 0, 'power': 0}

    # Dominant frequency
    peak_idx = np.argmax(fft)
    if freqs[peak_idx] > 0:
        period = 1.0 / freqs[peak_idx]
    else:
        period = np.inf

    # Spectral width (std of frequency weighted by power)
    power = fft**2
    total_power = np.sum(power)
    if total_power > 0:
        mean_freq = np.sum(freqs * power) / total_power
        var_freq = np.sum((freqs - mean_freq)**2 * power) / total_power
        width = np.sqrt(var_freq)
    else:
        width = 0

    return {'period': period, 'width': width, 'power': np.max(fft)}


# ============================================================
print("=" * 70)
print("SESSION 625: COHERENCE vs OSCILLATION")
print("=" * 70)

N = 256
np.random.seed(42)

# ============================================================
# Test 1: Correlation length vs temporal stability in edge-of-chaos
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: Coherence-Oscillation relationship at edge of chaos")
print("=" * 70)

configs = [
    ("Monotonic k=0.3", R_monotonic, {'n': 2}, 0.3),
    ("Monotonic k=0.5", R_monotonic, {'n': 2}, 0.5),
    ("Monotonic k=0.7", R_monotonic, {'n': 2}, 0.7),
    ("NonMono A=0.5 k=0.5", R_nonmonotonic, {'n': 2, 'A': 0.5}, 0.5),
    ("NonMono A=1.0 k=0.40", R_nonmonotonic, {'n': 2, 'A': 1.0}, 0.40),  # Class 4
    ("NonMono A=1.0 k=0.55", R_nonmonotonic, {'n': 2, 'A': 1.0}, 0.55),  # Class 3
    ("NonMono A=1.5 k=0.40", R_nonmonotonic, {'n': 2, 'A': 1.5}, 0.40),  # Class 4
    ("NonMono A=1.5 k=0.50", R_nonmonotonic, {'n': 2, 'A': 1.5}, 0.50),  # Class 3
    ("NonMono A=0.8 k=0.70", R_nonmonotonic, {'n': 2, 'A': 0.8}, 0.70),  # Class 3
]

print(f"\n{'Config':<30s} {'ξ_space':>8s} {'T_dom':>8s} {'Δf':>8s} {'ξ×Δf':>8s} {'Class':>8s}")
print("-" * 80)

for label, R_func, R_kwargs, k in configs:
    I0 = np.random.uniform(0.1, 0.9, N)
    I = I0.copy()

    # Burn-in
    for _ in range(500):
        I = step_1d(I, k, R_func, **R_kwargs)

    # Collect time series
    T_record = 2000
    time_series = np.zeros((N, T_record))
    for t in range(T_record):
        time_series[:, t] = I
        I = step_1d(I, k, R_func, **R_kwargs)

    # Spatial: correlation length (averaged over time)
    xi_vals = []
    for t in range(0, T_record, 50):
        xi = correlation_length(time_series[:, t])
        xi_vals.append(xi)
    xi_mean = np.mean(xi_vals)

    # Temporal: average spectral properties across cells
    periods = []
    widths = []
    for cell in range(0, N, 4):  # Sample every 4th cell
        spec = temporal_spectrum(time_series[cell, :])
        if spec['period'] < np.inf:
            periods.append(spec['period'])
            widths.append(spec['width'])

    if len(periods) > 0:
        T_dom = np.median(periods)
        Df = np.mean(widths)
        product = xi_mean * Df

        # Classify
        if Df < 0.01:
            wclass = "Static"
        elif T_dom < 3:
            wclass = "Checker"
        elif xi_mean > 10 and Df < 0.05:
            wclass = "Ordered"
        elif Df > 0.1:
            wclass = "Chaotic"
        else:
            wclass = "Complex"
    else:
        T_dom = Df = product = 0
        wclass = "Dead"

    print(f"  {label:<28s} {xi_mean:8.3f} {T_dom:8.1f} {Df:8.4f} {product:8.4f} {wclass:>8s}")


# ============================================================
# Test 2: Is there a trade-off? Scan across the edge-of-chaos boundary
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Trade-off scan across chaos boundary (A=1.0, varying k)")
print("=" * 70)

I0 = np.random.uniform(0.1, 0.9, N)

k_vals = np.arange(0.30, 0.65, 0.02)
results = []

for k in k_vals:
    I = I0.copy()

    # Burn-in
    for _ in range(500):
        I = step_1d(I, k, R_nonmonotonic, n=2, A=1.0)

    # Collect
    T_record = 1000
    time_series = np.zeros((N, T_record))
    for t in range(T_record):
        time_series[:, t] = I
        I = step_1d(I, k, R_nonmonotonic, n=2, A=1.0)

    # Spatial coherence
    xi_vals = []
    for t in range(0, T_record, 100):
        xi = correlation_length(time_series[:, t])
        xi_vals.append(xi)
    xi_mean = np.mean(xi_vals)

    # Temporal stability (inverse of spectral width)
    widths = []
    for cell in range(0, N, 8):
        spec = temporal_spectrum(time_series[cell, :])
        widths.append(spec['width'])
    Df = np.mean(widths)

    # Lyapunov-like: sensitivity
    I2 = I0.copy()
    I2[N//2] += 1e-10
    Ia = I0.copy()
    for _ in range(800):
        Ia = step_1d(Ia, k, R_nonmonotonic, n=2, A=1.0)
        I2 = step_1d(I2, k, R_nonmonotonic, n=2, A=1.0)
    diff = np.max(np.abs(Ia - I2))
    lyap = np.log10(diff / 1e-10) if diff > 0 else -20

    results.append((k, xi_mean, Df, lyap))

print(f"\n{'k':>6s} {'ξ_space':>8s} {'Δf':>8s} {'ξ×Δf':>8s} {'Lyapunov':>8s} {'Phase':>12s}")
print("-" * 60)

for k, xi, df, lyap in results:
    product = xi * df
    if lyap > 3:
        phase = "Chaotic"
    elif lyap > 0:
        phase = "Edge"
    elif xi < 2:
        phase = "Fixed"
    else:
        phase = "Periodic"
    print(f"  {k:5.2f} {xi:8.3f} {df:8.4f} {product:8.4f} {lyap:+8.2f} {phase:>12s}")


# ============================================================
# Test 3: The Heisenberg analogy
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: Is ξ × Δf ≥ constant (Heisenberg-like bound)?")
print("=" * 70)

products = [xi * df for _, xi, df, _ in results if df > 0.001]
if products:
    print(f"  Product ξ × Δf across all configs:")
    print(f"    min = {min(products):.4f}")
    print(f"    max = {max(products):.4f}")
    print(f"    mean = {np.mean(products):.4f}")
    print(f"    std = {np.std(products):.4f}")

    # Is there a lower bound?
    if min(products) > 0.01:
        print(f"\n  There IS a lower bound: ξ × Δf ≥ {min(products):.4f}")
        print(f"  But is this physical or trivial?")
        print(f"  Trivial if: Δf is bounded below by numerical precision")
        print(f"  Physical if: Δf → 0 as ξ → ∞ or vice versa (trade-off)")
    else:
        print(f"\n  No significant lower bound found.")
        print(f"  ξ and Δf are NOT in trade-off — they can both be small.")

    # Check correlation between ξ and Δf
    xi_arr = np.array([xi for _, xi, df, _ in results if df > 0.001])
    df_arr = np.array([df for _, xi, df, _ in results if df > 0.001])
    if len(xi_arr) > 3:
        corr = np.corrcoef(xi_arr, df_arr)[0, 1]
        print(f"\n  Correlation(ξ, Δf) = {corr:+.4f}")
        if corr < -0.5:
            print(f"  → TRADE-OFF: higher spatial coherence ↔ lower temporal stability")
        elif corr > 0.5:
            print(f"  → CO-OCCURRENCE: spatial and temporal structure rise together")
        else:
            print(f"  → INDEPENDENT: no systematic relationship")


# ============================================================
# Test 4: What C(ρ) actually measures
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: What does C(ρ) = tanh(ρ/ρ₀) actually measure?")
print("=" * 70)
print("C(ρ) is a function of density alone. It doesn't know about:")
print("  - Temporal dynamics (oscillation, stability)")
print("  - Spatial structure (correlation length)")
print("  - Phase relationships (coherence in the QM sense)")
print("It's a LOCAL, INSTANTANEOUS function. How does it relate to")
print("the DYNAMICAL properties measured above?")

# Compute C(ρ) for different density profiles and compare with ξ and Δf
rho_0 = 0.5  # Typical C(ρ) parameter

for k_test in [0.35, 0.40, 0.45, 0.55]:
    I = I0.copy()
    for _ in range(500):
        I = step_1d(I, k_test, R_nonmonotonic, n=2, A=1.0)

    # C(ρ) at each cell
    C_vals = np.tanh(I / rho_0)
    C_mean = np.mean(C_vals)
    C_std = np.std(C_vals)

    # Spatial coherence of C itself
    xi_C = correlation_length(C_vals)

    # Compare with xi of I
    xi_I = correlation_length(I)

    # Temporal stability
    ts = np.zeros(500)
    for t in range(500):
        ts[t] = I[N//2]
        I = step_1d(I, k_test, R_nonmonotonic, n=2, A=1.0)
    spec = temporal_spectrum(ts)

    print(f"\n  k={k_test:.2f}:")
    print(f"    C(ρ) mean={C_mean:.3f} std={C_std:.3f}")
    print(f"    ξ(I) = {xi_I:.3f}, ξ(C) = {xi_C:.3f}")
    print(f"    T_dom = {spec['period']:.1f}, Δf = {spec['width']:.4f}")
    print(f"    C(ρ) captures: {'density only' if abs(xi_C - xi_I) < 0.5 else 'different structure than I'}")


# ============================================================
# Test 5: The real question
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: DOES C(ρ) PREDICT OSCILLATION STABILITY?")
print("=" * 70)

# If C(ρ) and oscillation are "the same thing at different levels"
# (SESSION_FOCUS claim), then high C should correlate with stable oscillation.
# If they CONFLICT, high C should correlate with UNSTABLE oscillation.

# Use the edge-of-chaos regime for non-trivial dynamics
I = np.random.uniform(0.1, 0.9, N)
k_test = 0.42  # Near edge of chaos

# Burn-in
for _ in range(500):
    I = step_1d(I, k_test, R_nonmonotonic, n=2, A=1.0)

# Collect per-cell: mean C(ρ) and oscillation stability (spectral width)
T_collect = 2000
cell_data = {i: {'C_vals': [], 'I_vals': []} for i in range(N)}

I_current = I.copy()
for t in range(T_collect):
    C = np.tanh(I_current / rho_0)
    for i in range(N):
        cell_data[i]['C_vals'].append(C[i])
        cell_data[i]['I_vals'].append(I_current[i])
    I_current = step_1d(I_current, k_test, R_nonmonotonic, n=2, A=1.0)

# Per-cell: mean coherence vs temporal stability
C_means = []
delta_fs = []
for i in range(0, N, 2):
    C_mean = np.mean(cell_data[i]['C_vals'])
    spec = temporal_spectrum(np.array(cell_data[i]['I_vals']))
    C_means.append(C_mean)
    delta_fs.append(spec['width'])

C_arr = np.array(C_means)
Df_arr = np.array(delta_fs)

corr = np.corrcoef(C_arr, Df_arr)[0, 1]
print(f"\n  Correlation(mean C(ρ), spectral width Δf) = {corr:+.4f}")

if corr > 0.3:
    print(f"  → CONFLICT: Higher C(ρ) → WIDER spectrum → LESS stable oscillation")
    print(f"  → C(ρ) and oscillation stability are ANTI-CORRELATED")
elif corr < -0.3:
    print(f"  → HARMONY: Higher C(ρ) → NARROWER spectrum → MORE stable oscillation")
    print(f"  → C(ρ) predicts oscillation stability (framework claim supported)")
else:
    print(f"  → INDEPENDENT: C(ρ) tells you nothing about oscillation stability")
    print(f"  → They measure different things — combining them is a category error")

# Partition cells into high-C and low-C groups
median_C = np.median(C_arr)
high_C = Df_arr[C_arr > median_C]
low_C = Df_arr[C_arr <= median_C]
print(f"\n  High-C cells: Δf = {np.mean(high_C):.4f} ± {np.std(high_C):.4f}")
print(f"  Low-C cells:  Δf = {np.mean(low_C):.4f} ± {np.std(low_C):.4f}")
print(f"  Ratio: {np.mean(high_C)/max(np.mean(low_C), 1e-10):.3f}")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
