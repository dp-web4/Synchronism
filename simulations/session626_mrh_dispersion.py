"""
Session 626: MRH-Motivated Dispersion

FUNDAMENTALS.md commits to nearest-neighbor coupling.
MRH says the relevant scale depends on the entity.
These are in tension.

If MRH introduces scale-dependent coupling, the transfer rule
gains beyond-nearest-neighbor terms. In the continuum limit,
this adds dispersion (∇⁴I). The resulting PDE is:

  ∂I/∂t = D₁·∇·[R(I)·∇I] - D₂·∇⁴I

This is a Cahn-Hilliard-like equation. It CAN produce:
  - Spinodal decomposition (phase separation into domains)
  - Localized structures (humps and drops)
  - Oscillatory instabilities (when R is non-monotonic)

Tests:
1. Does dispersion + R(I) produce localized structures?
2. If so, do they oscillate?
3. What determines their size? (MRH scale prediction)
4. Is this actually Cahn-Hilliard or something new?

The key difference from S617-625: those sessions tested the
nearest-neighbor rule, which gives DIFFUSION (no localization).
This session tests what MRH IMPLIES about the coupling range.
"""

import numpy as np

I_MAX = 1.0


def R_monotonic(I, n=2):
    return np.maximum(0, 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n)


def R_nonmonotonic(I, n=2, A=1.0):
    x = np.clip(I, 0, I_MAX) / I_MAX
    base = 1.0 - x**n
    bump = 1.0 + A * np.sin(np.pi * x)
    return np.maximum(0, base * bump)


def step_dispersive(I, k1, k2, R_func, **R_kwargs):
    """Transfer rule with nearest + next-nearest neighbor coupling.

    ΔI = k1·Σ_nn(I_n - I)·R(I_n) + k2·Σ_nnn(I_n - I)·R(I_n)

    k1 > 0: diffusive (smoothing) — original term
    k2 < 0: dispersive (structure-forming) — MRH-motivated term

    In continuum: ∂I/∂t = k1·∇·[R·∇I] + k2·∇·[R·∇³I]
    """
    I_l1 = np.roll(I, 1)   # nearest left
    I_r1 = np.roll(I, -1)  # nearest right
    I_l2 = np.roll(I, 2)   # next-nearest left
    I_r2 = np.roll(I, -2)  # next-nearest right

    R_l1 = R_func(I_l1, **R_kwargs)
    R_r1 = R_func(I_r1, **R_kwargs)
    R_l2 = R_func(I_l2, **R_kwargs)
    R_r2 = R_func(I_r2, **R_kwargs)

    # Nearest neighbor (diffusion)
    dI_nn = k1 * ((I_l1 - I) * R_l1 + (I_r1 - I) * R_r1)

    # Next-nearest neighbor (dispersion)
    dI_nnn = k2 * ((I_l2 - I) * R_l2 + (I_r2 - I) * R_r2)

    I_new = I + dI_nn + dI_nnn
    return np.clip(I_new, 0.0, I_MAX)


def spatial_entropy(I, n_bins=20):
    hist, _ = np.histogram(I, bins=n_bins, range=(0, I_MAX))
    p = hist / hist.sum()
    p = p[p > 0]
    return -np.sum(p * np.log2(p))


def correlation_length(I, max_lag=50):
    I_centered = I - np.mean(I)
    var = np.var(I)
    if var < 1e-15:
        return 0
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(I_centered * np.roll(I_centered, lag)) / var
    for i in range(1, len(acf)):
        if acf[i] < 0:
            return np.sum(acf[:i])
    return np.sum(np.maximum(acf, 0))


def count_peaks(I, threshold=None):
    """Count localized peaks (potential entities)."""
    if threshold is None:
        threshold = np.mean(I) + 2 * np.std(I)
    above = I > threshold
    # Count connected regions
    peaks = 0
    in_peak = False
    for v in above:
        if v and not in_peak:
            peaks += 1
            in_peak = True
        elif not v:
            in_peak = False
    return peaks


# ============================================================
print("=" * 70)
print("SESSION 626: MRH-MOTIVATED DISPERSION")
print("=" * 70)
print()
print("Nearest-neighbor = diffusion (S617). MRH implies beyond-nearest-neighbor.")
print("Beyond-nearest = dispersion. Dispersion + nonlinearity = localization?")
print()

N = 256
np.random.seed(42)
I0 = np.random.uniform(0.1, 0.9, N)

# ============================================================
# Test 1: Monotonic R + dispersion — does structure form?
# ============================================================
print("=" * 70)
print("TEST 1: MONOTONIC R + DISPERSION")
print("=" * 70)
print("k1 > 0 (diffusion), k2 < 0 (dispersion/anti-diffusion at nnn scale)")
print()

print(f"{'k1':>5s} {'k2':>6s} {'S_final':>8s} {'ξ':>8s} {'peaks':>6s} {'std':>8s} {'range':>14s}")
print("-" * 65)

for k1 in [0.3, 0.5]:
    for k2 in [0.0, -0.05, -0.1, -0.15, -0.2, -0.25]:
        I = I0.copy()
        stable = True

        for t in range(5000):
            I_new = step_dispersive(I, k1, k2, R_monotonic, n=2)
            if np.any(np.isnan(I_new)):
                stable = False
                break
            I = I_new

        if stable:
            S = spatial_entropy(I)
            xi = correlation_length(I)
            pk = count_peaks(I)
            std_I = np.std(I)
            print(f"  {k1:4.1f} {k2:6.2f} {S:8.3f} {xi:8.3f} {pk:6d} {std_I:8.4f} [{np.min(I):.3f}, {np.max(I):.3f}]")
        else:
            print(f"  {k1:4.1f} {k2:6.2f}  UNSTABLE at t={t}")


# ============================================================
# Test 2: Non-monotonic R + dispersion — the full combination
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: NON-MONOTONIC R + DISPERSION")
print("=" * 70)
print("This combines S624's non-monotonic R (gives chaos) with")
print("MRH-motivated dispersion (gives structure). Both together?")
print()

print(f"{'k1':>5s} {'k2':>6s} {'A':>4s} {'S':>8s} {'ξ':>8s} {'peaks':>6s} {'std':>8s} {'range':>14s}")
print("-" * 70)

for A in [0.5, 1.0, 1.5]:
    for k1 in [0.3, 0.5]:
        for k2 in [-0.05, -0.1, -0.15, -0.2]:
            I = I0.copy()
            stable = True

            for t in range(5000):
                I_new = step_dispersive(I, k1, k2, R_nonmonotonic, n=2, A=A)
                if np.any(np.isnan(I_new)):
                    stable = False
                    break
                I = I_new

            if stable:
                S = spatial_entropy(I)
                xi = correlation_length(I)
                pk = count_peaks(I)
                std_I = np.std(I)
                print(f"  {k1:4.1f} {k2:6.2f} {A:4.1f} {S:8.3f} {xi:8.3f} {pk:6d} {std_I:8.4f} [{np.min(I):.3f}, {np.max(I):.3f}]")
            else:
                print(f"  {k1:4.1f} {k2:6.2f} {A:4.1f}  UNSTABLE at t={t}")


# ============================================================
# Test 3: Detailed analysis of best cases — do structures oscillate?
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: DO STRUCTURES OSCILLATE?")
print("=" * 70)

# Find the most structured stable case from above
# Run it for longer and track temporal dynamics

best_configs = [
    (0.3, -0.15, R_monotonic, {'n': 2}, "Mono k1=0.3 k2=-0.15"),
    (0.3, -0.15, R_nonmonotonic, {'n': 2, 'A': 1.0}, "NonMono A=1 k1=0.3 k2=-0.15"),
    (0.5, -0.10, R_monotonic, {'n': 2}, "Mono k1=0.5 k2=-0.10"),
    (0.5, -0.10, R_nonmonotonic, {'n': 2, 'A': 1.0}, "NonMono A=1 k1=0.5 k2=-0.10"),
]

for k1, k2, R_func, R_kwargs, label in best_configs:
    I = I0.copy()

    # Run to steady state
    for t in range(3000):
        I = step_dispersive(I, k1, k2, R_func, **R_kwargs)

    # Record time series at several cells
    T_record = 2000
    cell_ts = {N//4: [], N//2: [], 3*N//4: []}

    for t in range(T_record):
        for c in cell_ts:
            cell_ts[c].append(I[c])
        I = step_dispersive(I, k1, k2, R_func, **R_kwargs)

    print(f"\n  {label}:")

    for c in cell_ts:
        ts = np.array(cell_ts[c])
        ts_range = np.max(ts) - np.min(ts)
        ts_std = np.std(ts)

        # FFT
        ts_centered = ts - np.mean(ts)
        if np.std(ts_centered) > 1e-12:
            fft = np.abs(np.fft.rfft(ts_centered))
            fft[0] = 0
            freqs = np.fft.rfftfreq(len(ts_centered))
            peak_idx = np.argmax(fft)
            if freqs[peak_idx] > 0:
                period = 1.0 / freqs[peak_idx]
            else:
                period = np.inf
            oscillating = ts_range > 0.01
        else:
            period = np.inf
            oscillating = False

        status = f"OSCILLATING (T={period:.0f}, amp={ts_range:.4f})" if oscillating else "STATIC"
        print(f"    Cell {c:3d}: mean={np.mean(ts):.4f} range={ts_range:.6f} → {status}")

    # Spatial structure
    xi = correlation_length(I)
    pk = count_peaks(I)
    print(f"    Spatial: ξ={xi:.2f}, peaks={pk}")


# ============================================================
# Test 4: MRH scale prediction
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: MRH SCALE FROM DISPERSION")
print("=" * 70)
print("If dispersion introduces a length scale ℓ = sqrt(|k2|/k1),")
print("structures should have characteristic size ~ ℓ.")

for k1 in [0.3, 0.5]:
    for k2 in [-0.05, -0.1, -0.15, -0.2]:
        ell = np.sqrt(abs(k2) / k1)

        I = I0.copy()
        for t in range(5000):
            I = step_dispersive(I, k1, k2, R_monotonic, n=2)

        xi = correlation_length(I)
        pk = count_peaks(I)
        mean_size = N / max(pk, 1)

        print(f"  k1={k1:.1f} k2={k2:.2f}: ℓ_pred={ell:.3f} ξ_meas={xi:.3f} "
              f"peaks={pk} mean_size={mean_size:.1f}")


# ============================================================
# Test 5: THE FOUNDATIONAL QUESTION
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: IS THIS CAHN-HILLIARD OR SOMETHING NEW?")
print("=" * 70)

print("""
The dispersive transfer rule in continuum:
  ∂I/∂t = k1·∇·[R(I)·∇I] - |k2|·∇⁴I

This IS the Cahn-Hilliard equation if:
  - R(I) plays the role of chemical potential derivative
  - k2 is the gradient energy coefficient
  - The system undergoes spinodal decomposition

Cahn-Hilliard is standard physics (phase separation in binary mixtures,
polymer blends, metallic alloys). It was derived in 1958.

The Synchronism version differs in:
  1. R(I) = [1-(I/I_max)^n] — monotonic saturation (not a typical free energy)
  2. The equation is on a FIXED grid (Planck scale) — not a continuum
  3. Intent conservation — same as mass conservation in Cahn-Hilliard

If the structures that form are Cahn-Hilliard domains, they are:
  - STATIC (no oscillation)
  - MACROSCOPIC (size ~ MRH scale)
  - SEPARATED into high-I and low-I phases

These are NOT entities (no oscillation, no self-witnessing).
They are domain walls — borders between phases.
""")

# Check: are the structures static or dynamic?
print("Checking domain dynamics over 10k steps...")

I = I0.copy()
for t in range(5000):
    I = step_dispersive(I, 0.3, -0.15, R_monotonic, n=2)

I_snapshot_5k = I.copy()

for t in range(5000):
    I = step_dispersive(I, 0.3, -0.15, R_monotonic, n=2)

I_snapshot_10k = I.copy()

diff = np.max(np.abs(I_snapshot_10k - I_snapshot_5k))
print(f"  Max change between step 5000 and 10000: {diff:.2e}")
if diff < 1e-8:
    print(f"  → STATIC. Structures are frozen. Cahn-Hilliard domains, not entities.")
elif diff < 0.01:
    print(f"  → NEARLY STATIC. Very slow coarsening (Ostwald ripening).")
else:
    print(f"  → DYNAMIC. Structures are evolving. Potential for oscillation.")

# Also check conservation
total_5k = np.sum(I_snapshot_5k)
total_10k = np.sum(I_snapshot_10k)
print(f"  Conservation: total I changed by {abs(total_10k - total_5k):.2e}")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
