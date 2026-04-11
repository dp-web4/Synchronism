"""
Session 626 Part 2: Domain Wall Dynamics

Part 1 found that MRH-motivated dispersion produces Cahn-Hilliard
phase separation. Cell 192 showed oscillation with T=13.

Is this:
(a) Genuine domain-wall oscillation (periodic, self-sustaining)
(b) Ostwald ripening (slow, aperiodic domain coarsening)
(c) A clipping artifact

Also: what does the spatial structure look like?
Are there "entities" (localized oscillating peaks) or
"domain walls" (boundaries between uniform regions)?
"""

import numpy as np

I_MAX = 1.0


def R_monotonic(I, n=2):
    return np.maximum(0, 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n)


def step_dispersive(I, k1, k2, R_func=R_monotonic, **R_kwargs):
    I_l1 = np.roll(I, 1)
    I_r1 = np.roll(I, -1)
    I_l2 = np.roll(I, 2)
    I_r2 = np.roll(I, -2)

    R_l1 = R_func(I_l1, **R_kwargs)
    R_r1 = R_func(I_r1, **R_kwargs)
    R_l2 = R_func(I_l2, **R_kwargs)
    R_r2 = R_func(I_r2, **R_kwargs)

    dI_nn = k1 * ((I_l1 - I) * R_l1 + (I_r1 - I) * R_r1)
    dI_nnn = k2 * ((I_l2 - I) * R_l2 + (I_r2 - I) * R_r2)

    return np.clip(I + dI_nn + dI_nnn, 0.0, I_MAX)


N = 256
np.random.seed(42)
I0 = np.random.uniform(0.1, 0.9, N)

# ============================================================
print("=" * 70)
print("DOMAIN WALL ANALYSIS: k1=0.3, k2=-0.15, monotonic R")
print("=" * 70)

I = I0.copy()

# Run to phase separation
for t in range(5000):
    I = step_dispersive(I, 0.3, -0.15)

# Print spatial profile
print("\nSpatial profile at t=5000 (every 8th cell):")
for i in range(0, N, 8):
    bar = '#' * int(I[i] * 40)
    print(f"  [{i:3d}] {I[i]:.3f} {bar}")

# Find domain walls (cells where gradient is large)
grad = np.abs(np.diff(I))
wall_cells = np.where(grad > 0.1)[0]
print(f"\nDomain wall locations (|∇I| > 0.1): {wall_cells}")
print(f"Number of domain walls: {len(wall_cells)}")

# ============================================================
# Track domain wall dynamics for 10k more steps
# ============================================================
print("\n" + "=" * 70)
print("DOMAIN WALL TIME SERIES (10k steps)")
print("=" * 70)

# Record every cell every 10 steps
T_record = 10000
sample_rate = 10
n_samples = T_record // sample_rate

# Full space-time record (sampled)
spacetime = np.zeros((N, n_samples))

for t in range(T_record):
    if t % sample_rate == 0:
        spacetime[:, t // sample_rate] = I
    I = step_dispersive(I, 0.3, -0.15)

# Find cells that CHANGE (domain walls)
cell_range = np.max(spacetime, axis=1) - np.min(spacetime, axis=1)
dynamic_cells = np.where(cell_range > 0.01)[0]
print(f"Cells with range > 0.01: {len(dynamic_cells)} out of {N}")
if len(dynamic_cells) > 0:
    print(f"  Locations: {dynamic_cells[:20]}...")

    # For each dynamic cell, check if oscillation is periodic
    for c in dynamic_cells[:5]:
        ts = spacetime[c, :]
        ts_range = np.max(ts) - np.min(ts)
        ts_mean = np.mean(ts)

        # FFT
        ts_c = ts - ts_mean
        fft = np.abs(np.fft.rfft(ts_c))
        fft[0] = 0
        freqs = np.fft.rfftfreq(len(ts_c))

        if np.max(fft) > 0:
            peak_idx = np.argmax(fft)
            peak_freq = freqs[peak_idx]
            peak_power = fft[peak_idx]
            total_power = np.sum(fft)
            spectral_purity = peak_power / total_power if total_power > 0 else 0

            if peak_freq > 0:
                period = 1.0 / peak_freq * sample_rate  # Convert to actual steps
            else:
                period = np.inf

            print(f"\n  Cell {c}: range={ts_range:.4f}, period={period:.0f} steps")
            print(f"    Spectral purity: {spectral_purity:.3f} (1.0 = pure sine, 0 = noise)")
            print(f"    Is periodic: {'YES' if spectral_purity > 0.3 else 'NO — aperiodic drift'}")

            # Check direction: is I increasing or decreasing over time?
            # (Ostwald ripening = monotonic drift, not oscillation)
            first_half = np.mean(ts[:n_samples//2])
            second_half = np.mean(ts[n_samples//2:])
            drift = second_half - first_half
            print(f"    Mean drift: {drift:+.4f} ({'growing' if drift > 0.001 else 'shrinking' if drift < -0.001 else 'stable'})")
        else:
            print(f"\n  Cell {c}: no dynamics")
else:
    print("  ALL cells are static after phase separation.")


# ============================================================
# Check: how many domains, and are they coarsening?
# ============================================================
print("\n" + "=" * 70)
print("DOMAIN COARSENING CHECK")
print("=" * 70)

I = I0.copy()

checkpoints = [100, 500, 1000, 2000, 5000, 10000, 20000]
for target in checkpoints:
    while True:
        # Count steps already done
        break
    # Re-run from scratch to each checkpoint
I = I0.copy()
prev_n_domains = 0
for target in checkpoints:
    while True:
        break

# Actually just run forward and count domains at checkpoints
I = I0.copy()
step = 0
results = []
for target in checkpoints:
    while step < target:
        I = step_dispersive(I, 0.3, -0.15)
        step += 1

    # Count domains: number of transitions between "high" and "low"
    threshold = np.mean(I)
    above = I > threshold
    transitions = np.sum(np.abs(np.diff(above.astype(int))))
    n_domains = transitions // 2  # Each domain has two edges

    results.append((target, n_domains, np.mean(I), np.std(I)))

print(f"\n  {'Step':>8s} {'Domains':>8s} {'mean(I)':>8s} {'std(I)':>8s}")
print(f"  {'-'*36}")
for step, nd, mi, si in results:
    print(f"  {step:8d} {nd:8d} {mi:8.4f} {si:8.4f}")

if results[-1][1] < results[0][1]:
    print(f"\n  COARSENING: domains {results[0][1]} → {results[-1][1]}")
    print(f"  This IS Ostwald ripening (standard Cahn-Hilliard).")
    print(f"  The 'oscillation' at cell 192 was a domain wall passing through.")
elif results[-1][1] == results[0][1]:
    print(f"\n  STABLE: domain count unchanged ({results[-1][1]})")
else:
    print(f"\n  DOMAINS INCREASING: {results[0][1]} → {results[-1][1]}")


# ============================================================
# The tension: MRH vs nearest-neighbor
# ============================================================
print("\n" + "=" * 70)
print("THE MRH-COUPLING TENSION")
print("=" * 70)
print("""
FUNDAMENTALS.md Foundation 1: "all cells simultaneously evaluate
NEIGHBOR tensions and step forward together"

MRH definition: "The boundary beyond which correlations become
negligible for a given pattern's dynamics"

These are in conflict because:

1. Nearest-neighbor = FIXED coupling range (1 Planck length)
2. MRH = VARIABLE coupling range (depends on the entity)

If MRH is taken seriously:
  → coupling range varies by entity
  → some entities couple beyond nearest neighbor
  → dispersion terms appear in the dynamics
  → Cahn-Hilliard phase separation occurs
  → domains form (not entities)

If nearest-neighbor is taken literally:
  → MRH has no dynamical role at the grid scale
  → MRH is ONLY a coarse-graining concept
  → MRH doesn't affect fundamental dynamics
  → entities must form from nearest-neighbor coupling alone
  → but they can't (S617-625)

Either way, entities don't form:
  - With MRH-dispersion: domains, not entities (Cahn-Hilliard)
  - Without MRH-dispersion: diffusion, no structure (S617-625)

This is the NINTH structural impossibility, and the first that
is INTERNAL to the framework (not a comparison with external physics).
The framework's own concept (MRH) contradicts its own dynamics
(nearest-neighbor coupling), and NEITHER resolution produces entities.
""")

# ============================================================
# Total count of structural impossibilities
# ============================================================
print("=" * 70)
print("STRUCTURAL IMPOSSIBILITY SCORECARD (S617-626)")
print("=" * 70)
print("""
  #  Session  Impossibility                        Root cause
  1  S617     No oscillation from transfer rule     1-DOF diffusion
  2  S618     No wave propagation (c² < 0)          Inverted EOS
  3  S619     No gravity + waves from R(I)          P(ρ) no-go theorem
  4  S619     No cosmic acceleration                P ≥ 0 always
  5  S620     No phase dynamics (no synchronization) Real field, no phase
  6  S621     No novel prediction capacity           Unfalsifiable core + falsified math
  7  S622     No dark energy from I_max              Saturation duality theorem
  8  S625     No spatial + temporal coexistence      1-DOF exclusion
  9  S626     MRH vs coupling: no entities either way  Internal framework contradiction

  + 7 independent self-confinement failures (S19-S624)

  TOTAL: 9 structural impossibilities + 7 confinement failures = 16 independent
  proofs that the stated framework cannot produce entity-like structures.
""")
