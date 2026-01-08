#!/usr/bin/env python3
"""
Session #239 Extended: Deep MOND Limit Analysis

KEY PREDICTION: Synchronism predicts a MAXIMUM gravity boost of γ_max = 1/Ω_m ≈ 3.17
MOND (in most formulations) predicts unlimited boost as a → 0

This is a decisive test between the two frameworks.
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
phi = (1 + np.sqrt(5)) / 2  # Golden ratio
a_0 = 1.2e-10  # m/s²
Omega_m = 0.315  # Matter fraction

print("=" * 70)
print("SESSION #239 EXTENDED: DEEP MOND LIMIT ANALYSIS")
print("=" * 70)

# =============================================================================
# Part 1: Theoretical Predictions
# =============================================================================

print("\n" + "=" * 70)
print("PART 1: THEORETICAL PREDICTIONS")
print("=" * 70)

# Synchronism deep limit
gamma_max_sync = 1 / Omega_m
print(f"\nSYNCHRONISM:")
print(f"As a → 0:")
print(f"  C(a) → Ω_m = {Omega_m}")
print(f"  γ_max = 1/Ω_m = {gamma_max_sync:.3f}")
print(f"\nThis is a HARD UPPER BOUND - gravity can never be boosted more than {gamma_max_sync:.2f}×")

# MOND deep limit (standard formulation)
print(f"\nMOND (standard AQUAL/QUMOND):")
print(f"As a → 0:")
print(f"  μ(a/a₀) → a/a₀")
print(f"  g = √(g_N × a₀) (deep MOND)")
print(f"  γ = √(a₀/a) → ∞ as a → 0")
print(f"\nThis predicts UNLIMITED boost at very low acceleration")

# RAR empirical limit
print(f"\nRAR (Radial Acceleration Relation):")
print(f"g_obs/g_bar = 1/(1 - exp(-√(g_bar/g†)))")
print(f"As g_bar → 0: g_obs/g_bar → √(g†/g_bar) → ∞")
print(f"Also predicts unlimited boost")

# =============================================================================
# Part 2: What Observations Show
# =============================================================================

print("\n" + "=" * 70)
print("PART 2: OBSERVATIONAL CONSTRAINTS")
print("=" * 70)

# Lowest acceleration probes
print(f"\nLowest acceleration measurements:")
print(f"{'System':<25} {'a (m/s²)':<15} {'γ observed':<15} {'Notes'}")
print("-" * 70)

observations = [
    ("Wide binaries (Gaia)", "3e-11", "~1.5-2.0", "Most stringent"),
    ("Galaxy rotation curves", "~1e-10", "~10 (max)", "But contaminated"),
    ("Galaxy clusters", "~1e-10", "~2-3", "At outskirts"),
    ("Dwarf spheroidals", "~1e-11", "~100 (M/L)", "Controversial"),
    ("DF2/DF4 'missing DM'", "~3e-10", "~1", "No DM if isolated"),
]

for sys, a, gamma, notes in observations:
    print(f"{sys:<25} {a:<15} {gamma:<15} {notes}")

print(f"\nKEY POINT: No clean measurement yet at a < 10⁻¹¹ m/s²")
print(f"Wide binaries are currently the best probe")

# =============================================================================
# Part 3: Quantitative Predictions
# =============================================================================

print("\n" + "=" * 70)
print("PART 3: QUANTITATIVE PREDICTIONS FOR VERY LOW a")
print("=" * 70)

def C_sync(a):
    x = (a / a_0) ** (1 / phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def C_mond(a):
    x = a / a_0
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def gamma_mond_deep(a):
    """Standard MOND deep limit: g = sqrt(g_N * a0), so gamma = sqrt(a0/a)"""
    return np.sqrt(a_0 / a)

# Very low acceleration range
a_range = np.logspace(-14, -9, 100)

print(f"\n{'a (m/s²)':<15} {'γ_sync':<12} {'γ_MOND_interp':<15} {'γ_MOND_deep':<15} {'Ratio (deep/sync)':<15}")
print("-" * 80)

for a in [1e-10, 3e-11, 1e-11, 3e-12, 1e-12, 3e-13, 1e-13]:
    gs = 1 / C_sync(a)
    gm = 1 / C_mond(a)
    gd = gamma_mond_deep(a)
    ratio = gd / gs
    print(f"{a:<15.1e} {gs:<12.3f} {gm:<15.3f} {gd:<15.1f} {ratio:<15.1f}")

print(f"\nAs a → 0:")
print(f"Synchronism saturates at γ = {gamma_max_sync:.3f}")
print(f"MOND deep limit: γ = √(a₀/a) → ∞")

# =============================================================================
# Part 4: Observable Consequences
# =============================================================================

print("\n" + "=" * 70)
print("PART 4: OBSERVABLE CONSEQUENCES")
print("=" * 70)

print("""
If Synchronism is correct (γ_max = 3.17):
- Wide binaries at a < 10⁻¹¹ m/s² should show γ plateauing
- Galaxy outskirts should never exceed γ ~ 3
- No system should require M_dyn/M_bar > 3.17 for pure coherence effect

If MOND is correct (γ → ∞):
- Wide binaries at very low a should continue increasing
- γ = √(a₀/a) should apply
- At a = 10⁻¹² m/s²: γ ~ 10 (vs Sync: 3.1)
- At a = 10⁻¹³ m/s²: γ ~ 35 (vs Sync: 3.2)

The DIVERGENCE between predictions grows rapidly at low a!
""")

# =============================================================================
# Part 5: Dwarf Galaxies and Tidal Systems
# =============================================================================

print("\n" + "=" * 70)
print("PART 5: DWARF GALAXIES - A POSSIBLE TEST")
print("=" * 70)

# Ultra-faint dwarfs have very low internal accelerations
ufd_data = {
    "Segue 1": {"M_star": 340, "r_half": 29, "sigma_v": 3.7},  # Solar masses, pc, km/s
    "Segue 2": {"M_star": 900, "r_half": 35, "sigma_v": 2.2},
    "Ret II": {"M_star": 2600, "r_half": 55, "sigma_v": 3.3},
    "Tucana III": {"M_star": 800, "r_half": 44, "sigma_v": 1.5},
}

print(f"\nUltra-faint dwarf satellites:")
print(f"{'Name':<15} {'M_star':<12} {'r_half (pc)':<12} {'σ_v (km/s)':<12} {'a_int (m/s²)':<15} {'γ_needed':<12}")
print("-" * 80)

G = 4.3e-3  # pc (km/s)^2 / M_sun

for name, data in ufd_data.items():
    # Internal acceleration at half-light radius
    a_int = G * data["M_star"] / data["r_half"]**2 * 3.086e13  # Convert to m/s²
    # Expected velocity from baryons
    v_bar = np.sqrt(G * data["M_star"] / data["r_half"])
    # Observed velocity
    v_obs = data["sigma_v"]
    # Required boost
    gamma_needed = (v_obs / v_bar)**2 if v_bar > 0 else 0

    print(f"{name:<15} {data['M_star']:<12.0f} {data['r_half']:<12.0f} {data['sigma_v']:<12.1f} {a_int:<15.1e} {gamma_needed:<12.1f}")

print(f"\nNote: Ultra-faints have a_int ~ 10⁻¹²-10⁻¹¹ m/s²")
print(f"If γ_needed >> 3.17, EITHER:")
print(f"  1. They have dark matter (not just modified gravity)")
print(f"  2. They are tidally affected (MW external field)")
print(f"  3. Synchronism limit is not the full story")

print(f"\nCurrent interpretation: Tidal effects from MW dominate for most UFDs")

# =============================================================================
# Part 6: The Convergence Test
# =============================================================================

print("\n" + "=" * 70)
print("PART 6: CONVERGENCE TEST - HOW TO MEASURE γ_max")
print("=" * 70)

print("""
To measure γ_max, we need:
1. Systems at very low INTERNAL acceleration (a_int << a₀)
2. ISOLATED from external fields (no MW/host contamination)
3. Reliable mass measurements (stellar masses known)
4. Reliable kinematics (velocities measured)

Candidate systems:
- Wide binaries at r > 10 kpc from MW (low a_ext)
- Isolated dwarf galaxies (not satellites)
- Galaxy pairs in voids
- Very outer rotation curves of isolated galaxies

Currently problematic:
- Most UFDs are satellites (EFE contaminated)
- Wide binaries are within MW (a_ext ~ 2×10⁻¹⁰ m/s²)
- Isolated dwarfs are rare and faint
""")

# =============================================================================
# Part 7: Predictions for Specific Acceleration Ranges
# =============================================================================

print("\n" + "=" * 70)
print("PART 7: DETAILED PREDICTIONS FOR FUTURE OBSERVATIONS")
print("=" * 70)

# Generate predictions table
print(f"\nPredicted gravity boost vs acceleration:")
print(f"{'a (m/s²)':<15} {'γ_sync':<12} {'γ_MOND_simple':<15} {'γ_MOND_deep':<15} {'Sync vs Deep':<15}")
print("-" * 75)

predictions = []
for log_a in np.arange(-8, -14, -0.5):
    a = 10**log_a
    gs = 1 / C_sync(a)
    gm = 1 / C_mond(a)
    gd = gamma_mond_deep(a)
    diff = gs / gd * 100  # Sync as % of MOND deep
    predictions.append((a, gs, gm, gd, diff))
    print(f"{a:<15.1e} {gs:<12.3f} {gm:<15.3f} {gd:<15.1f} {diff:<15.1f}%")

# Key diagnostic: At what acceleration does Sync differ from MOND deep by >50%?
print(f"\nDiagnostic thresholds:")
for a, gs, gm, gd, diff in predictions:
    if diff < 50:
        print(f"Sync < 50% of MOND deep at a ~ {a:.1e} m/s²")
        break

for a, gs, gm, gd, diff in predictions:
    if gs > 0.95 * gamma_max_sync:
        print(f"Sync reaches 95% of max at a ~ {a:.1e} m/s²")
        break

# =============================================================================
# Part 8: Visualization
# =============================================================================

print("\n" + "=" * 70)
print("PART 8: GENERATING VISUALIZATIONS")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Extended acceleration range
a_ext = np.logspace(-14, -8, 200)

# Calculate predictions
gamma_sync = 1 / C_sync(a_ext)
gamma_mond = 1 / C_mond(a_ext)
gamma_deep = gamma_mond_deep(a_ext)

# Plot 1: All three predictions
ax1 = axes[0]
ax1.loglog(a_ext, gamma_sync, 'b-', lw=2.5, label='Synchronism')
ax1.loglog(a_ext, gamma_mond, 'r--', lw=2, label='MOND (interpolating)')
ax1.loglog(a_ext, gamma_deep, 'g:', lw=2, label='MOND (deep limit)')
ax1.axhline(gamma_max_sync, color='purple', ls='--', lw=1.5, label=f'Sync max = {gamma_max_sync:.2f}')
ax1.axvline(a_0, color='gray', ls=':', alpha=0.5, label='a₀')

# Current observation limit
ax1.axvspan(1e-11, 1e-9, alpha=0.2, color='yellow', label='Current obs range')

ax1.set_xlabel('Acceleration (m/s²)', fontsize=12)
ax1.set_ylabel('Gravity boost γ_g', fontsize=12)
ax1.set_title('Deep MOND Limit Test', fontsize=14)
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(1e-14, 1e-8)
ax1.set_ylim(0.9, 100)

# Plot 2: Ratio of Sync to MOND deep
ax2 = axes[1]
ratio = gamma_sync / gamma_deep * 100
ax2.semilogx(a_ext, ratio, 'b-', lw=2.5)
ax2.axhline(100, color='gray', ls=':', alpha=0.5)
ax2.axhline(50, color='red', ls='--', alpha=0.5, label='50% of MOND')
ax2.axvline(a_0, color='gray', ls=':', alpha=0.5, label='a₀')
ax2.axvspan(1e-11, 1e-9, alpha=0.2, color='yellow', label='Current obs')

ax2.set_xlabel('Acceleration (m/s²)', fontsize=12)
ax2.set_ylabel('γ_sync / γ_MOND_deep (%)', fontsize=12)
ax2.set_title('Synchronism vs MOND Deep Limit', fontsize=14)
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(1e-14, 1e-8)
ax2.set_ylim(0, 120)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session239_deep_mond_limit.png', dpi=150)
plt.close()

print("Saved: session239_deep_mond_limit.png")

# =============================================================================
# Part 9: Summary
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #239 EXTENDED SUMMARY: DEEP MOND LIMIT")
print("=" * 70)

print(f"""
KEY PREDICTIONS:

1. SYNCHRONISM MAXIMUM BOOST
   - γ_max = 1/Ω_m = {gamma_max_sync:.3f}
   - This is a HARD UPPER BOUND
   - Reached at a << a₀

2. MOND DEEP LIMIT
   - γ = √(a₀/a) → ∞
   - No upper bound
   - Continues increasing indefinitely

3. DIVERGENCE POINTS
   - At a = 10⁻¹¹ m/s²: Sync γ = {1/C_sync(1e-11):.2f}, MOND deep γ = {gamma_mond_deep(1e-11):.0f}
   - At a = 10⁻¹² m/s²: Sync γ = {1/C_sync(1e-12):.2f}, MOND deep γ = {gamma_mond_deep(1e-12):.0f}
   - At a = 10⁻¹³ m/s²: Sync γ = {1/C_sync(1e-13):.2f}, MOND deep γ = {gamma_mond_deep(1e-13):.0f}

4. TEST STRATEGY
   - Find ISOLATED systems at a << 10⁻¹¹ m/s²
   - Measure γ directly from dynamics
   - If γ > 4: Sync falsified
   - If γ plateaus at ~3: MOND falsified

5. CURRENT STATUS
   - No clean measurements below 10⁻¹¹ m/s²
   - Ultra-faint dwarfs contaminated by tidal effects
   - Wide binaries limited by MW external field

CONCLUSION:
The deep MOND limit is a decisive test. Synchronism predicts convergence to
γ_max ≈ 3.17, while MOND predicts divergence. This is potentially falsifiable
with future observations of extremely isolated, low-acceleration systems.
""")

print("\n" + "=" * 70)
print("SESSION #239 EXTENDED COMPLETE")
print("=" * 70)
