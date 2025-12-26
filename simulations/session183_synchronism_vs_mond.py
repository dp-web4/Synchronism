#!/usr/bin/env python3
"""
Session #183: Synchronism vs MOND Discrimination

Both Synchronism and MOND explain TDG mass discrepancies.
This session identifies observational tests that discriminate between them.

KEY DIFFERENCE:
- MOND: Mass discrepancy depends on ACCELERATION (a < a_0)
- Synchronism: Mass discrepancy depends on DENSITY (ρ < ρ_t)

This leads to different predictions for:
1. Objects at same acceleration but different densities
2. Objects at same density but different accelerations
3. The External Field Effect (MOND predicts it, Synchronism doesn't)
4. Environment-dependence vs acceleration-dependence of BTFR scatter

DISCRIMINATING TESTS:
A) TDGs vs isolated dwarfs at same V_rot
B) Cluster dwarfs (high ρ, low a) vs field dwarfs (low ρ, low a)
C) BTFR scatter correlation with environment vs acceleration
D) EFE detection (MOND-specific prediction)
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple

# Physical constants
G = 4.30e-6  # kpc (km/s)^2 / M_sun
a_0_mond = 1.2e-10  # m/s^2 - MOND acceleration scale

# Synchronism parameters
phi = (1 + np.sqrt(5)) / 2  # Golden ratio
Omega_m = 0.315

def mond_mu(a):
    """MOND interpolation function (simple form)."""
    x = a / a_0_mond
    return x / np.sqrt(1 + x**2)  # Simple interpolation

def mond_boost(a):
    """MOND boost factor: actual acceleration / Newtonian acceleration."""
    # In deep MOND: a_eff = sqrt(a_N * a_0)
    # Boost = a_eff / a_N = sqrt(a_0 / a_N) for a_N << a_0
    return 1.0 / mond_mu(a)

def synchronism_coherence(rho_ratio):
    """C(ρ) = Ω_m + (1 - Ω_m) × (ρ/ρ_t)^(1/φ) / [1 + (ρ/ρ_t)^(1/φ)]"""
    x = rho_ratio ** (1/phi)
    return Omega_m + (1 - Omega_m) * x / (1 + x)

def synchronism_G_eff_ratio(rho_ratio):
    """G_eff/G = 1/C(ρ)"""
    return 1.0 / synchronism_coherence(rho_ratio)

# Convert to M_dyn/M_bary predictions
def mond_mass_discrepancy(a):
    """M_dyn/M_bary from MOND at acceleration a."""
    return mond_boost(a)

def synchronism_mass_discrepancy(rho_ratio):
    """M_dyn/M_bary from Synchronism at density ratio ρ/ρ_t."""
    return synchronism_G_eff_ratio(rho_ratio)


# =============================================================================
# Discriminating Test A: Same Acceleration, Different Density
# =============================================================================

def test_A_same_acceleration_different_density():
    """
    Compare TDGs vs isolated dwarfs at same rotation velocity.

    TDGs: Low density (tidal stream) + Low acceleration → BOTH predict enhancement
    Isolated dwarfs: Higher density + Low acceleration → Different predictions

    At same V_rot (same a):
    - MOND: Same mass discrepancy (depends only on a)
    - Synchronism: TDGs have MORE discrepancy (lower ρ)
    """
    print("\n" + "="*70)
    print("TEST A: Same Acceleration, Different Density")
    print("="*70)

    # Typical parameters for V_rot = 30 km/s systems
    V_rot = 30  # km/s
    R = 3  # kpc

    # Centripetal acceleration
    a_cent = V_rot**2 / (R * 3.086e16)  # Convert R to meters, get m/s^2
    print(f"\nV_rot = {V_rot} km/s, R = {R} kpc")
    print(f"Centripetal acceleration: a = {a_cent:.2e} m/s² (a/a_0 = {a_cent/a_0_mond:.2f})")

    # MOND prediction (same for both)
    mond_pred = mond_mass_discrepancy(a_cent)
    print(f"\nMOND prediction: M_dyn/M_bary = {mond_pred:.2f} (acceleration-only)")

    # Synchronism predictions for different densities
    print("\nSynchronism predictions (density-dependent):")

    scenarios = [
        ("TDG in tidal stream", 0.1),      # Low density
        ("Isolated dwarf", 0.5),            # Medium density
        ("Dwarf in cluster", 2.0),          # Higher density
        ("Dwarf near cluster core", 10.0),  # High density
    ]

    print(f"{'Scenario':<30} {'ρ/ρ_t':>8} {'M_dyn/M_bary':>15}")
    print("-" * 55)

    for scenario, rho_ratio in scenarios:
        sync_pred = synchronism_mass_discrepancy(rho_ratio)
        print(f"{scenario:<30} {rho_ratio:>8.1f} {sync_pred:>15.2f}")

    print("\n" + "-"*70)
    print("DISCRIMINATING PREDICTION:")
    print("  At same V_rot:")
    print("  - MOND: ALL should have M_dyn/M_bary ~ {:.2f}".format(mond_pred))
    print("  - Synchronism: TDGs > isolated dwarfs > cluster dwarfs")
    print("-"*70)

    return True


# =============================================================================
# Discriminating Test B: BTFR Scatter Correlations
# =============================================================================

def test_B_btfr_scatter():
    """
    BTFR scatter should correlate with:
    - MOND: Nothing (scatter is measurement error)
    - Synchronism: Environment density

    Galaxies in voids should scatter HIGH on BTFR (too fast for their mass)
    Galaxies in clusters should scatter LOW on BTFR (too slow for their mass)
    """
    print("\n" + "="*70)
    print("TEST B: BTFR Scatter Correlations")
    print("="*70)

    print("""
The Baryonic Tully-Fisher Relation: M_bary ∝ V^4

MOND Prediction:
- BTFR is exact prediction of theory
- Scatter is purely measurement error
- No correlation with environment

Synchronism Prediction:
- BTFR holds on average
- But scatter correlates with environment density:
  * Void galaxies: V_observed > V_predicted (G_eff enhanced)
  * Cluster galaxies: V_observed ≈ V_predicted (G_eff ≈ G)
    """)

    # Simulate BTFR with Synchronism environment effect
    np.random.seed(42)
    N = 100

    # True baryonic masses (log scale)
    log_M_bary = np.random.uniform(8, 11, N)  # 10^8 to 10^11 M_sun

    # Environment densities (ρ/ρ_t)
    rho_ratios = np.random.lognormal(0, 1, N)  # Log-normal distribution
    rho_ratios = np.clip(rho_ratios, 0.1, 100)  # Clip to reasonable range

    # True V_flat from MOND-like relation (V^4 ∝ M)
    V_true = 10**(log_M_bary / 4) * 10  # Arbitrary scaling

    # MOND: No environment dependence
    V_mond = V_true * (1 + 0.05 * np.random.randn(N))  # Just measurement scatter

    # Synchronism: Environment-dependent G_eff
    G_eff_ratios = synchronism_G_eff_ratio(rho_ratios)
    # V_flat ∝ sqrt(G_eff * M) → V_sync = V_true * sqrt(G_eff/G)
    V_sync = V_true * np.sqrt(G_eff_ratios) * (1 + 0.05 * np.random.randn(N))

    # Calculate residuals from mean BTFR
    log_V_mond = np.log10(V_mond)
    log_V_sync = np.log10(V_sync)

    # Fit BTFR and get residuals
    from numpy.polynomial import polynomial as P

    # MOND residuals
    p_mond = np.polyfit(log_M_bary, log_V_mond, 1)
    resid_mond = log_V_mond - np.polyval(p_mond, log_M_bary)

    # Synchronism residuals
    p_sync = np.polyfit(log_M_bary, log_V_sync, 1)
    resid_sync = log_V_sync - np.polyval(p_sync, log_M_bary)

    # Correlation with environment
    log_rho = np.log10(rho_ratios)
    corr_mond = np.corrcoef(resid_mond, log_rho)[0,1]
    corr_sync = np.corrcoef(resid_sync, log_rho)[0,1]

    print(f"\nBTFR residual correlation with log(ρ/ρ_t):")
    print(f"  MOND:        r = {corr_mond:+.3f} (no correlation expected)")
    print(f"  Synchronism: r = {corr_sync:+.3f} (NEGATIVE correlation expected)")

    print("""
DISCRIMINATING PREDICTION:
- If BTFR scatter correlates with environment density → Synchronism
- If BTFR scatter is random → MOND (or measurement error)

Observable test:
1. Measure BTFR for SPARC galaxies
2. Estimate environment density for each (from galaxy catalogs)
3. Correlate BTFR residuals with environment
    """)

    return corr_mond, corr_sync


# =============================================================================
# Discriminating Test C: External Field Effect
# =============================================================================

def test_C_external_field_effect():
    """
    The External Field Effect (EFE) is a MOND-specific prediction.

    In MOND: Internal dynamics depend on external gravitational field
    - Dwarf in cluster feels cluster's gravity
    - This can suppress MOND boost

    In Synchronism: No EFE (only local density matters)
    """
    print("\n" + "="*70)
    print("TEST C: External Field Effect (EFE)")
    print("="*70)

    print("""
MOND's External Field Effect:
- Internal dynamics of a system depend on external gravitational field
- Dwarf in cluster outskirts: EFE suppresses MOND boost
- Same dwarf in isolation: Full MOND boost

Synchronism:
- NO External Field Effect
- Only local density matters
- Dwarf dynamics same regardless of external field
    """)

    # Example: Dwarf at R = 500 kpc from cluster center
    R_cluster = 500  # kpc from cluster center
    M_cluster = 1e14  # M_sun

    # External field from cluster
    g_ext = G * M_cluster / R_cluster**2  # (km/s)^2 / kpc
    g_ext_si = g_ext * 1e6 / 3.086e16  # Convert to m/s^2

    print(f"\nExample: Dwarf at {R_cluster} kpc from cluster (M = 10^14 M_sun)")
    print(f"External acceleration: g_ext = {g_ext_si:.2e} m/s^2")
    print(f"Ratio to a_0: g_ext/a_0 = {g_ext_si/a_0_mond:.2f}")

    # Internal acceleration of dwarf
    V_int = 20  # km/s internal velocity
    R_int = 2   # kpc internal radius
    g_int = V_int**2 / R_int  # (km/s)^2 / kpc
    g_int_si = g_int * 1e6 / 3.086e16

    print(f"\nDwarf internal acceleration: g_int = {g_int_si:.2e} m/s^2")
    print(f"Ratio to a_0: g_int/a_0 = {g_int_si/a_0_mond:.2f}")

    print(f"\nMOND prediction:")
    if g_ext_si > a_0_mond:
        print(f"  EFE regime: g_ext > a_0")
        print(f"  Internal dynamics become NEWTONIAN")
        print(f"  M_dyn/M_bary ≈ 1.0")
    elif g_ext_si > g_int_si:
        print(f"  EFE intermediate: g_int < g_ext < a_0")
        print(f"  MOND boost suppressed by EFE")
        mond_boost_suppressed = (g_int_si / g_ext_si)**0.5
        print(f"  M_dyn/M_bary ≈ {1/mond_boost_suppressed:.2f}")
    else:
        print(f"  Standard MOND: g_ext < g_int < a_0")
        mond_boost_val = mond_mass_discrepancy(g_int_si)
        print(f"  M_dyn/M_bary ≈ {mond_boost_val:.2f}")

    print(f"\nSynchronism prediction:")
    # Density at 500 kpc from cluster (NFW-ish profile)
    rho_ratio_cluster = 0.5  # Moderate overdensity in cluster outskirts
    sync_pred = synchronism_mass_discrepancy(rho_ratio_cluster)
    print(f"  ρ/ρ_t ≈ {rho_ratio_cluster} (cluster outskirts)")
    print(f"  M_dyn/M_bary = {sync_pred:.2f}")

    print("""
DISCRIMINATING PREDICTION:
- In MOND: Dwarfs in cluster outskirts show REDUCED mass discrepancy (EFE)
- In Synchronism: Dwarfs in cluster outskirts show INCREASED mass discrepancy
  (lower density than field due to cluster tidal effects)

Observable test:
1. Compare dwarf spheroidals at similar distances from cluster center
2. MOND predicts: Cluster dwarfs have LESS mass discrepancy than field dwarfs
3. Synchronism predicts: Depends on local density, not external field
    """)

    return True


# =============================================================================
# Discriminating Test D: Density vs Acceleration Plane
# =============================================================================

def test_D_density_acceleration_plane():
    """
    Map mass discrepancy in the (ρ, a) plane.

    MOND: Contours are horizontal (depend only on a)
    Synchronism: Contours are vertical (depend only on ρ)
    """
    print("\n" + "="*70)
    print("TEST D: Density vs Acceleration Plane")
    print("="*70)

    # Create 2D grid
    log_rho = np.linspace(-2, 2, 100)  # log10(ρ/ρ_t)
    log_a = np.linspace(-12, -8, 100)  # log10(a) in m/s^2

    RHO, A = np.meshgrid(10**log_rho, 10**log_a)

    # MOND prediction (depends only on a)
    Z_mond = mond_mass_discrepancy(A)

    # Synchronism prediction (depends only on ρ)
    Z_sync = synchronism_mass_discrepancy(RHO)

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # MOND
    ax1 = axes[0]
    c1 = ax1.contourf(log_rho, log_a, np.log10(Z_mond), levels=20, cmap='viridis')
    ax1.contour(log_rho, log_a, np.log10(Z_mond), levels=[0, 0.3, 0.5], colors='white', linewidths=1)
    ax1.axhline(np.log10(a_0_mond), color='red', linestyle='--', label='a = a₀')
    ax1.set_xlabel('log₁₀(ρ/ρ_t)', fontsize=12)
    ax1.set_ylabel('log₁₀(a) [m/s²]', fontsize=12)
    ax1.set_title('MOND: M_dyn/M_bary\n(Horizontal contours = acceleration-only)', fontsize=11)
    plt.colorbar(c1, ax=ax1, label='log₁₀(M_dyn/M_bary)')
    ax1.legend()

    # Synchronism
    ax2 = axes[1]
    c2 = ax2.contourf(log_rho, log_a, np.log10(Z_sync), levels=20, cmap='viridis')
    ax2.contour(log_rho, log_a, np.log10(Z_sync), levels=[0, 0.3, 0.5], colors='white', linewidths=1)
    ax2.axvline(0, color='red', linestyle='--', label='ρ = ρ_t')
    ax2.set_xlabel('log₁₀(ρ/ρ_t)', fontsize=12)
    ax2.set_ylabel('log₁₀(a) [m/s²]', fontsize=12)
    ax2.set_title('Synchronism: M_dyn/M_bary\n(Vertical contours = density-only)', fontsize=11)
    plt.colorbar(c2, ax=ax2, label='log₁₀(M_dyn/M_bary)')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session183_sync_vs_mond.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    print("\nFigure saved: session183_sync_vs_mond.png")
    print("""
The figure shows M_dyn/M_bary in the (ρ, a) plane:

MOND (left panel):
- Contours are HORIZONTAL
- Mass discrepancy depends ONLY on acceleration
- Same discrepancy at all densities for fixed a

Synchronism (right panel):
- Contours are VERTICAL
- Mass discrepancy depends ONLY on density
- Same discrepancy at all accelerations for fixed ρ

DISCRIMINATING OBSERVATION:
Find objects at SAME density but DIFFERENT acceleration:
- If same M_dyn/M_bary → Synchronism
- If different M_dyn/M_bary → MOND

Find objects at SAME acceleration but DIFFERENT density:
- If same M_dyn/M_bary → MOND
- If different M_dyn/M_bary → Synchronism
    """)

    return True


# =============================================================================
# Summary
# =============================================================================

def summarize_discrimination():
    """Print summary of discriminating tests."""
    print("\n" + "="*70)
    print("SESSION #183: SYNCHRONISM VS MOND DISCRIMINATION SUMMARY")
    print("="*70)

    print("""
FUNDAMENTAL DIFFERENCE:
- MOND: Mass discrepancy ∝ f(acceleration)
- Synchronism: Mass discrepancy ∝ g(density)

DISCRIMINATING TESTS:

┌─────────────────────────────────────────────────────────────────────┐
│ Test A: Same Acceleration, Different Density                       │
├─────────────────────────────────────────────────────────────────────┤
│ Compare TDGs vs isolated dwarfs at same V_rot                      │
│ MOND: Same M_dyn/M_bary                                             │
│ Synchronism: TDGs have MORE discrepancy (lower ρ)                   │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│ Test B: BTFR Scatter Correlations                                   │
├─────────────────────────────────────────────────────────────────────┤
│ Correlate BTFR residuals with environment density                   │
│ MOND: No correlation (scatter is random/measurement)                │
│ Synchronism: NEGATIVE correlation (voids scatter high)              │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│ Test C: External Field Effect                                       │
├─────────────────────────────────────────────────────────────────────┤
│ Compare dwarfs in cluster outskirts vs isolation                    │
│ MOND: Cluster dwarfs have LESS discrepancy (EFE suppression)        │
│ Synchronism: No EFE; depends only on local density                  │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│ Test D: Density-Acceleration Plane                                  │
├─────────────────────────────────────────────────────────────────────┤
│ Map M_dyn/M_bary vs (ρ, a)                                          │
│ MOND: Horizontal contours (a-dependence only)                       │
│ Synchronism: Vertical contours (ρ-dependence only)                  │
└─────────────────────────────────────────────────────────────────────┘

EXISTING DATA:
1. Chae et al. 2020: Detected EFE signal → Supports MOND
   BUT: Could also be environment-density correlation

2. TDG observations: Higher M_dyn/M_bary than isolated dwarfs
   → Could support Synchronism (lower density in tidal streams)

3. Fornax cluster dwarfs: Disturbed morphologies
   → MOND explains via EFE
   → Synchronism explains via low ρ in tidal environment

NEEDED OBSERVATIONS:
1. TDG vs isolated dwarf comparison at matched V_rot
2. BTFR residual correlation with environment density
3. Cluster dwarf velocity dispersions vs field dwarfs
4. Direct (ρ, a) mapping for sample of dwarfs
    """)


def main():
    """Run all discriminating tests."""
    print("="*70)
    print("Session #183: Synchronism vs MOND Discrimination")
    print("="*70)

    test_A_same_acceleration_different_density()
    test_B_btfr_scatter()
    test_C_external_field_effect()
    test_D_density_acceleration_plane()
    summarize_discrimination()

if __name__ == "__main__":
    main()
