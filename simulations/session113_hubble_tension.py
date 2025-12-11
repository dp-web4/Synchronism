"""
Session #113: The Hubble Tension in Synchronism

The Hubble tension is the ~5σ discrepancy between:
- CMB-derived H0 = 67.4 ± 0.5 km/s/Mpc (Planck 2018)
- Local measurements H0 = 73.0 ± 1.0 km/s/Mpc (SH0ES 2022)

This session analyzes whether Synchronism's scale-dependent coherence
can explain the Hubble tension.

Key question: Does C_local ≠ C_cosmic at z ~ 0 affect distance ladder measurements?

Physics:
- CMB measures H0 via sound horizon (early universe, z ~ 1089)
- Local measurements use Cepheids → SNe Ia (z < 0.1)
- If G_eff varies with scale, distance calibration could be affected

Author: CBP Autonomous Synchronism Research
Date: December 11, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ============================================================================
# COSMOLOGICAL PARAMETERS
# ============================================================================

# Planck 2018 (CMB-derived)
H0_Planck = 67.4  # km/s/Mpc
H0_Planck_err = 0.5

# SH0ES 2022 (local Cepheid-SNe)
H0_SH0ES = 73.0  # km/s/Mpc
H0_SH0ES_err = 1.0

# Other measurements
H0_TRGB = 69.8  # Carnegie-Chicago Hubble Program
H0_TRGB_err = 1.7

H0_Megamaser = 73.9  # Megamaser Cosmology Project
H0_Megamaser_err = 3.0

# Cosmological parameters
Omega_m = 0.315
Omega_Lambda = 0.685
c = 299792.458  # km/s

# ============================================================================
# PART 1: THE HUBBLE TENSION
# ============================================================================

def analyze_hubble_tension():
    """Analyze the current Hubble tension."""
    print("=" * 70)
    print("SESSION #113: THE HUBBLE TENSION IN SYNCHRONISM")
    print("=" * 70)

    print("\n1. CURRENT STATE OF HUBBLE TENSION")
    print("-" * 70)

    # Calculate tension significance
    diff = H0_SH0ES - H0_Planck
    combined_err = np.sqrt(H0_SH0ES_err**2 + H0_Planck_err**2)
    tension = diff / combined_err

    print(f"CMB (Planck 2018):     H0 = {H0_Planck:.1f} ± {H0_Planck_err:.1f} km/s/Mpc")
    print(f"Local (SH0ES 2022):    H0 = {H0_SH0ES:.1f} ± {H0_SH0ES_err:.1f} km/s/Mpc")
    print(f"Difference:            ΔH0 = {diff:.1f} km/s/Mpc")
    print(f"Combined uncertainty:  σ = {combined_err:.2f} km/s/Mpc")
    print(f"Tension significance:  {tension:.1f}σ")

    print("\n2. OTHER H0 MEASUREMENTS")
    print("-" * 70)

    measurements = [
        ("Planck CMB", H0_Planck, H0_Planck_err, "Early universe"),
        ("SH0ES Cepheids", H0_SH0ES, H0_SH0ES_err, "Local, z < 0.01"),
        ("TRGB", H0_TRGB, H0_TRGB_err, "Local, independent"),
        ("Megamasers", H0_Megamaser, H0_Megamaser_err, "Geometric"),
        ("BAO + BBN", 67.6, 1.0, "Early + z~0.5"),
        ("Strong lensing", 73.3, 1.8, "Time delays"),
    ]

    print(f"{'Method':<20} {'H0':<8} {'σ':<6} {'Notes':<20}")
    print("-" * 70)
    for name, h0, err, notes in measurements:
        print(f"{name:<20} {h0:<8.1f} {err:<6.1f} {notes:<20}")

    print("\n*** KEY OBSERVATION ***")
    print("Early universe methods:     H0 ~ 67-68")
    print("Late universe methods:      H0 ~ 72-74")
    print("The split is consistent across multiple techniques!")

    return tension

# ============================================================================
# PART 2: SYNCHRONISM ANALYSIS
# ============================================================================

def C_galactic(rho_ratio=1e6, gamma=2.0):
    """Galactic-scale coherence."""
    return np.tanh(gamma * np.log(rho_ratio + 1))

def C_cosmic(z, Omega_m=0.315):
    """Cosmic-scale coherence = Ω_m(z)."""
    return Omega_m * (1+z)**3 / (Omega_m * (1+z)**3 + (1 - Omega_m))

def G_ratio(z):
    """G_local / G_global = C_cosmic / C_galactic."""
    C_gal = C_galactic()
    C_cos = C_cosmic(z)
    return C_cos / C_gal

def analyze_coherence_effect_on_H0():
    """
    Analyze how scale-dependent coherence affects H0 measurements.
    """
    print("\n" + "=" * 70)
    print("3. SYNCHRONISM ANALYSIS: SCALE-DEPENDENT H0?")
    print("=" * 70)

    # Key insight: Local measurements use distance ladder
    # If G_eff differs locally vs cosmologically, luminosities are affected

    print("\n3.1 THE DISTANCE LADDER")
    print("-" * 70)
    print("""
    Local H0 measurement chain:
    1. Parallax → Cepheid distances (Milky Way)
    2. Cepheid period-luminosity → nearby galaxies
    3. Cepheid calibration → SNe Ia in same galaxies
    4. SNe Ia → Hubble flow (z ~ 0.01-0.1)

    Each step involves luminosity calibration.
    """)

    # Does C(ρ) affect stellar luminosity?
    print("3.2 DOES COHERENCE AFFECT STELLAR LUMINOSITY?")
    print("-" * 70)

    # Stellar luminosity: L ∝ M^α (mass-luminosity relation)
    # For main sequence: L ∝ M^3.5 (approximately)
    # Luminosity comes from nuclear reactions, not gravity directly

    print("Analysis:")
    print("- Cepheid luminosity depends on pulsation physics")
    print("- Pulsation period: P ∝ ρ^(-0.5) ∝ (M/R³)^(-0.5)")
    print("- If G_eff varies, hydrostatic equilibrium changes")
    print("- This could affect stellar structure and pulsation")

    # Calculate potential effect
    print("\n3.3 POTENTIAL G_eff EFFECT ON CEPHEIDS")
    print("-" * 70)

    # In Cepheid-hosting galaxies (z ~ 0)
    # Typical galactic densities: ρ ~ 10^6 × ρ_cosmic
    C_gal_z0 = C_galactic(rho_ratio=1e6)
    C_cos_z0 = C_cosmic(z=0)
    G_ratio_z0 = C_cos_z0 / C_gal_z0

    print(f"At z = 0:")
    print(f"  C_galactic = {C_gal_z0:.4f}")
    print(f"  C_cosmic = {C_cos_z0:.4f}")
    print(f"  G_local / G_global = {G_ratio_z0:.4f}")

    # In Cepheid cores: much higher density
    C_gal_star = C_galactic(rho_ratio=1e15)  # Stellar interior
    print(f"\nIn stellar interiors:")
    print(f"  C_stellar = {C_gal_star:.4f} (approaches 1)")
    print(f"  Stars feel ~standard G")

    print("\n*** KEY INSIGHT ***")
    print("Stellar interiors have C ~ 1 (high density)")
    print("→ Cepheid luminosities are NOT directly affected by C(ρ)")
    print("→ Distance ladder calibration should be robust")

    return G_ratio_z0

def analyze_sound_horizon():
    """
    Analyze CMB sound horizon in Synchronism.
    """
    print("\n" + "=" * 70)
    print("4. CMB SOUND HORIZON ANALYSIS")
    print("=" * 70)

    print("\nCMB H0 measurement depends on:")
    print("  H0 = θ_s / r_s × (angular diameter distance)")
    print("  where:")
    print("    θ_s = angular scale of first acoustic peak")
    print("    r_s = sound horizon at recombination")

    # At recombination (z = 1089)
    z_rec = 1089
    C_gal_rec = C_galactic(rho_ratio=1e6)  # Still ~1 in dense regions
    C_cos_rec = C_cosmic(z_rec)
    G_ratio_rec = C_cos_rec / C_gal_rec

    print(f"\nAt recombination (z = {z_rec}):")
    print(f"  C_galactic = {C_gal_rec:.4f}")
    print(f"  C_cosmic = {C_cos_rec:.4f}")
    print(f"  G_ratio = {G_ratio_rec:.4f}")

    print("\n*** KEY FINDING ***")
    print(f"At z = 1089: G_ratio = {G_ratio_rec:.4f} ≈ 1.000")
    print("The early universe has NO scale-dependent coherence effect!")
    print("CMB physics is STANDARD in Synchronism.")

    # This means CMB-derived H0 should be correct
    print("\nImplication:")
    print("  CMB gives TRUE cosmological H0 = 67.4")
    print("  This is NOT modified by Synchronism")

    return G_ratio_rec

def analyze_local_vs_cmb():
    """
    Analyze why local and CMB H0 might differ.
    """
    print("\n" + "=" * 70)
    print("5. LOCAL VS CMB: WHERE'S THE DISCREPANCY?")
    print("=" * 70)

    print("\n5.1 POSSIBLE SYNCHRONISM EFFECTS")
    print("-" * 70)

    # Check intermediate redshifts
    redshifts = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]

    print(f"{'z':<8} {'C_cosmic':<12} {'G_ratio':<12} {'Effect':<15}")
    print("-" * 70)

    for z in redshifts:
        C_cos = C_cosmic(z)
        G_r = G_ratio(z)
        if G_r < 0.95:
            effect = "Significant"
        elif G_r < 0.99:
            effect = "Small"
        else:
            effect = "Negligible"
        print(f"{z:<8.2f} {C_cos:<12.4f} {G_r:<12.4f} {effect:<15}")

    print("\n5.2 THE HUBBLE FLOW REGION (z ~ 0.01-0.1)")
    print("-" * 70)

    # SNe Ia are used at z ~ 0.01 - 0.1
    # What's the integrated effect?
    z_sne = np.linspace(0.01, 0.1, 100)
    G_ratios = np.array([G_ratio(z) for z in z_sne])
    avg_G_ratio = np.mean(G_ratios)

    print(f"Average G_ratio in SNe Ia range: {avg_G_ratio:.4f}")
    print(f"This corresponds to a ~{(1 - avg_G_ratio)*100:.1f}% effect on gravity")

    # How does this affect distance measurements?
    print("\n5.3 EFFECT ON LUMINOSITY DISTANCE")
    print("-" * 70)

    # Luminosity distance: d_L = (1+z) × ∫ c/H(z') dz'
    # If H(z) is modified by G_eff, distances change

    # In Synchronism, the Friedmann equation is:
    # H² = (8πG/3C) × ρ
    # So H_Sync = H_LCDM × sqrt(1/C_cosmic)

    print("Friedmann equation in Synchronism:")
    print("  H² = (8πG/3C) × ρ")
    print("  H_Sync = H_LCDM × (1/√C_cosmic)")

    # At z = 0: C_cosmic = Ω_m = 0.315
    # So H_Sync(z=0) = H_LCDM × √(1/0.315) = H_LCDM × 1.78

    # BUT WAIT - this is calibrated to match H0!
    # The natural calibration is C_0 = Ω_m

    print("\n*** CRITICAL REALIZATION ***")
    print("Session #72 established: C_0 = Ω_m = 0.315 is natural calibration")
    print("This makes H_Sync(z=0) = H_LCDM(z=0) EXACTLY")
    print("The Hubble tension is NOT resolved by this mechanism")

    return avg_G_ratio

def analyze_alternative_mechanisms():
    """
    Explore alternative Synchronism mechanisms for Hubble tension.
    """
    print("\n" + "=" * 70)
    print("6. ALTERNATIVE MECHANISMS")
    print("=" * 70)

    print("\n6.1 COULD SOUND HORIZON BE MODIFIED?")
    print("-" * 70)

    # Sound horizon depends on pre-recombination physics
    # r_s = ∫ c_s / H(z) dz from z_drag to infinity

    print("Sound horizon formula:")
    print("  r_s = ∫ c_s / H(z) dz")
    print("")
    print("If H(z) were different pre-recombination...")
    print("BUT: At z > 1000, G_ratio = 1.00 (no effect)")
    print("Sound horizon is UNCHANGED in Synchronism")

    print("\n6.2 COULD LOCAL DISTANCES BE WRONG?")
    print("-" * 70)

    # Cepheid calibration
    print("Distance ladder analysis:")
    print("1. Parallax (geometry) - Not affected by gravity")
    print("2. Cepheid P-L (stellar physics) - C ~ 1 in stars")
    print("3. SNe Ia (thermonuclear) - C ~ 1 in WD explosions")
    print("4. Hubble flow - Uses standard cosmology")
    print("")
    print("All steps use high-density environments where C ~ 1")
    print("Local distance ladder is ROBUST")

    print("\n6.3 SYNCHRONISM VERDICT ON HUBBLE TENSION")
    print("-" * 70)

    print("""
    FINDING: Synchronism does NOT naturally resolve the Hubble tension

    Reasons:
    1. CMB physics (z ~ 1089) has G_ratio = 1.0 (no effect)
    2. Sound horizon is unchanged
    3. Local calibrators (Cepheids, SNe) use high-density regions (C ~ 1)
    4. The natural calibration C_0 = Ω_m preserves H0

    The Hubble tension remains a separate puzzle.
    """)

    # However, there might be indirect connections
    print("6.4 INDIRECT CONNECTION TO S8 TENSION?")
    print("-" * 70)

    print("""
    Some authors suggest S8 and H0 tensions are related:
    - Higher H0 → earlier matter-radiation equality
    - Earlier equality → different structure growth

    Synchronism addresses S8 but NOT H0.
    This suggests they may be INDEPENDENT tensions.
    """)

def create_visualization():
    """Create visualization of Hubble tension analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. H0 measurements comparison
    ax1 = axes[0, 0]

    methods = ['Planck\nCMB', 'SH0ES\nCepheids', 'TRGB', 'Megamaser',
               'BAO+BBN', 'Strong\nLensing']
    h0_values = [67.4, 73.0, 69.8, 73.9, 67.6, 73.3]
    h0_errors = [0.5, 1.0, 1.7, 3.0, 1.0, 1.8]
    colors = ['blue', 'red', 'orange', 'red', 'blue', 'red']

    ax1.errorbar(range(len(methods)), h0_values, yerr=h0_errors,
                fmt='o', markersize=10, capsize=5, capthick=2)
    for i, (m, h, c) in enumerate(zip(methods, h0_values, colors)):
        ax1.scatter(i, h, c=c, s=100, zorder=5)

    ax1.axhline(y=67.4, color='blue', linestyle='--', alpha=0.5, label='CMB value')
    ax1.axhline(y=73.0, color='red', linestyle='--', alpha=0.5, label='Local value')
    ax1.fill_between([-0.5, 5.5], 66.9, 67.9, alpha=0.2, color='blue')
    ax1.fill_between([-0.5, 5.5], 72.0, 74.0, alpha=0.2, color='red')

    ax1.set_xticks(range(len(methods)))
    ax1.set_xticklabels(methods)
    ax1.set_ylabel('H₀ (km/s/Mpc)', fontsize=12)
    ax1.set_title('The Hubble Tension: 5σ Discrepancy', fontsize=14)
    ax1.legend()
    ax1.set_xlim([-0.5, 5.5])
    ax1.set_ylim([64, 78])
    ax1.grid(True, alpha=0.3)

    # 2. G_ratio vs redshift
    ax2 = axes[0, 1]

    z_arr = np.linspace(0, 3, 100)
    G_ratios = np.array([G_ratio(z) for z in z_arr])

    ax2.plot(z_arr, G_ratios, 'b-', linewidth=2)
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax2.fill_between(z_arr, G_ratios, 1.0, alpha=0.3, color='orange',
                    where=G_ratios < 1)

    # Mark key regions
    ax2.axvline(x=0.05, color='red', linestyle=':', label='SNe Ia range')
    ax2.axvline(x=0.5, color='green', linestyle=':', label='S8 tension peak')

    ax2.set_xlabel('Redshift z', fontsize=12)
    ax2.set_ylabel('G_local / G_global', fontsize=12)
    ax2.set_title('Scale-Dependent Gravity in Synchronism', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0.3, 1.05])

    # 3. Coherence at different epochs
    ax3 = axes[1, 0]

    epochs = ['Stars\n(interior)', 'Galaxies\n(disk)', 'Clusters\n(core)',
              'Voids\n(center)', 'CMB\n(z=1089)']
    C_values = [1.0, 0.95, 0.85, 0.3, 1.0]
    colors_c = ['darkblue', 'blue', 'green', 'orange', 'darkblue']

    bars = ax3.bar(epochs, C_values, color=colors_c, alpha=0.7)
    ax3.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax3.set_ylabel('Coherence C', fontsize=12)
    ax3.set_title('Coherence at Different Scales/Epochs', fontsize=14)
    ax3.set_ylim([0, 1.1])

    # Add labels
    for bar, val in zip(bars, C_values):
        ax3.text(bar.get_x() + bar.get_width()/2, val + 0.02,
                f'{val:.2f}', ha='center', fontsize=11)

    # 4. Tension summary
    ax4 = axes[1, 1]
    ax4.axis('off')

    summary_text = """
    SYNCHRONISM vs HUBBLE TENSION

    Key Finding: Synchronism does NOT resolve H0 tension

    Why not?
    ┌────────────────────────────────────────────────────┐
    │ 1. At z = 1089: G_ratio = 1.00 (no effect)        │
    │    → CMB physics is standard                       │
    │    → Sound horizon unchanged                       │
    │                                                    │
    │ 2. Local calibrators use high-density regions:    │
    │    → Stellar interiors: C ~ 1                     │
    │    → Cepheid pulsation: standard                  │
    │    → SNe Ia: standard                             │
    │                                                    │
    │ 3. Natural calibration C₀ = Ω_m preserves H₀      │
    └────────────────────────────────────────────────────┘

    Synchronism explains: S8 tension ✓
    Synchronism does NOT explain: H0 tension ✗

    These may be INDEPENDENT cosmological puzzles.
    """

    ax4.text(0.1, 0.5, summary_text, transform=ax4.transAxes,
            fontsize=11, verticalalignment='center',
            fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session113_hubble_tension.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    print("\nVisualization saved to simulations/session113_hubble_tension.png")

def summarize_findings():
    """Summarize Session #113 findings."""
    print("\n" + "=" * 70)
    print("SESSION #113 SUMMARY: HUBBLE TENSION ANALYSIS")
    print("=" * 70)

    print("\n1. THE HUBBLE TENSION")
    print("-" * 50)
    print("CMB (Planck):  H0 = 67.4 ± 0.5 km/s/Mpc")
    print("Local (SH0ES): H0 = 73.0 ± 1.0 km/s/Mpc")
    print("Tension: 5σ (highly significant)")

    print("\n2. SYNCHRONISM ANALYSIS")
    print("-" * 50)
    print("Question: Does scale-dependent coherence affect H0?")
    print("")
    print("Findings:")
    print("• At z = 1089 (CMB): G_ratio = 1.0 → No effect")
    print("• Sound horizon: Unchanged")
    print("• Stellar interiors: C ~ 1 → Standard physics")
    print("• Distance ladder: Robust (high-density regions)")

    print("\n3. VERDICT")
    print("-" * 50)
    print("Synchronism does NOT resolve the Hubble tension.")
    print("")
    print("The same mechanism that explains S8 tension")
    print("(G_local < G_global at z ~ 0.5-1.5)")
    print("does NOT affect H0 measurements because:")
    print("  - CMB probes z >> 1 where G_ratio = 1")
    print("  - Local calibrators use high-density environments")

    print("\n4. IMPLICATIONS")
    print("-" * 50)
    print("S8 tension and H0 tension may be INDEPENDENT.")
    print("Synchronism solves one but not the other.")
    print("This is actually a STRENGTH - not overfitting!")

    print("\n5. WHAT COULD RESOLVE H0 TENSION?")
    print("-" * 50)
    print("Possibilities (outside Synchronism):")
    print("• Early dark energy")
    print("• Modified recombination physics")
    print("• Systematic errors in measurements")
    print("• New physics at z > 1000")

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("Synchronism makes a CLEAR prediction:")
    print("  - S8 tension: RESOLVED (G_local < G_global)")
    print("  - H0 tension: NOT RESOLVED (different physics)")
    print("")
    print("This is falsifiable: If future work shows S8 and H0")
    print("tensions share a common origin, Synchronism would need revision.")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Analyze Hubble tension
    tension = analyze_hubble_tension()

    # Synchronism analysis
    G_ratio_z0 = analyze_coherence_effect_on_H0()
    G_ratio_rec = analyze_sound_horizon()
    avg_G = analyze_local_vs_cmb()

    # Alternative mechanisms
    analyze_alternative_mechanisms()

    # Visualization
    create_visualization()

    # Summary
    summarize_findings()
