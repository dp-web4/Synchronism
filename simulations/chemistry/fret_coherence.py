#!/usr/bin/env python3
"""
FRET (Förster Resonance Energy Transfer) Coherence Analysis
Chemistry Session #233

Explores how γ ~ 1 manifests in resonance energy transfer:
1. Förster radius R₀ - characteristic distance
2. FRET efficiency - E = R₀⁶/(R₀⁶ + r⁶)
3. Spectral overlap integral - donor-acceptor matching
4. Orientation factor κ² - dipole alignment
5. Quantum yield - donor emission efficiency
6. Distance measurements - molecular ruler
7. FRET pairs - optimal combinations

Key insight: FRET IS resonance energy transfer, with R/R₀ = 1
as the γ ~ 1 transition point between donor and acceptor dominance.
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("FRET (FÖRSTER RESONANCE ENERGY TRANSFER) COHERENCE ANALYSIS")
print("Chemistry Session #233: γ ~ 1 in Dipole-Dipole Energy Transfer")
print("=" * 70)

# =============================================================================
# 1. FÖRSTER RADIUS (R₀)
# =============================================================================
print("\n" + "=" * 70)
print("1. FÖRSTER RADIUS: THE CHARACTERISTIC DISTANCE")
print("=" * 70)

# R₀ is THE γ ~ 1 reference distance for FRET
# At r = R₀: E = 50% energy transfer

print("  Förster radius R₀: distance at 50% energy transfer")
print("  R₀⁶ = (9000 × ln(10) × κ² × Φ_D × J) / (128π⁵ × n⁴ × N_A)")
print("  Typical R₀ = 2-9 nm for common FRET pairs")

# Common FRET pairs and their R₀ values
fret_pairs = {
    'CFP-YFP': {'R0_nm': 4.9, 'donor': 'CFP', 'acceptor': 'YFP'},
    'BFP-GFP': {'R0_nm': 4.1, 'donor': 'BFP', 'acceptor': 'GFP'},
    'Cy3-Cy5': {'R0_nm': 5.4, 'donor': 'Cy3', 'acceptor': 'Cy5'},
    'FAM-TAMRA': {'R0_nm': 5.0, 'donor': 'FAM', 'acceptor': 'TAMRA'},
    'Alexa488-Alexa555': {'R0_nm': 6.4, 'donor': 'Alexa488', 'acceptor': 'Alexa555'},
    'EGFP-mCherry': {'R0_nm': 5.2, 'donor': 'EGFP', 'acceptor': 'mCherry'},
    'Trp-Dansyl': {'R0_nm': 2.4, 'donor': 'Tryptophan', 'acceptor': 'Dansyl'},
    'Fluorescein-Rhodamine': {'R0_nm': 5.5, 'donor': 'Fluorescein', 'acceptor': 'Rhodamine'},
}

print(f"\n  FRET Pair             | R₀ (nm) | Donor → Acceptor")
print("  " + "-" * 60)
for pair, data in fret_pairs.items():
    print(f"  {pair:22s} | {data['R0_nm']:.1f}     | {data['donor']} → {data['acceptor']}")

R0_mean = np.mean([d['R0_nm'] for d in fret_pairs.values()])
R0_std = np.std([d['R0_nm'] for d in fret_pairs.values()])
print(f"\n  Mean R₀: {R0_mean:.2f} ± {R0_std:.2f} nm")
print(f"  R₀ IS the γ ~ 1 reference distance for FRET!")

# =============================================================================
# 2. FRET EFFICIENCY
# =============================================================================
print("\n" + "=" * 70)
print("2. FRET EFFICIENCY: E = R₀⁶/(R₀⁶ + r⁶)")
print("=" * 70)

# FRET efficiency equation
# E = 1 / [1 + (r/R₀)⁶]
# At r = R₀: E = 1/2 = 50% (γ ~ 1!)

def fret_efficiency(r, R0):
    """Calculate FRET efficiency from distance ratio"""
    return 1 / (1 + (r/R0)**6)

print("  E = R₀⁶/(R₀⁶ + r⁶) = 1/[1 + (r/R₀)⁶]")
print("\n  Distance r/R₀ | FRET Efficiency E | Comment")
print("  " + "-" * 55)

r_ratios = [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]
for ratio in r_ratios:
    E = fret_efficiency(ratio, 1.0)
    comment = "50% - γ ~ 1!" if ratio == 1.0 else ""
    print(f"  {ratio:8.2f}      | {E:.4f}            | {comment}")

print(f"\n  At r = R₀: E = 0.500 EXACTLY (γ ~ 1!)")
print(f"  This IS the donor-acceptor transition point")
print(f"  r < R₀: donor-dominated (E > 0.5)")
print(f"  r > R₀: acceptor-dominated (E < 0.5)")

# =============================================================================
# 3. SPECTRAL OVERLAP INTEGRAL
# =============================================================================
print("\n" + "=" * 70)
print("3. SPECTRAL OVERLAP INTEGRAL (J)")
print("=" * 70)

# J = ∫ F_D(λ) × ε_A(λ) × λ⁴ dλ
# Normalized donor emission × acceptor absorption

print("  J = ∫ F_D(λ) × ε_A(λ) × λ⁴ dλ")
print("  where:")
print("    F_D(λ) = normalized donor emission spectrum")
print("    ε_A(λ) = acceptor molar extinction coefficient")

# Example overlap integrals
overlap_data = {
    'CFP-YFP': {'J': 1.5e15, 'overlap': 'Good'},
    'Cy3-Cy5': {'J': 7.6e15, 'overlap': 'Excellent'},
    'FITC-Rhodamine': {'J': 4.2e15, 'overlap': 'Good'},
    'Trp-Dansyl': {'J': 3.9e14, 'overlap': 'Moderate'},
    'BFP-GFP': {'J': 8.2e14, 'overlap': 'Moderate'},
}

print(f"\n  Pair              | J (M⁻¹cm⁻¹nm⁴) | Overlap Quality")
print("  " + "-" * 55)
for pair, data in overlap_data.items():
    print(f"  {pair:17s} | {data['J']:.2e}       | {data['overlap']}")

print(f"\n  Optimal spectral overlap:")
print(f"    Donor emission peak ≈ Acceptor absorption peak")
print(f"    This IS resonance matching (γ ~ 1!)")
print(f"    Maximum J when spectra coincide")

# =============================================================================
# 4. ORIENTATION FACTOR (κ²)
# =============================================================================
print("\n" + "=" * 70)
print("4. ORIENTATION FACTOR (κ²)")
print("=" * 70)

# κ² describes relative dipole orientations
# κ² = (cos θ_T - 3 cos θ_D cos θ_A)²

print("  κ² = (cos θ_T - 3 cos θ_D cos θ_A)²")
print("  where θ_T, θ_D, θ_A describe dipole orientations")

kappa_values = {
    'Parallel (aligned)': 4.0,
    'Perpendicular': 0.0,
    'Random static': 0.476,
    'Random dynamic': 2/3,  # 0.667 - THE isotropic average (γ ~ 1!)
    'Head-to-tail': 1.0,
}

print(f"\n  Orientation        | κ²    | Comment")
print("  " + "-" * 50)
for orient, k2 in kappa_values.items():
    comment = "isotropic average (γ ~ 1!)" if orient == 'Random dynamic' else ""
    print(f"  {orient:20s} | {k2:.3f} | {comment}")

print(f"\n  κ² = 2/3 (0.667) for randomly oriented, rapidly tumbling molecules")
print(f"  This IS the isotropic average - γ ~ 1 reference!")
print(f"  Most solution FRET experiments assume κ² = 2/3")

# =============================================================================
# 5. DONOR QUANTUM YIELD
# =============================================================================
print("\n" + "=" * 70)
print("5. DONOR QUANTUM YIELD (Φ_D)")
print("=" * 70)

# Φ_D = photons emitted / photons absorbed
# Ideal: Φ_D = 1.0 (γ ~ 1!)

print("  Φ_D = photons emitted / photons absorbed (in absence of acceptor)")
print("  Ideal: Φ_D = 1.0 (γ ~ 1)")

donor_qy = {
    'Fluorescein': 0.93,
    'Rhodamine 6G': 0.95,
    'CFP': 0.40,
    'GFP': 0.60,
    'YFP': 0.61,
    'mCherry': 0.22,
    'Cy3': 0.15,
    'Cy5': 0.28,
    'Alexa 488': 0.92,
    'Alexa 555': 0.10,
}

print(f"\n  Fluorophore   | Φ_D   | Φ_D/1.0")
print("  " + "-" * 40)
for fluor, qy in sorted(donor_qy.items(), key=lambda x: -x[1]):
    print(f"  {fluor:13s} | {qy:.2f}  | {qy:.3f}")

avg_qy = np.mean(list(donor_qy.values()))
print(f"\n  Average Φ_D: {avg_qy:.2f}")
print(f"  Fluorescein, Rhodamine 6G, Alexa 488 approach Φ_D ~ 1")
print(f"  High Φ_D essential for sensitive FRET detection")

# =============================================================================
# 6. DISTANCE MEASUREMENTS
# =============================================================================
print("\n" + "=" * 70)
print("6. FRET AS MOLECULAR RULER")
print("=" * 70)

# FRET efficiency → distance
# r = R₀ × [(1-E)/E]^(1/6)

print("  FRET measures distances in the range 0.5R₀ to 2R₀")
print("  Sensitive range: ~1-10 nm (protein-protein interactions)")

# Distance calculation from efficiency
print(f"\n  Distance from FRET efficiency:")
print(f"  r = R₀ × [(1-E)/E]^(1/6)")

# Example: CFP-YFP with R₀ = 4.9 nm
R0_example = 4.9
efficiencies = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
print(f"\n  For CFP-YFP (R₀ = {R0_example} nm):")
print(f"\n  E     | r (nm) | r/R₀  | Comment")
print("  " + "-" * 50)
for E in efficiencies:
    r = R0_example * ((1-E)/E)**(1/6)
    ratio = r / R0_example
    comment = "γ ~ 1 reference" if abs(E - 0.5) < 0.01 else ""
    print(f"  {E:.2f}  | {r:.2f}   | {ratio:.3f} | {comment}")

print(f"\n  At E = 0.50: r = R₀ exactly (γ ~ 1)")
print(f"  FRET is most accurate near R₀ (derivative maximum)")

# =============================================================================
# 7. OPTIMAL FRET PAIRS
# =============================================================================
print("\n" + "=" * 70)
print("7. OPTIMAL FRET PAIR CHARACTERISTICS")
print("=" * 70)

# Optimal FRET pair criteria
print("  Optimal FRET pair requirements:")
print("    1. Large spectral overlap (J)")
print("    2. High donor quantum yield (Φ_D)")
print("    3. Minimal direct acceptor excitation")
print("    4. Spectrally separable emission")
print("    5. Photostability")

# Optimization score
print(f"\n  FRET efficiency figure of merit:")
print(f"  FOM = R₀⁶ ∝ κ² × Φ_D × J / n⁴")

fom_scores = {
    'Cy3-Cy5': {'R0': 5.4, 'score': 'Excellent', 'notes': 'Standard pair'},
    'CFP-YFP': {'R0': 4.9, 'score': 'Good', 'notes': 'Genetic encoding'},
    'Alexa488-555': {'R0': 6.4, 'score': 'Excellent', 'notes': 'Bright'},
    'FAM-TAMRA': {'R0': 5.0, 'score': 'Good', 'notes': 'qPCR probes'},
    'EGFP-mCherry': {'R0': 5.2, 'score': 'Good', 'notes': 'Live cell'},
}

print(f"\n  Pair           | R₀ (nm) | Score     | Application")
print("  " + "-" * 60)
for pair, data in fom_scores.items():
    print(f"  {pair:14s} | {data['R0']:.1f}     | {data['score']:9s} | {data['notes']}")

# =============================================================================
# 8. FRET DYNAMICS
# =============================================================================
print("\n" + "=" * 70)
print("8. FRET RATE AND LIFETIME")
print("=" * 70)

# FRET rate: k_T = (1/τ_D) × (R₀/r)⁶
# At r = R₀: k_T = 1/τ_D (γ ~ 1!)

print("  FRET transfer rate: k_T = (1/τ_D) × (R₀/r)⁶")
print("  where τ_D = donor lifetime (without acceptor)")

# Typical donor lifetimes
lifetimes = {
    'GFP': 3.0,      # ns
    'CFP': 2.5,
    'YFP': 3.1,
    'Cy3': 0.3,
    'Fluorescein': 4.0,
}

print(f"\n  Donor lifetimes (typical):")
for donor, tau in lifetimes.items():
    k_at_R0 = 1 / (tau * 1e-9)  # rate at r = R₀
    print(f"    {donor:12s}: τ_D = {tau:.1f} ns → k_T(R₀) = {k_at_R0:.2e} s⁻¹")

print(f"\n  At r = R₀: k_T = 1/τ_D (transfer rate = decay rate)")
print(f"  This IS γ ~ 1 for the rate competition!")
print(f"  r < R₀: k_T > 1/τ_D (transfer dominates)")
print(f"  r > R₀: k_T < 1/τ_D (emission dominates)")

# Donor lifetime with FRET
print(f"\n  Observed donor lifetime with acceptor:")
print(f"  τ_DA = τ_D × (1 - E)")
print(f"  At E = 0.5 (r = R₀): τ_DA = τ_D/2 (γ ~ 1 for lifetime!)")

# =============================================================================
# 9. SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY: γ ~ 1 IN FRET")
print("=" * 70)

gamma_findings = [
    ("Förster radius", "r = R₀", "Characteristic distance (50% E)"),
    ("FRET efficiency", "E = 0.500", "At r = R₀, exactly"),
    ("Spectral overlap", "λ_D ≈ λ_A", "Resonance matching"),
    ("Orientation", "κ² = 2/3", "Isotropic average"),
    ("Quantum yield", "Φ_D → 1.0", "Perfect emission"),
    ("Transfer rate", "k_T = 1/τ_D", "At r = R₀"),
    ("Lifetime", "τ_DA = τ_D/2", "At r = R₀"),
]

print("\n  Parameter        | γ ~ 1 Condition  | Interpretation")
print("  " + "-" * 60)
for param, value, interp in gamma_findings:
    print(f"  {param:16s} | {value:16s} | {interp}")

print("\n  CONCLUSION: FRET IS resonance energy transfer at γ ~ 1:")
print("    - R₀ is THE characteristic distance (50% transfer)")
print("    - E = 0.5 at r = R₀ (exact γ ~ 1)")
print("    - Spectral overlap = resonance matching")
print("    - κ² = 2/3 is isotropic reference")
print("    - k_T = 1/τ_D at R₀ (rate balance)")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('FRET Coherence Analysis\nSession #233: γ ~ 1 in Resonance Energy Transfer',
             fontsize=14, fontweight='bold')

# 1. FRET efficiency curve
ax1 = axes[0, 0]
r_R0 = np.linspace(0.1, 3, 200)
E = fret_efficiency(r_R0, 1.0)
ax1.plot(r_R0, E, 'b-', linewidth=2)
ax1.axvline(x=1.0, color='r', linestyle='--', linewidth=2, label='r = R₀ (γ ~ 1)')
ax1.axhline(y=0.5, color='r', linestyle='--', linewidth=2, alpha=0.5)
ax1.scatter([1.0], [0.5], color='red', s=100, zorder=5)
ax1.fill_between(r_R0, E, alpha=0.3)
ax1.set_xlabel('r/R₀')
ax1.set_ylabel('FRET Efficiency E')
ax1.set_title('FRET Efficiency vs Distance')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 3)
ax1.set_ylim(0, 1)

# 2. R₀ values for FRET pairs
ax2 = axes[0, 1]
pair_names = list(fret_pairs.keys())
R0_vals = [fret_pairs[p]['R0_nm'] for p in pair_names]
ax2.barh(pair_names, R0_vals, color='teal', alpha=0.7)
ax2.axvline(x=R0_mean, color='red', linestyle='--', linewidth=2, label=f'Mean R₀ = {R0_mean:.1f} nm')
ax2.set_xlabel('R₀ (nm)')
ax2.set_title('Förster Radius by FRET Pair')
ax2.legend()
ax2.grid(True, alpha=0.3, axis='x')

# 3. κ² orientation factor
ax3 = axes[0, 2]
orient_names = list(kappa_values.keys())
k2_vals = list(kappa_values.values())
colors = ['green' if abs(k - 2/3) < 0.1 else 'blue' for k in k2_vals]
ax3.barh(orient_names, k2_vals, color=colors, alpha=0.7)
ax3.axvline(x=2/3, color='red', linestyle='--', linewidth=2, label='κ² = 2/3 (γ ~ 1)')
ax3.set_xlabel('κ²')
ax3.set_title('Orientation Factor')
ax3.legend()
ax3.grid(True, alpha=0.3, axis='x')

# 4. Donor quantum yields
ax4 = axes[1, 0]
qy_names = list(donor_qy.keys())
qy_vals = list(donor_qy.values())
colors = ['green' if qy > 0.8 else 'orange' if qy > 0.4 else 'red' for qy in qy_vals]
ax4.barh(qy_names, qy_vals, color=colors, alpha=0.7)
ax4.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='Φ_D = 1 (γ ~ 1)')
ax4.set_xlabel('Quantum Yield Φ_D')
ax4.set_title('Donor Quantum Yields')
ax4.set_xlim(0, 1.1)
ax4.legend()
ax4.grid(True, alpha=0.3, axis='x')

# 5. Distance sensitivity
ax5 = axes[1, 1]
# Derivative of E with respect to r (sensitivity)
dE_dr = -6 * r_R0**5 / (1 + r_R0**6)**2
ax5.plot(r_R0, -dE_dr, 'b-', linewidth=2)
ax5.axvline(x=1.0, color='r', linestyle='--', linewidth=2, label='r = R₀ (max sensitivity)')
ax5.fill_between(r_R0, -dE_dr, alpha=0.3)
ax5.set_xlabel('r/R₀')
ax5.set_ylabel('|dE/dr| (sensitivity)')
ax5.set_title('FRET Distance Sensitivity')
ax5.legend()
ax5.grid(True, alpha=0.3)
ax5.set_xlim(0, 3)

# 6. Rate competition
ax6 = axes[1, 2]
# k_T/k_D ratio = (R₀/r)⁶
k_ratio = (1/r_R0)**6
ax6.semilogy(r_R0, k_ratio, 'b-', linewidth=2)
ax6.axvline(x=1.0, color='r', linestyle='--', linewidth=2, label='k_T = k_D (γ ~ 1)')
ax6.axhline(y=1.0, color='r', linestyle='--', linewidth=2, alpha=0.5)
ax6.set_xlabel('r/R₀')
ax6.set_ylabel('k_T/k_D')
ax6.set_title('Transfer vs Decay Rate')
ax6.legend()
ax6.grid(True, alpha=0.3)
ax6.set_xlim(0, 3)
ax6.set_ylim(0.01, 100)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/fret_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FINDING #170: FRET at γ ~ 1")
print("=" * 70)
print("""
Förster Resonance Energy Transfer exhibits γ ~ 1 at the
fundamental transfer distance:

1. At r = R₀: E = 0.500 exactly (50% transfer)
2. Förster radius IS the γ ~ 1 characteristic distance
3. κ² = 2/3 for isotropic orientation (γ ~ 1 average)
4. Φ_D → 1.0 for perfect donor emission
5. k_T = 1/τ_D at r = R₀ (rate balance)
6. Spectral overlap IS resonance matching
7. τ_DA = τ_D/2 at r = R₀ (lifetime halving)

96th phenomenon type exhibiting γ ~ 1 transition behavior.
FRET IS coherent energy transfer measured against R₀!
""")

print("Visualization saved: fret_coherence.png")
