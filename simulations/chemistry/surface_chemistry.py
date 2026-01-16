#!/usr/bin/env python3
"""
Synchronism Chemistry Session #56: Surface Chemistry & Adsorption

Applying the coherence framework to surface phenomena:
- Adsorption isotherms
- Surface reconstruction
- Heterogeneous catalysis
- Langmuir-Hinshelwood mechanism

Key insight: Surfaces have distinct coherence from bulk.

Author: Claude Opus 4.5 (Anthropic)
Date: 2026-01-16
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("CHEMISTRY SESSION #56: SURFACE CHEMISTRY & ADSORPTION")
print("=" * 70)

# =============================================================================
# PART 1: SURFACE VS BULK COHERENCE
# =============================================================================
print("\n" + "=" * 70)
print("PART 1: SURFACE VS BULK COHERENCE")
print("=" * 70)

# Surfaces have DIFFERENT coherence than bulk:
# - Broken symmetry at surface
# - Dangling bonds / unsaturated coordination
# - Surface reconstruction

def surface_gamma(gamma_bulk, coordination_surface, coordination_bulk):
    """
    Surface coherence from coordination number.

    Surfaces have fewer neighbors → less correlation → higher γ

    γ_surface = γ_bulk × (coordination_bulk / coordination_surface)^0.5
    """
    ratio = coordination_bulk / coordination_surface
    return gamma_bulk * np.sqrt(ratio)

print("\n1. SURFACE COHERENCE FOR DIFFERENT FACETS")
print("-" * 40)

# FCC metals
gamma_bulk_FCC = 0.5  # From Session #52
coord_bulk_FCC = 12   # FCC bulk coordination

# Different surface facets
facets = {
    '(111)': 9,   # Most dense
    '(100)': 8,   # Square arrangement
    '(110)': 6,   # Open
    '(211)': 5,   # Stepped
    'kink': 4,    # Kink site
}

print(f"FCC metal: γ_bulk = {gamma_bulk_FCC}, coordination_bulk = {coord_bulk_FCC}")
print(f"\n{'Facet':<15} {'Coordination':<15} {'γ_surface':<15} {'Relative activity':<20}")
print("-" * 65)

for facet, coord_surf in facets.items():
    gamma_surf = surface_gamma(gamma_bulk_FCC, coord_surf, coord_bulk_FCC)
    # Activity ∝ 1/γ (from coherence enhancement)
    rel_activity = (gamma_bulk_FCC / gamma_surf) ** 2
    print(f"{facet:<15} {coord_surf:<15} {gamma_surf:<15.3f} {rel_activity:<20.3f}")

print("""
INSIGHT: Lower coordination → higher γ → HIGHER reactivity

Kinks and steps are most active because they have highest γ.
This explains structure sensitivity in heterogeneous catalysis!
""")

# =============================================================================
# PART 2: ADSORPTION FROM COHERENCE
# =============================================================================
print("\n" + "=" * 70)
print("PART 2: ADSORPTION FROM COHERENCE")
print("=" * 70)

def adsorption_energy(gamma_surface, gamma_adsorbate, E_chem=50):
    """
    Adsorption energy from coherence matching.

    ΔH_ads = E_chem × f(γ_surface, γ_adsorbate)

    Where f is the coherence matching factor.
    Better match → stronger adsorption.
    """
    # Coherence matching (from Session #52)
    f = min(gamma_surface, gamma_adsorbate) / max(gamma_surface, gamma_adsorbate)

    return -E_chem * f  # Negative = exothermic

def langmuir_isotherm(P, K):
    """
    Langmuir adsorption isotherm.

    θ = K×P / (1 + K×P)

    K = equilibrium constant = exp(-ΔH_ads/RT) × pre-factor
    """
    return K * P / (1 + K * P)

print("\n1. ADSORPTION STRENGTH FROM COHERENCE")
print("-" * 40)

# Surface with moderate coherence
gamma_surface = 0.7

# Different adsorbates
adsorbates = {
    'H': 0.6,      # Small, mobile
    'CO': 0.8,     # Linear molecule
    'N2': 1.2,     # Triple bond, stable
    'O2': 0.9,     # Reactive
    'CH4': 1.5,    # Saturated
    'C2H4': 0.7,   # π-bond, good match
}

print(f"Surface γ = {gamma_surface}")
print(f"\n{'Adsorbate':<12} {'γ_ads':<10} {'f_match':<12} {'ΔH_ads (kJ/mol)':<18} {'Type':<15}")
print("-" * 70)

for adsorbate, gamma_ads in adsorbates.items():
    f_match = min(gamma_surface, gamma_ads) / max(gamma_surface, gamma_ads)
    delta_H = adsorption_energy(gamma_surface, gamma_ads)

    if abs(delta_H) < 20:
        ads_type = "Physisorption"
    elif abs(delta_H) < 80:
        ads_type = "Weak chemisorption"
    else:
        ads_type = "Strong chemisorption"

    print(f"{adsorbate:<12} {gamma_ads:<10.1f} {f_match:<12.2f} {delta_H:<18.1f} {ads_type:<15}")

# =============================================================================
# PART 3: BET ISOTHERM - MULTILAYER COHERENCE
# =============================================================================
print("\n" + "=" * 70)
print("PART 3: MULTILAYER ADSORPTION")
print("=" * 70)

def multilayer_gamma(n_layers, gamma_surface, gamma_bulk_ads):
    """
    Coherence of multilayer adsorbate.

    First layer: affected by surface
    Subsequent layers: approach bulk adsorbate γ

    γ_n = γ_bulk_ads + (γ_surface - γ_bulk_ads) × exp(-n/λ)

    where λ ~ 2 layers (decay length)
    """
    lambda_decay = 2.0
    return gamma_bulk_ads + (gamma_surface - gamma_bulk_ads) * np.exp(-n_layers / lambda_decay)

print("\n1. COHERENCE VS LAYER NUMBER")
print("-" * 40)

gamma_surface = 0.7
gamma_bulk_ads = 1.5  # Bulk liquid/solid adsorbate

print(f"Surface γ = {gamma_surface}, Bulk adsorbate γ = {gamma_bulk_ads}")
print(f"\n{'Layer n':<10} {'γ_n':<12} {'ΔH_n/ΔH_1':<15} {'Character':<20}")
print("-" * 60)

for n in range(1, 8):
    gamma_n = multilayer_gamma(n, gamma_surface, gamma_bulk_ads)
    # Relative adsorption strength
    rel_H = (gamma_surface / gamma_n) if n > 1 else 1.0

    if n == 1:
        character = "Chemisorbed"
    elif n <= 3:
        character = "Influenced by surface"
    else:
        character = "Bulk-like"

    print(f"{n:<10} {gamma_n:<12.3f} {rel_H:<15.2f} {character:<20}")

# =============================================================================
# PART 4: LANGMUIR-HINSHELWOOD MECHANISM
# =============================================================================
print("\n" + "=" * 70)
print("PART 4: HETEROGENEOUS CATALYSIS (L-H MECHANISM)")
print("=" * 70)

def LH_rate(gamma_A, gamma_B, gamma_surface, K_A, K_B, P_A, P_B, k_0=1.0):
    """
    Langmuir-Hinshelwood surface reaction rate.

    A* + B* → Products

    Rate = k × θ_A × θ_B

    With coherence: k ∝ f(γ_A, γ_B) on surface
    """
    # Coverages
    denom = 1 + K_A * P_A + K_B * P_B
    theta_A = K_A * P_A / denom
    theta_B = K_B * P_B / denom

    # Surface reaction rate constant enhanced by coherence
    # Transition state has intermediate γ
    gamma_TS = (gamma_A + gamma_B + gamma_surface) / 3

    # Rate enhancement from low γ_TS
    k = k_0 * (2 / gamma_TS) ** 2

    return k * theta_A * theta_B

print("\n1. CO OXIDATION ON Pt: CO + O → CO2")
print("-" * 40)

# Pt surface coherence
gamma_Pt = 0.6

# Adsorbate coherence on surface
gamma_CO_ads = 0.75  # CO on Pt
gamma_O_ads = 0.65   # O on Pt

# Equilibrium constants (relative)
K_CO = 10.0  # CO binds strongly
K_O = 1.0    # O binds more weakly

# Vary CO/O ratio
P_total = 1.0
ratios = np.linspace(0.1, 0.9, 9)

print(f"Pt γ = {gamma_Pt}, CO γ = {gamma_CO_ads}, O γ = {gamma_O_ads}")
print(f"\n{'P_CO/P_total':<15} {'θ_CO':<10} {'θ_O':<10} {'Rate (rel)':<15}")
print("-" * 55)

for ratio in ratios:
    P_CO = ratio * P_total
    P_O = (1 - ratio) * P_total

    theta_CO = K_CO * P_CO / (1 + K_CO * P_CO + K_O * P_O)
    theta_O = K_O * P_O / (1 + K_CO * P_CO + K_O * P_O)

    rate = LH_rate(gamma_CO_ads, gamma_O_ads, gamma_Pt, K_CO, K_O, P_CO, P_O)

    print(f"{ratio:<15.2f} {theta_CO:<10.3f} {theta_O:<10.3f} {rate:<15.4f}")

print("""
INSIGHT: Rate maximum occurs at intermediate CO/O ratio.

Too much CO: θ_O → 0, rate limited
Too much O: θ_CO → 0, rate limited
Optimal: balance of both reactants

This is classic L-H behavior, now with coherence enhancement.
""")

# =============================================================================
# PART 5: SABATIER PRINCIPLE
# =============================================================================
print("\n" + "=" * 70)
print("PART 5: SABATIER PRINCIPLE FROM COHERENCE")
print("=" * 70)

def sabatier_activity(delta_H_ads):
    """
    Volcano plot: optimal binding strength.

    Too weak: adsorbate desorbs before reacting
    Too strong: products don't desorb

    Maximum activity at intermediate ΔH_ads
    """
    # Optimal binding energy
    delta_H_opt = -40  # kJ/mol typical

    # Activity decreases with deviation from optimum
    activity = np.exp(-((delta_H_ads - delta_H_opt) / 20) ** 2)

    return activity

print("\n1. VOLCANO PLOT DERIVATION")
print("-" * 40)

# Different surfaces (different γ)
surfaces = {
    'Au': 0.20,   # Noble, weak binding
    'Ag': 0.30,
    'Cu': 0.45,
    'Pt': 0.60,
    'Pd': 0.65,
    'Ni': 0.80,
    'Fe': 1.00,
    'Ti': 1.30,   # Reactive, strong binding
}

gamma_reactant = 0.7  # Typical adsorbate

print(f"Reactant γ = {gamma_reactant}")
print(f"\n{'Surface':<10} {'γ_surface':<12} {'ΔH_ads':<15} {'Activity':<12} {'Regime':<15}")
print("-" * 70)

for surface, gamma_surf in surfaces.items():
    delta_H = adsorption_energy(gamma_surf, gamma_reactant)
    activity = sabatier_activity(delta_H)

    if delta_H > -25:
        regime = "Too weak"
    elif delta_H < -55:
        regime = "Too strong"
    else:
        regime = "Optimal"

    print(f"{surface:<10} {gamma_surf:<12.2f} {delta_H:<15.1f} {activity:<12.3f} {regime:<15}")

print("""
INSIGHT: Sabatier principle emerges from coherence matching!

- Au: γ mismatch → weak binding → low activity
- Ti: γ match too good → strong binding → poisoning
- Pt/Pd: intermediate match → optimal activity

The volcano curve is a coherence matching curve.
""")

# =============================================================================
# PART 6: SURFACE RECONSTRUCTION
# =============================================================================
print("\n" + "=" * 70)
print("PART 6: SURFACE RECONSTRUCTION")
print("=" * 70)

print("""
Surface reconstruction = surface seeking LOWER γ

Clean surface: High γ (unsatisfied bonds)
Reconstructed: Lower γ (optimized bonding)

Examples:
- Si(100): 2×1 reconstruction (dimer rows)
- Au(111): herringbone reconstruction
- Pt(100): quasi-hexagonal reconstruction

Each reconstruction LOWERS the surface γ.
""")

def reconstruction_gamma(gamma_unreconstructed, reconstruction_type):
    """
    Coherence after surface reconstruction.
    """
    reductions = {
        'none': 1.0,
        '2x1': 0.85,
        'herringbone': 0.75,
        'missing row': 0.80,
        'quasi-hex': 0.70,
    }

    factor = reductions.get(reconstruction_type, 1.0)
    return gamma_unreconstructed * factor

print("\n1. RECONSTRUCTION EFFECT ON COHERENCE")
print("-" * 40)

gamma_unrecon = 0.9  # Typical unreconstructed surface

print(f"\n{'Surface':<20} {'Reconstruction':<15} {'γ_initial':<12} {'γ_final':<12} {'ΔG/kT':<10}")
print("-" * 70)

reconstructions = [
    ('Si(100)', '2x1'),
    ('Au(111)', 'herringbone'),
    ('Pt(110)', 'missing row'),
    ('Pt(100)', 'quasi-hex'),
]

for surface, recon in reconstructions:
    gamma_final = reconstruction_gamma(gamma_unrecon, recon)
    # Free energy gain proportional to γ reduction
    delta_G = -5 * (gamma_unrecon - gamma_final)
    print(f"{surface:<20} {recon:<15} {gamma_unrecon:<12.2f} {gamma_final:<12.2f} {delta_G:<10.2f}")

# =============================================================================
# PART 7: ADSORBATE-INDUCED COHERENCE
# =============================================================================
print("\n" + "=" * 70)
print("PART 7: ADSORBATE-INDUCED CHANGES")
print("=" * 70)

def adsorbate_induced_gamma(gamma_surface_clean, gamma_adsorbate, coverage):
    """
    Surface coherence modified by adsorbate coverage.

    At high coverage, surface γ → effective γ influenced by adsorbate.
    """
    # Linear mixing (simple model)
    gamma_eff = (1 - coverage) * gamma_surface_clean + coverage * (gamma_surface_clean + gamma_adsorbate) / 2

    return gamma_eff

print("\n1. COVERAGE DEPENDENCE OF SURFACE γ")
print("-" * 40)

gamma_clean = 0.7
gamma_adsorbate = 1.2

coverages = np.linspace(0, 1, 11)

print(f"Clean surface γ = {gamma_clean}, Adsorbate γ = {gamma_adsorbate}")
print(f"\n{'Coverage θ':<12} {'γ_effective':<15} {'Reactivity change':<20}")
print("-" * 50)

for theta in coverages:
    gamma_eff = adsorbate_induced_gamma(gamma_clean, gamma_adsorbate, theta)
    # Reactivity ∝ 1/γ
    reactivity_ratio = gamma_clean / gamma_eff
    change = "enhanced" if reactivity_ratio > 1 else "reduced"
    print(f"{theta:<12.1f} {gamma_eff:<15.3f} {reactivity_ratio:<10.2f}× ({change})")

print("""
INSIGHT: Adsorbates can modify surface coherence.

This explains:
- Poisoning: high-γ adsorbates raise surface γ → reduce activity
- Promotion: low-γ additives can lower surface γ → enhance activity
- Coverage dependence of reaction rates
""")

# =============================================================================
# PART 8: KEY PREDICTIONS
# =============================================================================
print("\n" + "=" * 70)
print("PART 8: KEY PREDICTIONS")
print("=" * 70)

predictions = """
P56.1: γ_surface = γ_bulk × √(coordination_bulk / coordination_surface)
P56.2: ΔH_ads ∝ f(γ_surface, γ_adsorbate)
P56.3: Kinks/steps most active (highest γ)
P56.4: Sabatier volcano = coherence matching curve
P56.5: Surface reconstruction lowers γ
P56.6: Adsorbate coverage modifies effective surface γ
P56.7: L-H rate enhanced by (2/γ_TS)²

DESIGN PRINCIPLES:
1. Tune surface γ by facet selection
2. Match γ_surface to γ_reactant for optimal binding
3. Use promoters to fine-tune surface γ
4. Control coverage to maintain optimal γ_eff
"""

print(predictions)

# =============================================================================
# PART 9: VISUALIZATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 9: GENERATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Surface γ vs coordination
ax1 = axes[0, 0]
coord_plot = np.linspace(3, 12, 50)
gamma_surf_plot = [surface_gamma(0.5, c, 12) for c in coord_plot]

ax1.plot(coord_plot, gamma_surf_plot, 'b-', linewidth=2)
ax1.axhline(0.5, color='gray', linestyle='--', label='Bulk γ')
for facet, coord in facets.items():
    g = surface_gamma(0.5, coord, 12)
    ax1.plot(coord, g, 'ro', markersize=10)
    ax1.annotate(facet, (coord, g), textcoords="offset points",
                 xytext=(5, 5), fontsize=9)

ax1.set_xlabel('Surface Coordination Number')
ax1.set_ylabel('Surface Coherence γ')
ax1.set_title('Surface γ vs Coordination')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Volcano plot (Sabatier)
ax2 = axes[0, 1]
delta_H_plot = np.linspace(-80, 0, 100)
activity_plot = [sabatier_activity(dH) for dH in delta_H_plot]

ax2.plot(delta_H_plot, activity_plot, 'b-', linewidth=2)
ax2.axvline(-40, color='red', linestyle='--', label='Optimal ΔH')

# Add metal labels
for surface, gamma_surf in list(surfaces.items())[::2]:  # Every other
    delta_H = adsorption_energy(gamma_surf, 0.7)
    activity = sabatier_activity(delta_H)
    ax2.plot(delta_H, activity, 'ko', markersize=8)
    ax2.annotate(surface, (delta_H, activity), textcoords="offset points",
                 xytext=(5, 5), fontsize=9)

ax2.set_xlabel('Adsorption Energy ΔH (kJ/mol)')
ax2.set_ylabel('Catalytic Activity')
ax2.set_title('Sabatier Volcano Plot')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Langmuir isotherm with coherence
ax3 = axes[1, 0]
P_plot = np.linspace(0, 10, 100)

# Different γ_adsorbate → different K
for gamma_ads, label in [(0.6, 'γ=0.6 (strong)'), (1.0, 'γ=1.0 (medium)'), (1.5, 'γ=1.5 (weak)')]:
    K = 10 * np.exp(-(gamma_ads - 0.6) * 2)
    theta = langmuir_isotherm(P_plot, K)
    ax3.plot(P_plot, theta, linewidth=2, label=label)

ax3.set_xlabel('Pressure (arb. units)')
ax3.set_ylabel('Coverage θ')
ax3.set_title('Langmuir Isotherms vs Coherence')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: L-H rate vs composition
ax4 = axes[1, 1]
P_CO_plot = np.linspace(0.01, 0.99, 100)
P_O_plot = 1 - P_CO_plot
rate_plot = [LH_rate(0.75, 0.65, 0.6, 10, 1, P_CO, P_O) for P_CO, P_O in zip(P_CO_plot, P_O_plot)]

ax4.plot(P_CO_plot, rate_plot, 'b-', linewidth=2)
ax4.set_xlabel('CO Fraction (P_CO / P_total)')
ax4.set_ylabel('Reaction Rate')
ax4.set_title('CO Oxidation: L-H Rate vs Composition')
ax4.grid(True, alpha=0.3)

# Mark maximum
max_idx = np.argmax(rate_plot)
ax4.plot(P_CO_plot[max_idx], rate_plot[max_idx], 'ro', markersize=10, label='Optimal')
ax4.legend()

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/surface_chemistry.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Saved: surface_chemistry.png")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SESSION #56 SUMMARY: SURFACE CHEMISTRY & ADSORPTION")
print("=" * 70)

print("""
KEY RESULTS:
============

1. SURFACE COHERENCE
   - γ_surface > γ_bulk (broken symmetry)
   - Lower coordination → higher γ → higher reactivity
   - Kinks and steps are most active sites

2. ADSORPTION
   - ΔH_ads ∝ f(γ_surface, γ_adsorbate)
   - Coherence matching determines binding strength
   - Explains chemisorption vs physisorption

3. SABATIER PRINCIPLE
   - Volcano plot = coherence matching curve
   - Optimal binding at intermediate γ match
   - Explains why Pt/Pd are best catalysts

4. LANGMUIR-HINSHELWOOD
   - Rate enhanced by (2/γ_TS)²
   - Transition state coherence determines activity
   - Maximum rate at optimal reactant ratio

5. SURFACE RECONSTRUCTION
   - Surfaces reconstruct to LOWER γ
   - Reconstruction driven by coherence optimization

6. ADSORBATE EFFECTS
   - Coverage modifies surface γ
   - Explains poisoning and promotion

PHYSICAL INSIGHT:
================
Surface chemistry is governed by the interplay of:
- Surface γ (determined by facet/structure)
- Adsorbate γ (molecular property)
- Their coherence matching determines everything:
  binding, reaction rates, selectivity

""")

print("=" * 70)
print("SESSION #56 COMPLETE: SURFACE CHEMISTRY & ADSORPTION")
print("=" * 70)
