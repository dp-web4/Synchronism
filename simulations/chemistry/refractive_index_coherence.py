#!/usr/bin/env python3
"""
Chemistry Session #76: Refractive Index & Coherence
Test whether coherence framework predicts optical properties.

Refractive index (n) measures how light slows in a medium:
- n = c/v (ratio of vacuum to medium light speed)
- Related to dielectric constant: n² ≈ ε_r (non-magnetic materials)

The Lorentz-Lorenz equation:
(n² - 1)/(n² + 2) = (4π/3) × N × α

Where α is molecular polarizability.

Coherence interpretation:
- Polarizability α measures how easily electrons are displaced
- More coherent electrons (lower γ) → harder to displace → lower α
- But also: higher coherence → more collective response → could increase n

This is nuanced! Let's explore the relationship.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

print("=" * 70)
print("CHEMISTRY SESSION #76: REFRACTIVE INDEX & COHERENCE")
print("=" * 70)

# ==============================================================================
# DATASET: REFRACTIVE INDICES
# ==============================================================================

print("\n" + "=" * 70)
print("DATASET: REFRACTIVE INDICES AT 589nm (sodium D line)")
print("=" * 70)

# Refractive index data at 25°C
# Format: (n, Hildebrand δ, band gap eV for semiconductors)
# n measured at sodium D line (589nm)

liquids_n = {
    # Low n (nonpolar)
    'hexane': (1.375, 14.9),
    'pentane': (1.358, 14.4),
    'cyclohexane': (1.426, 16.8),
    'diethyl_ether': (1.353, 15.1),

    # Medium n
    'water': (1.333, 47.8),
    'methanol': (1.329, 29.6),
    'ethanol': (1.361, 26.5),
    'acetone': (1.359, 20.0),

    # Higher n
    'benzene': (1.501, 18.8),
    'toluene': (1.497, 18.2),
    'chloroform': (1.446, 19.0),
    'carbon_tetrachloride': (1.460, 17.8),
    'carbon_disulfide': (1.627, 20.5),
    'nitrobenzene': (1.553, 22.1),

    # High n
    'bromoform': (1.595, 21.4),
    'methylene_iodide': (1.742, 24.0),
}

# Transparent solids (glasses, crystals)
solids_n = {
    # Glasses
    'fused_silica': (1.458, 9.0),    # SiO2, E_g ~ 9 eV
    'crown_glass': (1.520, 7.0),      # approximate
    'flint_glass': (1.620, 5.0),      # heavy metals increase n

    # Crystals
    'NaCl': (1.544, 8.6),            # E_g ~ 8.6 eV
    'CaF2': (1.434, 12.0),           # E_g ~ 12 eV
    'MgF2': (1.378, 13.0),           # Very high gap
    'Al2O3': (1.768, 8.8),           # Sapphire

    # Semiconductors (below band gap)
    'Si': (3.42, 1.12),              # E_g = 1.12 eV
    'Ge': (4.00, 0.67),              # E_g = 0.67 eV
    'GaAs': (3.30, 1.42),            # E_g = 1.42 eV
    'ZnS': (2.36, 3.6),              # E_g ~ 3.6 eV
    'ZnSe': (2.67, 2.7),             # E_g ~ 2.7 eV
    'diamond': (2.417, 5.5),         # E_g = 5.5 eV
    'TiO2': (2.613, 3.2),            # Rutile
    'ZrO2': (2.16, 5.0),             # Zirconia
}

print(f"Liquids: {len(liquids_n)}")
print(f"Solids: {len(solids_n)}")

# Print sorted by n
print("\nLiquids sorted by refractive index:")
print("-" * 50)
for name, (n, delta) in sorted(liquids_n.items(), key=lambda x: x[1][0]):
    print(f"{name:<25}: n = {n:.3f}, δ = {delta:.1f}")

print("\nSolids sorted by refractive index:")
print("-" * 50)
for name, (n, E_g) in sorted(solids_n.items(), key=lambda x: x[1][0]):
    print(f"{name:<25}: n = {n:.3f}, E_g = {E_g:.2f} eV")

# ==============================================================================
# SEMICONDUCTOR BAND GAP RELATIONSHIP
# ==============================================================================

print("\n" + "=" * 70)
print("BAND GAP - REFRACTIVE INDEX RELATIONSHIP")
print("=" * 70)

print("""
Moss's rule (empirical):
E_g × n^4 ≈ constant (≈ 95 eV for many semiconductors)

Or: n ∝ E_g^(-1/4)

Coherence interpretation:
- Smaller band gap → electrons more easily excited
- Lower E_g → higher polarizability → higher n
- E_g ∝ 2/γ (from Session #60)
- Therefore: n ∝ γ^(1/4)
""")

# Extract semiconductor data
sc_names = []
sc_n = []
sc_Eg = []

for name, (n, E_g) in solids_n.items():
    if name in ['Si', 'Ge', 'GaAs', 'ZnS', 'ZnSe', 'diamond', 'TiO2', 'ZrO2']:
        sc_names.append(name)
        sc_n.append(n)
        sc_Eg.append(E_g)

sc_n_arr = np.array(sc_n)
sc_Eg_arr = np.array(sc_Eg)

# Test Moss's rule
moss_product = sc_Eg_arr * sc_n_arr**4
print(f"\nMoss's rule test (E_g × n^4):")
print("-" * 50)
for i, name in enumerate(sc_names):
    print(f"{name:<12}: n = {sc_n_arr[i]:.2f}, E_g = {sc_Eg_arr[i]:.2f}, E_g×n^4 = {moss_product[i]:.1f}")

print(f"\nMean E_g×n^4 = {moss_product.mean():.1f} ± {moss_product.std():.1f}")

# Correlation: n vs E_g
r_n_Eg, p_n_Eg = stats.pearsonr(sc_n_arr, sc_Eg_arr)
print(f"\nn vs E_g: r = {r_n_Eg:.3f}")

# n vs E_g^(-1/4)
Eg_inv_fourth = sc_Eg_arr**(-0.25)
r_n_Eg_inv, _ = stats.pearsonr(sc_n_arr, Eg_inv_fourth)
print(f"n vs E_g^(-1/4): r = {r_n_Eg_inv:.3f}")

# ==============================================================================
# COHERENCE INTERPRETATION
# ==============================================================================

print("\n" + "=" * 70)
print("COHERENCE INTERPRETATION")
print("=" * 70)

def gamma_from_Eg(E_g, E_ref=3.0):
    """
    Estimate γ from band gap.
    Larger E_g → more tightly bound electrons → lower γ (more coherent).
    """
    # Using relationship from Session #60: E_g ∝ 2/γ
    # So γ ∝ 2/E_g
    gamma = 2.0 * E_ref / E_g if E_g > 0 else 2.0
    return np.clip(gamma, 0.5, 2.0)

def gamma_from_delta(delta, delta_ref=25.0):
    """
    Estimate γ from Hildebrand parameter for liquids.
    """
    gamma = 2.0 - 1.5 * (delta / 50.0)
    return np.clip(gamma, 0.5, 2.0)

# Calculate γ for semiconductors
print("\nγ from band gap (semiconductors):")
print("-" * 50)
gamma_sc_list = []
for i, name in enumerate(sc_names):
    gamma = gamma_from_Eg(sc_Eg_arr[i])
    gamma_sc_list.append(gamma)
    print(f"{name:<12}: E_g = {sc_Eg_arr[i]:.2f}, γ = {gamma:.2f}")

gamma_sc_arr = np.array(gamma_sc_list)

# Correlation: n vs γ
r_n_gamma, _ = stats.pearsonr(sc_n_arr, gamma_sc_arr)
print(f"\nn vs γ (semiconductors): r = {r_n_gamma:.3f}")

# n vs 2/γ
coh_factor_sc = 2.0 / gamma_sc_arr
r_n_coh, _ = stats.pearsonr(sc_n_arr, coh_factor_sc)
print(f"n vs 2/γ (semiconductors): r = {r_n_coh:.3f}")

# ==============================================================================
# LIQUID ANALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("LIQUID REFRACTIVE INDEX ANALYSIS")
print("=" * 70)

# Extract liquid data
liq_names = list(liquids_n.keys())
liq_n = np.array([liquids_n[name][0] for name in liq_names])
liq_delta = np.array([liquids_n[name][1] for name in liq_names])

# γ from Hildebrand
liq_gamma = np.array([gamma_from_delta(d) for d in liq_delta])

# Correlations
r_liq_delta, _ = stats.pearsonr(liq_n, liq_delta)
r_liq_gamma, _ = stats.pearsonr(liq_n, liq_gamma)

print(f"n vs δ (liquids): r = {r_liq_delta:.3f}")
print(f"n vs γ (liquids): r = {r_liq_gamma:.3f}")

# Actually, refractive index correlates with POLARIZABILITY
# which is related to electron mobility, not δ directly
# Let's check molecular properties

print("""
Note: For liquids, refractive index relates to:
1. Molecular polarizability (α)
2. Molecular density
3. Electron delocalization (aromatic compounds have higher n)

Aromatic compounds (benzene, toluene) have high n due to
π-electron delocalization, NOT high δ.

Carbon disulfide has high n due to large, polarizable S atoms.
Iodomethane has very high n due to large I atom.

The relationship n-δ is weak because they measure different things:
- δ = cohesive energy density (intermolecular)
- n = electronic polarizability (intramolecular)
""")

# ==============================================================================
# CLAUSIUS-MOSSOTTI / LORENTZ-LORENZ
# ==============================================================================

print("\n" + "=" * 70)
print("LORENTZ-LORENZ EQUATION + COHERENCE")
print("=" * 70)

print("""
Lorentz-Lorenz equation:
(n² - 1)/(n² + 2) = (4π/3) × N × α = R_m/V_m

Where R_m is molar refractivity.

This separates:
- Intrinsic molecular property (α, polarizability)
- Concentration effect (N, number density)

For coherence interpretation:
α ∝ (electron mobility) ∝ 1/(binding strength) ∝ γ/2

So:
- High coherence (low γ) → low α → low n
- Low coherence (high γ) → high α → high n

This is OPPOSITE to many other properties!
""")

# Calculate Lorentz-Lorenz function for liquids
def lorentz_lorenz(n):
    return (n**2 - 1) / (n**2 + 2)

liq_LL = lorentz_lorenz(liq_n)
sc_LL = lorentz_lorenz(sc_n_arr)

print("\nLorentz-Lorenz function (n² - 1)/(n² + 2):")
print("-" * 40)
print("Liquids:")
for i, name in enumerate(liq_names[:5]):
    print(f"  {name:<20}: n = {liq_n[i]:.3f}, LL = {liq_LL[i]:.3f}")

print("\nSemiconductors:")
for i, name in enumerate(sc_names):
    print(f"  {name:<12}: n = {sc_n_arr[i]:.2f}, LL = {sc_LL[i]:.3f}")

# ==============================================================================
# DISPERSION & COHERENCE
# ==============================================================================

print("\n" + "=" * 70)
print("DISPERSION (WAVELENGTH DEPENDENCE)")
print("=" * 70)

print("""
Cauchy's equation:
n(λ) = A + B/λ² + C/λ^4 + ...

Near absorption edge:
n increases strongly (anomalous dispersion)

Coherence interpretation:
- Near resonance (band edge), coherent coupling increases
- This enhances refractive index
- Dispersion reflects γ changing with photon energy

At ℏω → E_g:
- Strong resonant enhancement
- γ_optical → 0 (resonant coherence)
- n → maximum (before absorption onset)
""")

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SESSION #76 SUMMARY: REFRACTIVE INDEX & COHERENCE")
print("=" * 70)

print(f"""
Correlations Found:

SEMICONDUCTORS:
- n vs E_g: r = {r_n_Eg:.3f} (strong negative - expected)
- n vs E_g^(-1/4): r = {r_n_Eg_inv:.3f} (Moss's rule)
- n vs γ: r = {r_n_gamma:.3f}
- n vs 2/γ: r = {r_n_coh:.3f}
- Moss's rule: E_g×n^4 = {moss_product.mean():.0f} ± {moss_product.std():.0f}

LIQUIDS:
- n vs δ: r = {r_liq_delta:.3f} ({"WEAK" if abs(r_liq_delta) < 0.4 else "MODERATE"})
- n vs γ: r = {r_liq_gamma:.3f}

Key Findings:
1. For semiconductors, n correlates with γ (r = {r_n_gamma:.3f})
   - Smaller E_g → larger γ → higher n
   - This is mediated through band gap relationship

2. For liquids, n-δ correlation is weak
   - n relates to polarizability (intramolecular)
   - δ relates to cohesion (intermolecular)
   - Different physics!

3. Moss's rule (E_g×n^4 ≈ constant) approximately holds
   - Through coherence: n ∝ γ^(1/4) and E_g ∝ 2/γ
   - Product E_g × n^4 ∝ (2/γ) × γ = constant (approximately)

Physical Interpretation:
- Refractive index measures how strongly electrons respond to light
- More tightly bound electrons (low γ, high E_g) → lower polarizability → lower n
- n ∝ γ^(1/4) approximately (through Moss's rule)
""")

# ==============================================================================
# PREDICTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print("""
P76.1: n ∝ γ^(1/4) for semiconductors
From Moss's rule: E_g × n^4 ≈ constant, E_g ∝ 2/γ.

P76.2: High-n materials have small E_g (large γ)
Ge (n=4.0, E_g=0.67) vs diamond (n=2.4, E_g=5.5).

P76.3: Dispersion reflects γ(ω)
Near absorption: γ → 0 (resonance), n → maximum.

P76.4: For liquids, n correlates weakly with δ
Polarizability ≠ cohesive energy.

P76.5: Aromatic/conjugated → high n
π-electron delocalization increases polarizability.
""")

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: n vs E_g for semiconductors
ax1 = axes[0, 0]
ax1.scatter(sc_Eg_arr, sc_n_arr, s=100, alpha=0.7, c='blue')
for i, name in enumerate(sc_names):
    ax1.annotate(name, (sc_Eg_arr[i], sc_n_arr[i]), fontsize=9)
# Fit line
z = np.polyfit(sc_Eg_arr, sc_n_arr, 1)
p = np.poly1d(z)
x_line = np.linspace(sc_Eg_arr.min(), sc_Eg_arr.max(), 100)
ax1.plot(x_line, p(x_line), 'r--', label=f'r = {r_n_Eg:.3f}')
ax1.set_xlabel('Band Gap (eV)', fontsize=12)
ax1.set_ylabel('Refractive Index n', fontsize=12)
ax1.set_title('Semiconductor n vs Band Gap', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend()

# Plot 2: Moss's rule test
ax2 = axes[0, 1]
ax2.bar(sc_names, moss_product, color='green', alpha=0.7)
ax2.axhline(y=moss_product.mean(), color='red', linestyle='--', label=f'Mean = {moss_product.mean():.0f}')
ax2.set_ylabel('E_g × n^4', fontsize=12)
ax2.set_title("Moss's Rule: E_g × n^4 ≈ constant", fontsize=14)
ax2.tick_params(axis='x', rotation=45)
ax2.legend()
ax2.grid(True, alpha=0.3, axis='y')

# Plot 3: n vs γ for semiconductors
ax3 = axes[1, 0]
ax3.scatter(gamma_sc_arr, sc_n_arr, s=100, alpha=0.7, c='purple')
for i, name in enumerate(sc_names):
    ax3.annotate(name, (gamma_sc_arr[i], sc_n_arr[i]), fontsize=9)
ax3.set_xlabel('γ (from E_g)', fontsize=12)
ax3.set_ylabel('Refractive Index n', fontsize=12)
ax3.set_title(f'Semiconductor n vs γ\n(r = {r_n_gamma:.3f})', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: Liquid n vs δ
ax4 = axes[1, 1]
ax4.scatter(liq_delta, liq_n, s=80, alpha=0.7, c='orange')
# Highlight some
for i, name in enumerate(liq_names):
    if name in ['water', 'benzene', 'carbon_disulfide', 'methylene_iodide']:
        ax4.annotate(name, (liq_delta[i], liq_n[i]), fontsize=9)
ax4.set_xlabel('Hildebrand δ (MPa^0.5)', fontsize=12)
ax4.set_ylabel('Refractive Index n', fontsize=12)
ax4.set_title(f'Liquid n vs δ\n(r = {r_liq_delta:.3f} - weak!)', fontsize=14)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/refractive_index_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved to: simulations/chemistry/refractive_index_coherence.png")

print("\n" + "=" * 70)
print("SESSION #76 COMPLETE: REFRACTIVE INDEX & COHERENCE")
print("=" * 70)
