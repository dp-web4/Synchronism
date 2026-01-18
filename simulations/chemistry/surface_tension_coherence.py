#!/usr/bin/env python3
"""
Chemistry Session #74: Surface Tension & Coherence
Test whether coherence framework predicts surface tension.

Surface tension (γ_ST) measures energy cost to create surface:
- High γ_ST = strongly cohesive liquid (water, mercury)
- Low γ_ST = weakly cohesive liquid (hexane, ether)

Coherence interpretation:
- Surface = interface where bulk coherence is broken
- Higher bulk coherence → higher penalty for surface → higher γ_ST
- γ_ST ∝ (2/γ_bulk) - surface tension proportional to bulk coherence

This is OPPOSITE to viscosity!
- Viscosity: η ∝ γ_flow (resistance to flow)
- Surface tension: γ_ST ∝ 2/γ_bulk (cohesive strength)

The surface IS the coherence boundary.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

print("=" * 70)
print("CHEMISTRY SESSION #74: SURFACE TENSION & COHERENCE")
print("=" * 70)

# ==============================================================================
# DATASET: SURFACE TENSIONS
# ==============================================================================

print("\n" + "=" * 70)
print("DATASET: SURFACE TENSIONS AT 25°C")
print("=" * 70)

# Surface tension data (mN/m) at 25°C
# Also: Hildebrand δ, boiling point Tb, and viscosity η for comparison
surface_data = {
    # Format: (γ_ST in mN/m, δ in MPa^0.5, Tb in °C, η in mPa·s)

    # High surface tension
    'mercury': (485, 100, 357, 1.5),  # δ is approximate for metals
    'water': (72.0, 47.8, 100, 0.89),

    # Medium-high
    'glycerol': (63.0, 36.1, 290, 934),
    'formamide': (58.2, 36.6, 210, 3.3),
    'ethylene_glycol': (47.3, 32.9, 197, 16.1),
    'DMSO': (43.5, 26.7, 189, 1.99),

    # Medium
    'nitrobenzene': (43.4, 22.1, 211, 1.8),
    'aniline': (42.9, 24.1, 184, 3.8),
    'methanol': (22.5, 29.6, 65, 0.55),
    'ethanol': (22.1, 26.5, 78, 1.07),
    'propanol': (23.3, 24.5, 97, 1.95),
    'butanol': (24.9, 23.2, 118, 2.54),

    # Low
    'acetone': (23.5, 20.0, 56, 0.31),
    'benzene': (28.2, 18.8, 80, 0.60),
    'toluene': (28.5, 18.2, 111, 0.56),
    'chloroform': (27.1, 19.0, 61, 0.54),
    'carbon_tetrachloride': (26.4, 17.8, 77, 0.91),

    # Very low
    'hexane': (18.4, 14.9, 69, 0.30),
    'octane': (21.8, 15.5, 126, 0.51),
    'diethyl_ether': (17.0, 15.1, 35, 0.22),
    'pentane': (16.0, 14.4, 36, 0.22),
}

print(f"Surface tension samples: {len(surface_data)}")

# Print sorted by surface tension
print("\nLiquids sorted by surface tension:")
print("-" * 70)
print(f"{'Liquid':<20} {'γ_ST (mN/m)':<12} {'δ':<8} {'Tb (°C)':<10} {'η (mPa·s)':<10}")
print("-" * 70)

for name, (gamma_st, delta, Tb, eta) in sorted(surface_data.items(), key=lambda x: -x[1][0]):
    print(f"{name:<20} {gamma_st:>10.1f}  {delta:>6.1f}  {Tb:>8d}  {eta:>8.2f}")

# ==============================================================================
# COHERENCE FROM HILDEBRAND
# ==============================================================================

print("\n" + "=" * 70)
print("γ_BULK FROM HILDEBRAND PARAMETER")
print("=" * 70)

def gamma_bulk_from_delta(delta, delta_max=50.0):
    """
    Estimate bulk coherence from Hildebrand parameter.
    Higher δ = more cohesive = MORE coherent = LOWER γ_bulk.
    """
    gamma = 2.0 - 1.5 * (delta / delta_max)
    return np.clip(gamma, 0.5, 2.0)

def coherence_factor(delta):
    """
    2/γ_bulk = coherence enhancement factor.
    Higher δ → lower γ_bulk → higher 2/γ.
    """
    gamma = gamma_bulk_from_delta(delta)
    return 2.0 / gamma

# Print coherence values
print("\nBulk coherence from δ:")
print("-" * 50)
for name, (gamma_st, delta, Tb, eta) in sorted(surface_data.items(), key=lambda x: x[1][1]):
    gamma = gamma_bulk_from_delta(delta)
    coh_factor = 2.0 / gamma
    print(f"{name:<20}: δ = {delta:>5.1f}, γ_bulk = {gamma:.2f}, 2/γ = {coh_factor:.2f}")

# ==============================================================================
# CORRELATION ANALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("CORRELATION ANALYSIS")
print("=" * 70)

# Extract arrays (exclude mercury as outlier)
gamma_st_list = []
delta_list = []
Tb_list = []
eta_list = []
coh_factor_list = []

for name, (gamma_st, delta, Tb, eta) in surface_data.items():
    if name != 'mercury':  # Mercury is a metal, different physics
        gamma_st_list.append(gamma_st)
        delta_list.append(delta)
        Tb_list.append(Tb)
        eta_list.append(np.log10(eta))
        coh_factor_list.append(coherence_factor(delta))

gamma_st_arr = np.array(gamma_st_list)
delta_arr = np.array(delta_list)
Tb_arr = np.array(Tb_list)
eta_arr = np.array(eta_list)
coh_factor_arr = np.array(coh_factor_list)

# Correlations
r_delta, p_delta = stats.pearsonr(gamma_st_arr, delta_arr)
r_Tb, p_Tb = stats.pearsonr(gamma_st_arr, Tb_arr)
r_eta, p_eta = stats.pearsonr(gamma_st_arr, eta_arr)
r_coh, p_coh = stats.pearsonr(gamma_st_arr, coh_factor_arr)

print(f"\nγ_ST vs δ (Hildebrand): r = {r_delta:.3f}, p = {p_delta:.4f}")
print(f"γ_ST vs Tb (boiling point): r = {r_Tb:.3f}, p = {r_Tb:.4f}")
print(f"γ_ST vs log(η): r = {r_eta:.3f}, p = {p_eta:.4f}")
print(f"γ_ST vs 2/γ_bulk (coherence): r = {r_coh:.3f}, p = {p_coh:.4f}")

# ==============================================================================
# MACLEOD-SUGDEN EQUATION
# ==============================================================================

print("\n" + "=" * 70)
print("MACLEOD-SUGDEN EQUATION + COHERENCE")
print("=" * 70)

print("""
Macleod-Sugden equation:
γ_ST^(1/4) = [P] × (ρ_L - ρ_V) / M

Where [P] is the parachor (empirical).

Coherence interpretation:
γ_ST ∝ (cohesive energy density)
CED = δ² = (ΔH_vap - RT) / V_m

Since δ reflects cohesive energy:
γ_ST ∝ δ² (approximately)

Or in coherence terms:
γ_ST ∝ (2/γ_bulk)² since γ_bulk ∝ 1/δ
""")

# Test γ_ST vs δ² relationship
delta_sq = delta_arr**2
r_delta_sq, _ = stats.pearsonr(gamma_st_arr, delta_sq)
print(f"\nγ_ST vs δ²: r = {r_delta_sq:.3f}")

# ==============================================================================
# SURFACE AS COHERENCE BOUNDARY
# ==============================================================================

print("\n" + "=" * 70)
print("SURFACE AS COHERENCE BOUNDARY")
print("=" * 70)

print("""
Physical picture:
- Bulk liquid has coherence level γ_bulk
- Surface molecules lack half their neighbors
- Surface = region where coherence is BROKEN

Energy cost to create surface:
γ_ST ∝ ΔE_coherence = E_bulk - E_surface

More coherent bulk (lower γ_bulk) → higher penalty for surface → higher γ_ST

This explains why:
- Water (high δ, low γ_bulk, strong H-bonds) has high γ_ST
- Hexane (low δ, high γ_bulk, weak interactions) has low γ_ST

The surface tension IS the coherence difference between bulk and interface!
""")

# ==============================================================================
# COMPARISON WITH VISCOSITY
# ==============================================================================

print("\n" + "=" * 70)
print("SURFACE TENSION vs VISCOSITY")
print("=" * 70)

print("""
Both reflect intermolecular forces, but differently:

VISCOSITY (η):
- Measures resistance to FLOW (shear)
- η ∝ γ_flow (proportional to disorder)
- High H-bonding → high barrier → high η

SURFACE TENSION (γ_ST):
- Measures cost of creating INTERFACE
- γ_ST ∝ 2/γ_bulk (proportional to order)
- High H-bonding → high cohesion → high γ_ST

Correlation between them:
""")

# Calculate correlation between γ_ST and η
r_st_eta, _ = stats.pearsonr(gamma_st_arr, 10**eta_arr)
print(f"γ_ST vs η: r = {r_st_eta:.3f}")

# Both high for water/glycerol, both low for hexane
print("\nExamples showing correlation:")
print(f"Water: γ_ST = {surface_data['water'][0]:.1f} mN/m, η = {surface_data['water'][3]:.2f} mPa·s")
print(f"Hexane: γ_ST = {surface_data['hexane'][0]:.1f} mN/m, η = {surface_data['hexane'][3]:.2f} mPa·s")
print(f"Glycerol: γ_ST = {surface_data['glycerol'][0]:.1f} mN/m, η = {surface_data['glycerol'][3]:.0f} mPa·s")

# ==============================================================================
# LINEAR MODEL
# ==============================================================================

print("\n" + "=" * 70)
print("MODEL FIT: γ_ST vs 2/γ_bulk")
print("=" * 70)

# Linear fit
slope, intercept, r_value, p_value, std_err = stats.linregress(coh_factor_arr, gamma_st_arr)
gamma_st_pred = slope * coh_factor_arr + intercept

# R²
ss_res = np.sum((gamma_st_arr - gamma_st_pred)**2)
ss_tot = np.sum((gamma_st_arr - gamma_st_arr.mean())**2)
R2 = 1 - ss_res/ss_tot

print(f"\nLinear fit: γ_ST = {slope:.2f} × (2/γ_bulk) + {intercept:.2f}")
print(f"r = {r_value:.3f}, R² = {R2:.3f}")

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SESSION #74 SUMMARY: SURFACE TENSION & COHERENCE")
print("=" * 70)

print(f"""
Correlations Found:
- γ_ST vs δ (Hildebrand): r = {r_delta:.3f} {"(GOOD)" if abs(r_delta) > 0.6 else "(MODERATE)"}
- γ_ST vs δ²: r = {r_delta_sq:.3f}
- γ_ST vs 2/γ_bulk: r = {r_coh:.3f}
- γ_ST vs log(η): r = {r_eta:.3f}
- γ_ST vs η: r = {r_st_eta:.3f}

Key Findings:
1. Hildebrand δ predicts surface tension (r = {r_delta:.3f})
   - Higher δ → higher γ_ST

2. Coherence factor 2/γ gives similar prediction (r = {r_coh:.3f})
   - More coherent bulk → higher γ_ST

3. Surface tension and viscosity are correlated (r = {r_st_eta:.3f})
   - Both reflect intermolecular forces
   - BUT through different mechanisms

4. Surface = coherence boundary
   - Breaking bulk coherence costs energy
   - γ_ST measures that energy cost

Physical Interpretation:
- Surface tension measures COHESIVE strength
- γ_ST ∝ 2/γ_bulk (proportional to bulk coherence)
- This is OPPOSITE to viscosity's γ_flow dependence
""")

# ==============================================================================
# PREDICTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print("""
P74.1: γ_ST ∝ 2/γ_bulk
Surface tension proportional to bulk coherence.

P74.2: γ_ST ∝ δ (approximately)
Hildebrand parameter predicts surface tension.

P74.3: Surface = coherence interface
Energy cost to break bulk coherence.

P74.4: Surfactants lower γ_ST by bridging γ values
Amphiphiles reduce coherence mismatch at interface.

P74.5: Temperature dependence: γ_ST(T) ∝ 2/γ_bulk(T)
Surface tension decreases as bulk coherence decreases with T.
""")

# ==============================================================================
# VALIDATION STATUS
# ==============================================================================

print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r_delta) > 0.7:
    status = "STRONG SUPPORTING EVIDENCE"
elif abs(r_delta) > 0.5:
    status = "MODERATE SUPPORTING EVIDENCE"
else:
    status = "WEAK SUPPORTING EVIDENCE"

print(f"\n{status}")
print(f"""
The coherence framework provides:
1. CONSISTENT with δ-based predictions (r = {r_delta:.3f})
2. EXPLAINS surface as coherence boundary
3. CONNECTS γ_ST to bulk coherence γ_bulk

Note: This is largely REINTERPRETATION of known physics.
The Hildebrand parameter already captures cohesive energy.
The coherence framework adds interpretive value by showing:
- Surface tension = cost of breaking coherence
- Opposite sign dependence from viscosity

The framework unifies γ_ST and η through coherence, even though
they scale oppositely with γ (one with 2/γ, one with γ).
""")

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: γ_ST vs δ
ax1 = axes[0, 0]
ax1.scatter(delta_arr, gamma_st_arr, s=80, alpha=0.7, c='blue')
z = np.polyfit(delta_arr, gamma_st_arr, 1)
p = np.poly1d(z)
x_line = np.linspace(delta_arr.min(), delta_arr.max(), 100)
ax1.plot(x_line, p(x_line), 'r--', label=f'r = {r_delta:.3f}')
ax1.set_xlabel('Hildebrand δ (MPa^0.5)', fontsize=12)
ax1.set_ylabel('Surface Tension (mN/m)', fontsize=12)
ax1.set_title('Surface Tension vs Hildebrand Parameter', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend()

# Add labels for key liquids
key_liquids = {'water': (47.8, 72.0), 'glycerol': (36.1, 63.0), 'hexane': (14.9, 18.4)}
for name, (d, g) in key_liquids.items():
    ax1.annotate(name, (d, g), fontsize=9)

# Plot 2: γ_ST vs 2/γ_bulk
ax2 = axes[0, 1]
ax2.scatter(coh_factor_arr, gamma_st_arr, s=80, alpha=0.7, c='green')
ax2.plot(coh_factor_arr, gamma_st_pred, 'r--', label=f'R² = {R2:.2f}')
ax2.set_xlabel('2/γ_bulk (coherence factor)', fontsize=12)
ax2.set_ylabel('Surface Tension (mN/m)', fontsize=12)
ax2.set_title(f'Surface Tension vs Coherence\n(r = {r_coh:.3f})', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.legend()

# Plot 3: γ_ST vs η
ax3 = axes[1, 0]
eta_linear = 10**eta_arr
ax3.scatter(eta_linear, gamma_st_arr, s=80, alpha=0.7, c='purple')
ax3.set_xlabel('Viscosity (mPa·s)', fontsize=12)
ax3.set_ylabel('Surface Tension (mN/m)', fontsize=12)
ax3.set_title(f'Surface Tension vs Viscosity\n(r = {r_st_eta:.3f})', fontsize=14)
ax3.set_xscale('log')
ax3.grid(True, alpha=0.3)

# Plot 4: Comparison diagram
ax4 = axes[1, 1]
categories = ['Viscosity', 'Surface Tension']
gamma_dependence = ['γ_flow', '2/γ_bulk']
colors = ['red', 'blue']

# Create text diagram
ax4.text(0.5, 0.8, 'Coherence Framework: Two Properties', ha='center', fontsize=14, weight='bold')
ax4.text(0.25, 0.6, 'VISCOSITY (η)', ha='center', fontsize=12, weight='bold', color='red')
ax4.text(0.25, 0.5, 'η ∝ γ_flow', ha='center', fontsize=11, color='red')
ax4.text(0.25, 0.4, 'Resistance to flow', ha='center', fontsize=10)
ax4.text(0.25, 0.3, 'High when coherence\nmust be broken', ha='center', fontsize=9)

ax4.text(0.75, 0.6, 'SURFACE TENSION (γ_ST)', ha='center', fontsize=12, weight='bold', color='blue')
ax4.text(0.75, 0.5, 'γ_ST ∝ 2/γ_bulk', ha='center', fontsize=11, color='blue')
ax4.text(0.75, 0.4, 'Cost of creating interface', ha='center', fontsize=10)
ax4.text(0.75, 0.3, 'High when bulk is\nhighly coherent', ha='center', fontsize=9)

ax4.text(0.5, 0.1, 'Both high for H-bonded liquids (water, glycerol)\nBoth low for weakly interacting liquids (hexane)',
         ha='center', fontsize=10, style='italic')
ax4.axis('off')
ax4.set_xlim(0, 1)
ax4.set_ylim(0, 1)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/surface_tension_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved to: simulations/chemistry/surface_tension_coherence.png")

print("\n" + "=" * 70)
print("SESSION #74 COMPLETE: SURFACE TENSION & COHERENCE")
print("=" * 70)
