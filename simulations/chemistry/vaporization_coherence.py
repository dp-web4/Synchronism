"""
Chemistry Session #201: Heat of Vaporization at γ ~ 1
Analyzing vaporization thermodynamics through the coherence framework.

Key hypothesis: Trouton's rule ΔS_vap ≈ 85 J/(mol·K) IS γ ~ 1
- ΔH_vap / T_b ≈ 85 J/(mol·K) for many liquids
- This represents the entropy of vaporization
- Deviations indicate H-bonding or molecular complexity

γ parameters:
- ΔS_vap / 85 (Trouton deviation)
- ΔH_vap / ΔH_sub (sublimation ratio)
- Clausius-Clapeyron slope
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 60)
print("CHEMISTRY SESSION #201: VAPORIZATION AT γ ~ 1")
print("=" * 60)

# ============================================================
# 1. TROUTON'S RULE
# ============================================================
print("\n" + "=" * 60)
print("1. TROUTON'S RULE")
print("=" * 60)

print("""
Trouton's rule:
  ΔS_vap = ΔH_vap / T_b ≈ 85 J/(mol·K)

This is remarkably universal for non-polar liquids!

γ_Trouton = ΔS_vap / 85

At γ = 1: normal vaporization behavior
γ > 1: strong H-bonding (water, alcohols)
γ < 1: dimerization in liquid (acetic acid)

The constant ~85 J/(mol·K) IS the γ ~ 1 reference.
""")

# Vaporization data: (ΔH_vap kJ/mol, T_b K)
vap_data = {
    # Non-polar - should follow Trouton
    'Ar': (6.5, 87.3),
    'N2': (5.6, 77.4),
    'O2': (6.8, 90.2),
    'CH4': (8.2, 111.7),
    'C2H6': (14.7, 184.6),
    'C3H8': (19.0, 231.0),
    'n-Hexane': (28.9, 341.9),
    'n-Octane': (34.4, 398.8),
    'Benzene': (30.8, 353.3),
    'CCl4': (30.0, 349.9),
    'Cyclohexane': (29.9, 353.9),
    # Polar - slight deviation
    'Acetone': (29.1, 329.4),
    'Diethyl ether': (26.5, 307.6),
    'Chloroform': (29.4, 334.3),
    # H-bonding - high deviation
    'Water': (40.7, 373.2),
    'Methanol': (35.2, 337.8),
    'Ethanol': (38.6, 351.4),
    'Ammonia': (23.4, 239.8),
    'HF': (25.2, 292.7),
    # Dimerization - low deviation
    'Acetic acid': (23.7, 391.1),
    'Formic acid': (22.7, 373.7),
}

print("\nTrouton's Rule Analysis:")
print("-" * 70)
print(f"{'Liquid':<15} {'ΔH_vap (kJ)':<12} {'T_b (K)':<10} {'ΔS_vap':<12} {'γ = ΔS/85'}")
print("-" * 70)

gamma_values = []
delta_S = []
categories = []

for liquid, (dH, T_b) in vap_data.items():
    dS = dH * 1000 / T_b  # Convert to J/(mol·K)
    delta_S.append(dS)
    gamma = dS / 85
    gamma_values.append(gamma)

    if liquid in ['Water', 'Methanol', 'Ethanol', 'Ammonia', 'HF']:
        cat = 'H-bond'
    elif liquid in ['Acetic acid', 'Formic acid']:
        cat = 'dimer'
    else:
        cat = 'normal'
    categories.append(cat)

    print(f"{liquid:<15} {dH:<12.1f} {T_b:<10.1f} {dS:<12.1f} {gamma:.2f}")

print("-" * 70)
mean_gamma = np.mean(gamma_values)
std_gamma = np.std(gamma_values)
print(f"Mean γ = {mean_gamma:.2f} ± {std_gamma:.2f}")
print(f"Mean ΔS_vap = {np.mean(delta_S):.1f} J/(mol·K)")

# By category
normal = [g for g, c in zip(gamma_values, categories) if c == 'normal']
h_bond = [g for g, c in zip(gamma_values, categories) if c == 'H-bond']
dimer = [g for g, c in zip(gamma_values, categories) if c == 'dimer']

print(f"\nBy category:")
print(f"  Normal liquids: γ = {np.mean(normal):.2f} ± {np.std(normal):.2f}")
print(f"  H-bonding: γ = {np.mean(h_bond):.2f} ± {np.std(h_bond):.2f} (HIGH)")
print(f"  Dimerizing: γ = {np.mean(dimer):.2f} ± {np.std(dimer):.2f} (LOW)")

# ============================================================
# 2. PICTET-TROUTON VARIATION
# ============================================================
print("\n" + "=" * 60)
print("2. PICTET-TROUTON (MODIFIED)")
print("=" * 60)

print("""
Pictet-Trouton modification:
  ΔS_vap = 85 + 4.2 × ln(T_b/298)  (temperature correction)

For low-boiling: ΔS_vap < 85
For high-boiling: ΔS_vap > 85

This accounts for temperature dependence of liquid entropy.

γ_PT = ΔS_obs / ΔS_PT_predicted
""")

print("\nPictet-Trouton Correction:")
print("-" * 60)
print(f"{'Liquid':<15} {'ΔS_obs':<12} {'ΔS_PT_pred':<12} {'γ_PT'}")
print("-" * 60)

gamma_PT = []
for liquid, (dH, T_b) in vap_data.items():
    dS_obs = dH * 1000 / T_b
    dS_PT = 85 + 4.2 * np.log(T_b / 298)
    gamma = dS_obs / dS_PT
    gamma_PT.append(gamma)
    print(f"{liquid:<15} {dS_obs:<12.1f} {dS_PT:<12.1f} {gamma:.2f}")

print("-" * 60)
print(f"Mean γ_PT = {np.mean(gamma_PT):.2f} ± {np.std(gamma_PT):.2f}")

# ============================================================
# 3. WATSON CORRELATION
# ============================================================
print("\n" + "=" * 60)
print("3. WATSON CORRELATION (ΔH_vap vs T)")
print("=" * 60)

print("""
Watson correlation for temperature dependence:
  ΔH_vap(T) = ΔH_vap(T_b) × [(T_c - T)/(T_c - T_b)]^0.38

At T = T_c: ΔH_vap = 0 (critical point)
At T = T_b: normal vaporization

γ_Watson = (T_c - T) / (T_c - T_b)

The exponent 0.38 is nearly universal!
""")

# Critical temperatures
T_c_data = {
    'Water': 647.3,
    'Methanol': 512.6,
    'Ethanol': 514.0,
    'Benzene': 562.2,
    'n-Hexane': 507.6,
    'Ammonia': 405.5,
}

print("\nWatson Correlation Check:")
print("-" * 50)
print(f"{'Liquid':<15} {'T_b (K)':<10} {'T_c (K)':<10} {'T_b/T_c'}")
print("-" * 50)

T_ratio_list = []
for liquid in T_c_data:
    if liquid in vap_data:
        T_b = vap_data[liquid][1]
        T_c = T_c_data[liquid]
        ratio = T_b / T_c
        T_ratio_list.append(ratio)
        print(f"{liquid:<15} {T_b:<10.1f} {T_c:<10.1f} {ratio:.3f}")

print("-" * 50)
print(f"Mean T_b/T_c = {np.mean(T_ratio_list):.3f} ± {np.std(T_ratio_list):.3f}")

# ============================================================
# 4. CLAUSIUS-CLAPEYRON
# ============================================================
print("\n" + "=" * 60)
print("4. CLAUSIUS-CLAPEYRON EQUATION")
print("=" * 60)

print("""
Clausius-Clapeyron:
  d(ln P) / d(1/T) = -ΔH_vap / R

At the boiling point (P = 1 atm):
  γ_CC = ΔH_vap / (R × T_b)

This is Trouton's rule in different form!
  γ_CC = ΔS_vap / R ≈ 85/8.314 ≈ 10.2

At γ ~ 10: normal vaporization
""")

R = 8.314  # J/(mol·K)

print("\nClausius-Clapeyron γ = ΔH/(RT_b):")
print("-" * 50)
gamma_CC = []
for liquid, (dH, T_b) in vap_data.items():
    gamma = (dH * 1000) / (R * T_b)
    gamma_CC.append(gamma)
    if liquid in ['Water', 'Methanol', 'Benzene', 'n-Hexane', 'Acetic acid']:
        print(f"{liquid:<15} γ_CC = {gamma:.1f}")

print(f"\nMean γ_CC = {np.mean(gamma_CC):.1f} ± {np.std(gamma_CC):.1f}")
print(f"Trouton predicts: γ_CC = 85/8.314 = {85/8.314:.1f}")

# ============================================================
# 5. ENTHALPY-ENTROPY COMPENSATION
# ============================================================
print("\n" + "=" * 60)
print("5. ENTHALPY-ENTROPY COMPENSATION")
print("=" * 60)

print("""
For vaporization across substances:
  ΔH_vap vs T_b × ΔS_vap should correlate

At equilibrium: ΔG = ΔH - TΔS = 0
So: ΔH_vap = T_b × ΔS_vap (exactly at T_b!)

This IS γ ~ 1 equilibrium.
""")

# Calculate correlation
T_b_values = [vap_data[l][1] for l in vap_data]
dH_values = [vap_data[l][0] for l in vap_data]

# ΔH should equal T_b × ΔS at boiling point
dH_predicted = [T_b * dS / 1000 for T_b, dS in zip(T_b_values, delta_S)]

correlation = np.corrcoef(dH_values, dH_predicted)[0, 1]
print(f"\nCorrelation ΔH_obs vs T_b×ΔS: r = {correlation:.4f}")
print("(Should be 1.0000 by definition at T_b!)")

# ============================================================
# 6. VISUALIZATION
# ============================================================
print("\n" + "=" * 60)
print("6. GENERATING VISUALIZATION")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Trouton γ distribution
ax1 = axes[0, 0]
colors = ['green' if c == 'normal' else ('red' if c == 'H-bond' else 'blue')
          for c in categories]
ax1.bar(range(len(gamma_values)), gamma_values, color=colors, edgecolor='black', alpha=0.7)
ax1.axhline(y=1, color='black', linestyle='--', linewidth=2, label='γ = 1 (Trouton)')
ax1.set_xticks(range(0, len(gamma_values), 3))
ax1.set_xticklabels([list(vap_data.keys())[i][:6] for i in range(0, len(gamma_values), 3)],
                    rotation=45, ha='right', fontsize=8)
ax1.set_ylabel('γ = ΔS_vap / 85', fontsize=12)
ax1.set_title("Trouton's Rule: ΔS_vap ≈ 85 J/(mol·K)", fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Add legend for colors
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='green', label='Normal'),
                   Patch(facecolor='red', label='H-bonding'),
                   Patch(facecolor='blue', label='Dimerizing')]
ax1.legend(handles=legend_elements, loc='upper right')

# Plot 2: ΔS_vap histogram
ax2 = axes[0, 1]
ax2.hist(delta_S, bins=10, edgecolor='black', alpha=0.7, color='coral')
ax2.axvline(x=85, color='red', linestyle='--', linewidth=2, label='Trouton (85)')
ax2.axvline(x=np.mean(delta_S), color='green', linestyle='-', linewidth=2,
            label=f'Mean = {np.mean(delta_S):.1f}')
ax2.set_xlabel('ΔS_vap [J/(mol·K)]', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.set_title('Entropy of Vaporization Distribution', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: ΔH vs T_b
ax3 = axes[1, 0]
ax3.scatter(T_b_values, dH_values, s=60, c=colors, edgecolors='black', alpha=0.7)
# Trouton prediction line
T_fit = np.linspace(50, 450, 100)
dH_fit = 85 * T_fit / 1000  # Trouton prediction
ax3.plot(T_fit, dH_fit, 'r--', linewidth=2, label='Trouton: ΔH = 85×T_b')
ax3.set_xlabel('Boiling Point T_b (K)', fontsize=12)
ax3.set_ylabel('ΔH_vap (kJ/mol)', fontsize=12)
ax3.set_title('ΔH_vap vs T_b: Trouton Correlation', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Pictet-Trouton γ by category
ax4 = axes[1, 1]
cat_names = ['Normal', 'H-bonding', 'Dimerizing']
cat_means = [np.mean(normal), np.mean(h_bond), np.mean(dimer)]
cat_stds = [np.std(normal), np.std(h_bond), np.std(dimer)]
cat_colors = ['green', 'red', 'blue']

bars = ax4.bar(cat_names, cat_means, yerr=cat_stds, capsize=5,
               color=cat_colors, edgecolor='black', alpha=0.7)
ax4.axhline(y=1, color='black', linestyle='--', linewidth=2, label='γ = 1')
ax4.set_ylabel('Mean γ = ΔS_vap / 85', fontsize=12)
ax4.set_title('Trouton Deviation by Category', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_ylim(0, 1.5)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/vaporization_coherence.png', dpi=150)
print("Saved: vaporization_coherence.png")
plt.close()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SESSION #201 SUMMARY: VAPORIZATION AT γ ~ 1")
print("=" * 60)

print(f"""
KEY FINDINGS:

1. TROUTON'S RULE
   ΔS_vap ≈ 85 J/(mol·K) for normal liquids
   γ_Trouton = ΔS_vap / 85
   Mean γ = {mean_gamma:.2f} ± {std_gamma:.2f}

2. CATEGORY DEVIATIONS
   Normal liquids: γ = {np.mean(normal):.2f} (γ ~ 1!)
   H-bonding: γ = {np.mean(h_bond):.2f} (higher order in liquid)
   Dimerizing: γ = {np.mean(dimer):.2f} (partial order in vapor)

3. CLAUSIUS-CLAPEYRON
   γ_CC = ΔH/(RT_b) = ΔS/R ≈ 10.2
   Mean γ_CC = {np.mean(gamma_CC):.1f}
   Trouton IS Clausius-Clapeyron at T_b

4. BOILING POINT RATIO
   Mean T_b/T_c = {np.mean(T_ratio_list):.3f}
   Boiling point ≈ 60% of critical temperature

5. ENTHALPY-ENTROPY COMPENSATION
   ΔG = 0 at equilibrium (T = T_b)
   ΔH = T×ΔS exactly at boiling point
   r = {correlation:.4f} (perfect by definition)

CENTRAL INSIGHT:
Trouton's rule ΔS_vap ≈ 85 J/(mol·K) IS γ ~ 1:
- Normal liquids have ΔS_vap ≈ 85 (γ = 1)
- H-bonding increases order → higher ΔS_vap
- Dimerization carries order to vapor → lower ΔS_vap

The universal entropy of vaporization reflects
the coherence difference between liquid (ordered)
and gas (disordered) phases.

This is the 64th phenomenon type at γ ~ 1!
""")

print("=" * 60)
print("SESSION #201 COMPLETE")
print("=" * 60)
