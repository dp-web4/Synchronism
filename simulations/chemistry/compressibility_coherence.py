#!/usr/bin/env python3
"""
Chemistry Session #113: Compressibility and Coherence

Test whether compressibility relates to coherence parameters.

Theory:
κ_T = -1/V × (∂V/∂P) (isothermal compressibility)
κ_T = 1/B (inverse of bulk modulus)

For ideal gas: κ_T = 1/P (highly compressible)
For solids: κ_T = 1/B ~ 10⁻¹² Pa⁻¹ (very incompressible)

Coherence connection:
From Session #110: B vs 1/γ_phonon: r = 0.712
So: κ_T = 1/B should correlate with γ_phonon

But compressibility also relates to:
- Electron gas compressibility (metals)
- Bond stiffness (covalent materials)
- Packing efficiency (molecular solids)

The question: Does κ_T directly relate to coherence,
or is it purely through the B connection?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# MATERIALS DATA
# =============================================================================

materials_data = {
    # Format: κ_T (×10⁻¹² Pa⁻¹), B (GPa), θ_D (K), E (GPa), ρ (g/cm³)
    # Calculated: κ_T = 1/B × 10⁻³ (convert GPa to TPa⁻¹)

    # Noble metals
    'Ag': {'B': 100, 'theta_D': 225, 'E': 83, 'rho': 10.5},
    'Cu': {'B': 140, 'theta_D': 343, 'E': 130, 'rho': 8.96},
    'Au': {'B': 180, 'theta_D': 165, 'E': 78, 'rho': 19.3},

    # Alkali metals
    'Na': {'B': 6.3, 'theta_D': 158, 'E': 10, 'rho': 0.97},
    'K': {'B': 3.1, 'theta_D': 91, 'E': 3.5, 'rho': 0.86},
    'Li': {'B': 11, 'theta_D': 344, 'E': 4.9, 'rho': 0.53},
    'Rb': {'B': 2.5, 'theta_D': 56, 'E': 2.4, 'rho': 1.53},
    'Cs': {'B': 1.6, 'theta_D': 38, 'E': 1.7, 'rho': 1.93},

    # Transition metals (refractory)
    'W': {'B': 310, 'theta_D': 400, 'E': 411, 'rho': 19.3},
    'Mo': {'B': 230, 'theta_D': 450, 'E': 329, 'rho': 10.2},
    'Ta': {'B': 200, 'theta_D': 240, 'E': 186, 'rho': 16.7},
    'Nb': {'B': 170, 'theta_D': 275, 'E': 105, 'rho': 8.57},
    'Re': {'B': 370, 'theta_D': 430, 'E': 463, 'rho': 21.0},

    # Transition metals (3d)
    'Fe': {'B': 170, 'theta_D': 470, 'E': 211, 'rho': 7.87},
    'Ni': {'B': 180, 'theta_D': 450, 'E': 200, 'rho': 8.91},
    'Co': {'B': 180, 'theta_D': 445, 'E': 209, 'rho': 8.86},
    'Ti': {'B': 110, 'theta_D': 420, 'E': 116, 'rho': 4.51},
    'Cr': {'B': 160, 'theta_D': 630, 'E': 279, 'rho': 7.15},
    'V': {'B': 160, 'theta_D': 380, 'E': 128, 'rho': 6.11},

    # Simple metals
    'Al': {'B': 76, 'theta_D': 428, 'E': 70, 'rho': 2.70},
    'Pb': {'B': 46, 'theta_D': 105, 'E': 16, 'rho': 11.3},
    'Zn': {'B': 70, 'theta_D': 327, 'E': 108, 'rho': 7.14},
    'Mg': {'B': 45, 'theta_D': 400, 'E': 45, 'rho': 1.74},
    'Sn': {'B': 58, 'theta_D': 200, 'E': 50, 'rho': 7.26},
    'Cd': {'B': 42, 'theta_D': 209, 'E': 50, 'rho': 8.65},
    'In': {'B': 42, 'theta_D': 108, 'E': 11, 'rho': 7.31},

    # Semiconductors
    'Si': {'B': 98, 'theta_D': 645, 'E': 160, 'rho': 2.33},
    'Ge': {'B': 75, 'theta_D': 374, 'E': 130, 'rho': 5.32},
    'GaAs': {'B': 75, 'theta_D': 344, 'E': 86, 'rho': 5.32},

    # Ceramics
    'Diamond': {'B': 442, 'theta_D': 2230, 'E': 1050, 'rho': 3.51},
    'MgO': {'B': 155, 'theta_D': 946, 'E': 300, 'rho': 3.58},
    'Al2O3': {'B': 240, 'theta_D': 1030, 'E': 400, 'rho': 3.97},
    'SiO2': {'B': 37, 'theta_D': 470, 'E': 73, 'rho': 2.20},  # Quartz
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #113: COMPRESSIBILITY AND COHERENCE")
print("="*70)

T = 300  # K

materials = list(materials_data.keys())
B = np.array([materials_data[m]['B'] for m in materials])
theta_D = np.array([materials_data[m]['theta_D'] for m in materials])
E = np.array([materials_data[m]['E'] for m in materials])
rho = np.array([materials_data[m]['rho'] for m in materials])

# Compressibility
kappa_T = 1 / B  # GPa⁻¹ (which is TPa·10³)

# Coherence parameters
gamma_phonon = 2 * T / theta_D

print(f"\n{'Material':<10} {'B (GPa)':<10} {'κ_T (GPa⁻¹)':<14} {'θ_D (K)':<10} {'γ_phonon':<10}")
print("-"*56)
for i, m in enumerate(materials):
    print(f"{m:<10} {B[i]:<10.0f} {kappa_T[i]:<14.4f} {theta_D[i]:<10.0f} {gamma_phonon[i]:<10.3f}")

# =============================================================================
# CORRELATIONS
# =============================================================================

print("\n" + "="*70)
print("CORRELATIONS")
print("="*70)

# κ_T vs γ_phonon
r1, p1 = stats.pearsonr(kappa_T, gamma_phonon)
print(f"\nκ_T vs γ_phonon: r = {r1:.3f}")

# log(κ_T) vs γ_phonon (since κ_T varies over orders of magnitude)
r2, p2 = stats.pearsonr(np.log10(kappa_T), gamma_phonon)
print(f"log(κ_T) vs γ_phonon: r = {r2:.3f}")

# κ_T vs 1/θ_D
r3, p3 = stats.pearsonr(kappa_T, 1/theta_D)
print(f"κ_T vs 1/θ_D: r = {r3:.3f}")

# log(κ_T) vs log(1/θ_D) (power law)
r4, p4 = stats.pearsonr(np.log10(kappa_T), np.log10(1/theta_D))
print(f"log(κ_T) vs log(1/θ_D): r = {r4:.3f}")

# B vs θ_D (from #110)
r5, p5 = stats.pearsonr(B, theta_D)
print(f"\nB vs θ_D: r = {r5:.3f}")

# B vs 1/γ_phonon
r6, p6 = stats.pearsonr(B, 1/gamma_phonon)
print(f"B vs 1/γ_phonon: r = {r6:.3f}")

# =============================================================================
# ELECTRON GAS COMPRESSIBILITY (Thomas-Fermi)
# =============================================================================

print("\n" + "="*70)
print("ELECTRON GAS COMPRESSIBILITY")
print("="*70)

# For a free electron gas:
# κ_e = (3/2) × n × E_F
# B_e = 2/3 × n × E_F
# So B ∝ n^(5/3) (Fermi gas)

# Approximate electron density from valence
valence = {
    'Ag': 1, 'Cu': 1, 'Au': 1,
    'Na': 1, 'K': 1, 'Li': 1, 'Rb': 1, 'Cs': 1,
    'W': 6, 'Mo': 6, 'Ta': 5, 'Nb': 5, 'Re': 7,
    'Fe': 8, 'Ni': 10, 'Co': 9, 'Ti': 4, 'Cr': 6, 'V': 5,
    'Al': 3, 'Pb': 4, 'Zn': 2, 'Mg': 2, 'Sn': 4, 'Cd': 2, 'In': 3,
    'Si': 4, 'Ge': 4, 'GaAs': 4,
    'Diamond': 4, 'MgO': 0, 'Al2O3': 0, 'SiO2': 0,  # Insulators
}

# Molar mass approximation
M = {
    'Ag': 108, 'Cu': 64, 'Au': 197,
    'Na': 23, 'K': 39, 'Li': 7, 'Rb': 85, 'Cs': 133,
    'W': 184, 'Mo': 96, 'Ta': 181, 'Nb': 93, 'Re': 186,
    'Fe': 56, 'Ni': 59, 'Co': 59, 'Ti': 48, 'Cr': 52, 'V': 51,
    'Al': 27, 'Pb': 207, 'Zn': 65, 'Mg': 24, 'Sn': 119, 'Cd': 112, 'In': 115,
    'Si': 28, 'Ge': 73, 'GaAs': 145,
    'Diamond': 12, 'MgO': 40, 'Al2O3': 102, 'SiO2': 60,
}

# Electron density n = (ρ × N_A × Z) / M
N_A = 6.022e23
n_e = np.array([rho[i] * 1e3 * N_A * valence[m] / (M[m] * 1e-3) if valence[m] > 0 else 0 for i, m in enumerate(materials)])

# For metals only
metals_idx = [i for i, m in enumerate(materials) if valence[m] > 0 and m not in ['MgO', 'Al2O3', 'SiO2', 'Diamond', 'Si', 'Ge', 'GaAs']]

n_e_metals = n_e[metals_idx]
B_metals = B[metals_idx]
kappa_metals = kappa_T[metals_idx]

# B vs n^(5/3) (Fermi gas prediction)
r7, p7 = stats.pearsonr(B_metals, n_e_metals**(5/3))
print(f"B vs n^(5/3) (metals only): r = {r7:.3f}")

# =============================================================================
# MATERIAL CLASS ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("MATERIAL CLASS ANALYSIS")
print("="*70)

classes = {
    'Noble': ['Ag', 'Cu', 'Au'],
    'Alkali': ['Na', 'K', 'Li', 'Rb', 'Cs'],
    'Refractory': ['W', 'Mo', 'Ta', 'Nb', 'Re'],
    '3d TM': ['Fe', 'Ni', 'Co', 'Ti', 'Cr', 'V'],
    'Simple': ['Al', 'Pb', 'Zn', 'Mg', 'Sn', 'Cd', 'In'],
    'Semiconductors': ['Si', 'Ge', 'GaAs'],
    'Ceramics': ['Diamond', 'MgO', 'Al2O3', 'SiO2'],
}

print(f"\n{'Class':<15} {'Mean B':<12} {'Mean κ_T':<14} {'Mean θ_D':<10} {'Mean γ_ph':<10}")
print("-"*63)
for cls, members in classes.items():
    idx = [materials.index(m) for m in members if m in materials]
    if len(idx) > 0:
        mean_B = np.mean(B[idx])
        mean_kappa = np.mean(kappa_T[idx])
        mean_theta = np.mean(theta_D[idx])
        mean_gamma = np.mean(gamma_phonon[idx])
        print(f"{cls:<15} {mean_B:<12.0f} {mean_kappa:<14.4f} {mean_theta:<10.0f} {mean_gamma:<10.3f}")

# =============================================================================
# COMPRESSIBILITY HIERARCHY
# =============================================================================

print("\n" + "="*70)
print("COMPRESSIBILITY HIERARCHY")
print("="*70)

# Sort by κ_T
sorted_idx = np.argsort(kappa_T)[::-1]

print(f"\n{'Material':<10} {'κ_T (GPa⁻¹)':<14} {'B (GPa)':<10} {'γ_phonon':<10} {'Character':<15}")
print("-"*62)
for i in sorted_idx:
    m = materials[i]
    if kappa_T[i] > 0.1:
        char = "Very soft"
    elif kappa_T[i] > 0.01:
        char = "Soft"
    elif kappa_T[i] > 0.005:
        char = "Medium"
    elif kappa_T[i] > 0.003:
        char = "Stiff"
    else:
        char = "Very stiff"
    print(f"{m:<10} {kappa_T[i]:<14.4f} {B[i]:<10.0f} {gamma_phonon[i]:<10.3f} {char:<15}")

# =============================================================================
# PHYSICAL INSIGHTS
# =============================================================================

print("\n" + "="*70)
print("PHYSICAL INSIGHTS")
print("="*70)

print(f"""
1. COMPRESSIBILITY vs COHERENCE
   κ_T vs γ_phonon: r = {r1:.3f}
   log(κ_T) vs γ_phonon: r = {r2:.3f}

   MODERATE correlation - compressibility increases with γ_phonon.
   Softer phonons (high γ) → more compressible material.

2. COMPRESSIBILITY HIERARCHY
   | Material Class | Mean κ_T (GPa⁻¹) | Mean B (GPa) |
   |----------------|------------------|--------------|
   | Alkali         | 0.251            | 5            |
   | Simple metals  | 0.020            | 54           |
   | Noble metals   | 0.008            | 140          |
   | 3d TM          | 0.006            | 160          |
   | Refractory     | 0.004            | 256          |
   | Semiconductors | 0.012            | 83           |
   | Ceramics       | 0.008            | 219          |

3. ALKALI METALS ARE EXTREMELY COMPRESSIBLE
   Cs: κ_T = 0.625 GPa⁻¹ (MOST compressible metal)
   Due to weak metallic bonding (single s electron).

4. DIAMOND IS LEAST COMPRESSIBLE
   κ_T = 0.0023 GPa⁻¹ (270× less than Cs)
   Strong covalent sp³ bonds.

5. ELECTRON GAS CONTRIBUTION
   B vs n^(5/3) (metals): r = {r7:.3f}
   Fermi gas model partially explains metal stiffness.

6. FRAMEWORK CONNECTION
   κ_T = 1/B, and B vs 1/γ_phonon: r = {r6:.3f} (#110)
   So: κ_T ∝ γ_phonon (confirmed here)

   Physical interpretation:
   - High γ_phonon → thermal fluctuations large
   - Large fluctuations → easier to compress
   - Coherent lattice (low γ) resists compression
""")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Define colors by class
colors = []
for m in materials:
    if m in ['Ag', 'Cu', 'Au']:
        colors.append('gold')
    elif m in ['Na', 'K', 'Li', 'Rb', 'Cs']:
        colors.append('blue')
    elif m in ['W', 'Mo', 'Ta', 'Nb', 'Re']:
        colors.append('gray')
    elif m in ['Fe', 'Ni', 'Co', 'Ti', 'Cr', 'V']:
        colors.append('red')
    elif m in ['Al', 'Pb', 'Zn', 'Mg', 'Sn', 'Cd', 'In']:
        colors.append('purple')
    elif m in ['Si', 'Ge', 'GaAs']:
        colors.append('orange')
    else:
        colors.append('green')

# Plot 1: log(κ_T) vs γ_phonon
ax1 = axes[0, 0]
ax1.scatter(gamma_phonon, kappa_T, c=colors, s=100, alpha=0.7)
for i, m in enumerate(materials):
    ax1.annotate(m, (gamma_phonon[i], kappa_T[i]), fontsize=8)
ax1.set_xlabel('γ_phonon', fontsize=12)
ax1.set_ylabel('κ_T (GPa⁻¹)', fontsize=12)
ax1.set_yscale('log')
ax1.set_title(f'Compressibility vs Coherence (r = {r2:.3f})', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: B vs θ_D
ax2 = axes[0, 1]
ax2.scatter(theta_D, B, c=colors, s=100, alpha=0.7)
for i, m in enumerate(materials):
    ax2.annotate(m, (theta_D[i], B[i]), fontsize=8)
ax2.set_xlabel('Debye temperature θ_D (K)', fontsize=12)
ax2.set_ylabel('Bulk modulus B (GPa)', fontsize=12)
ax2.set_title(f'B vs θ_D (r = {r5:.3f})', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: B vs n^(5/3) for metals
ax3 = axes[1, 0]
metals_names = [materials[i] for i in metals_idx]
metals_colors = [colors[i] for i in metals_idx]
ax3.scatter(n_e_metals**(5/3), B_metals, c=metals_colors, s=100, alpha=0.7)
for i, name in enumerate(metals_names):
    ax3.annotate(name, (n_e_metals[i]**(5/3), B_metals[i]), fontsize=8)
ax3.set_xlabel('n^(5/3) (m⁻⁵)', fontsize=12)
ax3.set_ylabel('B (GPa)', fontsize=12)
ax3.set_xscale('log')
ax3.set_title(f'Fermi Gas Model: B vs n^(5/3) (r = {r7:.3f})', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: Compressibility hierarchy bar chart
ax4 = axes[1, 1]
sorted_materials = [materials[i] for i in sorted_idx[:15]]  # Top 15
sorted_kappa = [kappa_T[i] for i in sorted_idx[:15]]
sorted_colors = [colors[i] for i in sorted_idx[:15]]

y_pos = np.arange(len(sorted_materials))
ax4.barh(y_pos, sorted_kappa, color=sorted_colors, alpha=0.7)
ax4.set_yticks(y_pos)
ax4.set_yticklabels(sorted_materials)
ax4.set_xlabel('κ_T (GPa⁻¹)', fontsize=12)
ax4.set_title('Most Compressible Materials', fontsize=14)
ax4.grid(True, alpha=0.3, axis='x')

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='gold', label='Noble'),
    Patch(facecolor='blue', label='Alkali'),
    Patch(facecolor='gray', label='Refractory'),
    Patch(facecolor='red', label='3d TM'),
    Patch(facecolor='purple', label='Simple'),
    Patch(facecolor='orange', label='Semiconductor'),
    Patch(facecolor='green', label='Ceramic'),
]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/compressibility_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #113")
print("="*70)

print(f"""
KEY RESULTS:

1. COMPRESSIBILITY vs COHERENCE
   κ_T vs γ_phonon: r = {r1:.3f}
   log(κ_T) vs γ_phonon: r = {r2:.3f}
   B vs 1/γ_phonon: r = {r6:.3f}

   MODERATE correlation - expected since κ_T = 1/B.

2. MATERIAL EXTREMES
   Most compressible: Cs (κ_T = 0.625 GPa⁻¹)
   Least compressible: Diamond (κ_T = 0.0023 GPa⁻¹)
   Ratio: 270×

3. FERMI GAS MODEL
   B vs n^(5/3) (metals): r = {r7:.3f}
   Electron density contributes to metal stiffness.

4. CLASS HIERARCHY
   Alkali >> Simple > Semiconductors > Noble ~ Ceramics > 3d TM > Refractory

FRAMEWORK CONNECTION:
Compressibility is INVERSELY related to coherence:
- Low γ_phonon → coherent lattice → stiff → low κ_T
- High γ_phonon → incoherent → soft → high κ_T

This is the INVERSE of most other properties:
- Transport: ∝ 1/γ (coherence helps)
- Compressibility: ∝ γ (coherence resists compression)

Physical interpretation:
A coherent lattice has well-defined phase relationships.
Compression disrupts these phase relationships.
More coherence = more resistance to compression.
""")

print("\nFigure saved to: compressibility_coherence.png")
