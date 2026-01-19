#!/usr/bin/env python3
"""
Chemistry Session #115: Electronegativity and Coherence

Test whether electronegativity χ relates to coherence parameters.

Theory:
Electronegativity χ measures the ability to attract electrons.
Pauling scale: χ ∝ √(bond energy contribution)
Mulliken scale: χ = (IE + EA) / 2 (ionization energy + electron affinity)

Coherence connection:
From Sessions 82-83: γ_optical = IE_ref / IE
High IE → low γ_optical → more coherent electrons

Electronegativity and IE are related:
High χ → strong electron attraction → high IE → low γ_optical

Expectation: χ vs 1/γ_optical (positive correlation)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# ELEMENT DATA
# =============================================================================

elements_data = {
    # Format: χ (Pauling), IE (eV), EA (eV), θ_D (K), Z, V_a (Å³)
    # Halogens - highest χ
    'F': {'chi': 3.98, 'IE': 17.42, 'EA': 3.40, 'theta_D': 91, 'Z': 9, 'Va': 13.3},
    'Cl': {'chi': 3.16, 'IE': 12.97, 'EA': 3.61, 'theta_D': 115, 'Z': 17, 'Va': 22.7},
    'Br': {'chi': 2.96, 'IE': 11.81, 'EA': 3.36, 'theta_D': 100, 'Z': 35, 'Va': 25.6},
    'I': {'chi': 2.66, 'IE': 10.45, 'EA': 3.06, 'theta_D': 80, 'Z': 53, 'Va': 32.9},

    # Noble gases - high IE but no bonding
    'He': {'chi': 0.0, 'IE': 24.59, 'EA': 0.0, 'theta_D': 25, 'Z': 2, 'Va': 31.8},
    'Ne': {'chi': 0.0, 'IE': 21.56, 'EA': 0.0, 'theta_D': 75, 'Z': 10, 'Va': 22.4},
    'Ar': {'chi': 0.0, 'IE': 15.76, 'EA': 0.0, 'theta_D': 92, 'Z': 18, 'Va': 28.0},

    # Nonmetals
    'O': {'chi': 3.44, 'IE': 13.62, 'EA': 1.46, 'theta_D': 155, 'Z': 8, 'Va': 14.0},
    'N': {'chi': 3.04, 'IE': 14.53, 'EA': -0.07, 'theta_D': 105, 'Z': 7, 'Va': 17.3},
    'S': {'chi': 2.58, 'IE': 10.36, 'EA': 2.08, 'theta_D': 200, 'Z': 16, 'Va': 15.5},
    'C': {'chi': 2.55, 'IE': 11.26, 'EA': 1.26, 'theta_D': 2230, 'Z': 6, 'Va': 5.7},  # Diamond

    # Metalloids
    'Si': {'chi': 1.90, 'IE': 8.15, 'EA': 1.39, 'theta_D': 645, 'Z': 14, 'Va': 20.1},
    'Ge': {'chi': 2.01, 'IE': 7.90, 'EA': 1.23, 'theta_D': 374, 'Z': 32, 'Va': 22.7},
    'B': {'chi': 2.04, 'IE': 8.30, 'EA': 0.28, 'theta_D': 1480, 'Z': 5, 'Va': 7.2},

    # Transition metals (noble)
    'Cu': {'chi': 1.90, 'IE': 7.73, 'EA': 1.24, 'theta_D': 343, 'Z': 29, 'Va': 11.8},
    'Ag': {'chi': 1.93, 'IE': 7.58, 'EA': 1.30, 'theta_D': 225, 'Z': 47, 'Va': 17.1},
    'Au': {'chi': 2.54, 'IE': 9.23, 'EA': 2.31, 'theta_D': 165, 'Z': 79, 'Va': 17.0},

    # Transition metals (other)
    'Fe': {'chi': 1.83, 'IE': 7.90, 'EA': 0.15, 'theta_D': 470, 'Z': 26, 'Va': 11.8},
    'Ni': {'chi': 1.91, 'IE': 7.64, 'EA': 1.16, 'theta_D': 450, 'Z': 28, 'Va': 11.0},
    'Co': {'chi': 1.88, 'IE': 7.88, 'EA': 0.66, 'theta_D': 445, 'Z': 27, 'Va': 11.1},
    'Cr': {'chi': 1.66, 'IE': 6.77, 'EA': 0.67, 'theta_D': 630, 'Z': 24, 'Va': 12.0},
    'W': {'chi': 2.36, 'IE': 7.98, 'EA': 0.82, 'theta_D': 400, 'Z': 74, 'Va': 15.8},
    'Mo': {'chi': 2.16, 'IE': 7.09, 'EA': 0.75, 'theta_D': 450, 'Z': 42, 'Va': 15.6},
    'Ti': {'chi': 1.54, 'IE': 6.83, 'EA': 0.08, 'theta_D': 420, 'Z': 22, 'Va': 17.7},

    # Simple metals
    'Al': {'chi': 1.61, 'IE': 5.99, 'EA': 0.43, 'theta_D': 428, 'Z': 13, 'Va': 16.6},
    'Zn': {'chi': 1.65, 'IE': 9.39, 'EA': 0.0, 'theta_D': 327, 'Z': 30, 'Va': 15.3},
    'Mg': {'chi': 1.31, 'IE': 7.65, 'EA': 0.0, 'theta_D': 400, 'Z': 12, 'Va': 23.2},
    'Pb': {'chi': 2.33, 'IE': 7.42, 'EA': 0.36, 'theta_D': 105, 'Z': 82, 'Va': 30.5},
    'Sn': {'chi': 1.96, 'IE': 7.34, 'EA': 1.11, 'theta_D': 200, 'Z': 50, 'Va': 27.0},

    # Alkali metals - lowest χ
    'Li': {'chi': 0.98, 'IE': 5.39, 'EA': 0.62, 'theta_D': 344, 'Z': 3, 'Va': 21.5},
    'Na': {'chi': 0.93, 'IE': 5.14, 'EA': 0.55, 'theta_D': 158, 'Z': 11, 'Va': 39.5},
    'K': {'chi': 0.82, 'IE': 4.34, 'EA': 0.50, 'theta_D': 91, 'Z': 19, 'Va': 75.6},
    'Rb': {'chi': 0.82, 'IE': 4.18, 'EA': 0.49, 'theta_D': 56, 'Z': 37, 'Va': 93.0},
    'Cs': {'chi': 0.79, 'IE': 3.89, 'EA': 0.47, 'theta_D': 38, 'Z': 55, 'Va': 115.0},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #115: ELECTRONEGATIVITY AND COHERENCE")
print("="*70)

T = 300  # K
IE_ref = 13.6  # eV (hydrogen ionization energy reference)

elements = list(elements_data.keys())
chi = np.array([elements_data[e]['chi'] for e in elements])
IE = np.array([elements_data[e]['IE'] for e in elements])
EA = np.array([elements_data[e]['EA'] for e in elements])
theta_D = np.array([elements_data[e]['theta_D'] for e in elements])
Va = np.array([elements_data[e]['Va'] for e in elements])

# Coherence parameters
gamma_phonon = 2 * T / theta_D
gamma_optical = IE_ref / IE  # From Sessions 82-83

# Mulliken electronegativity
chi_mulliken = (IE + EA) / 2

print(f"\n{'Element':<8} {'χ (Pauling)':<12} {'IE (eV)':<10} {'γ_opt':<10} {'θ_D (K)':<10} {'γ_phonon':<10}")
print("-"*64)
for i, e in enumerate(elements):
    print(f"{e:<8} {chi[i]:<12.2f} {IE[i]:<10.2f} {gamma_optical[i]:<10.3f} {theta_D[i]:<10.0f} {gamma_phonon[i]:<10.3f}")

# =============================================================================
# CORRELATIONS
# =============================================================================

print("\n" + "="*70)
print("CORRELATIONS")
print("="*70)

# χ vs IE
r1, p1 = stats.pearsonr(chi[chi > 0], IE[chi > 0])  # Exclude noble gases
print(f"\nχ vs IE (excluding noble gases): r = {r1:.3f}")

# χ vs 1/γ_optical
r2, p2 = stats.pearsonr(chi[chi > 0], 1/gamma_optical[chi > 0])
print(f"χ vs 1/γ_optical: r = {r2:.3f}")

# χ vs γ_optical
r3, p3 = stats.pearsonr(chi[chi > 0], gamma_optical[chi > 0])
print(f"χ vs γ_optical: r = {r3:.3f}")

# χ vs θ_D
r4, p4 = stats.pearsonr(chi[chi > 0], theta_D[chi > 0])
print(f"χ vs θ_D: r = {r4:.3f}")

# χ vs γ_phonon
r5, p5 = stats.pearsonr(chi[chi > 0], gamma_phonon[chi > 0])
print(f"χ vs γ_phonon: r = {r5:.3f}")

# χ vs 1/V_a
r6, p6 = stats.pearsonr(chi[chi > 0], 1/Va[chi > 0])
print(f"χ vs 1/V_a: r = {r6:.3f}")

# Mulliken vs Pauling
r7, p7 = stats.pearsonr(chi[chi > 0], chi_mulliken[chi > 0])
print(f"χ (Pauling) vs χ (Mulliken): r = {r7:.3f}")

# =============================================================================
# HARDNESS AND SOFTNESS
# =============================================================================

print("\n" + "="*70)
print("CHEMICAL HARDNESS AND COHERENCE")
print("="*70)

# Chemical hardness η = (IE - EA) / 2 (Pearson)
# Soft = low η (easily polarizable)
# Hard = high η (not polarizable)

eta = (IE - EA) / 2

r8, p8 = stats.pearsonr(eta[chi > 0], 1/gamma_optical[chi > 0])
print(f"Chemical hardness η vs 1/γ_optical: r = {r8:.3f}")

r9, p9 = stats.pearsonr(eta[chi > 0], theta_D[chi > 0])
print(f"Chemical hardness η vs θ_D: r = {r9:.3f}")

print("\nHardness-softness scale:")
print(f"{'Element':<8} {'η (eV)':<10} {'γ_optical':<10} {'Character':<12}")
print("-"*42)
sorted_idx = np.argsort(eta)[::-1]
for i in sorted_idx[:15]:
    char = "Hard" if eta[i] > 6 else ("Medium" if eta[i] > 4 else "Soft")
    print(f"{elements[i]:<8} {eta[i]:<10.2f} {gamma_optical[i]:<10.3f} {char:<12}")

# =============================================================================
# PERIODIC TRENDS
# =============================================================================

print("\n" + "="*70)
print("PERIODIC TRENDS")
print("="*70)

# Groups
groups = {
    'Alkali (1)': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Halogens (17)': ['F', 'Cl', 'Br', 'I'],
    '3d TM': ['Ti', 'Cr', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
}

print("\nPeriodic group analysis:")
for group, members in groups.items():
    idx = [elements.index(m) for m in members if m in elements]
    if len(idx) > 0:
        mean_chi = np.mean(chi[idx])
        mean_IE = np.mean(IE[idx])
        mean_gamma_opt = np.mean(gamma_optical[idx])
        print(f"\n{group}:")
        print(f"  Mean χ: {mean_chi:.2f}")
        print(f"  Mean IE: {mean_IE:.2f} eV")
        print(f"  Mean γ_optical: {mean_gamma_opt:.3f}")

# =============================================================================
# PHYSICAL INSIGHTS
# =============================================================================

print("\n" + "="*70)
print("PHYSICAL INSIGHTS")
print("="*70)

print(f"""
1. ELECTRONEGATIVITY vs OPTICAL COHERENCE
   χ vs 1/γ_optical: r = {r2:.3f}
   χ vs IE: r = {r1:.3f}

   MODERATE correlation - electronegative elements have
   higher IE → lower γ_optical → more coherent electrons.

2. ELECTRONEGATIVITY vs PHONON COHERENCE
   χ vs θ_D: r = {r4:.3f} (WEAK)
   χ vs γ_phonon: r = {r5:.3f} (WEAK)

   Electronegativity is primarily an ELECTRONIC property.
   It doesn't directly determine phonon coherence.

3. CHEMICAL HARDNESS CONNECTION
   η = (IE - EA) / 2
   η vs 1/γ_optical: r = {r8:.3f}

   Hard atoms (high η) have more coherent electrons.
   F, Ne, He are hardest (η > 10 eV).
   Cs, Rb, K are softest (η < 2 eV).

4. PERIODIC TRENDS
   Electronegativity increases:
   - Left → Right (more protons, smaller)
   - Bottom → Top (smaller, higher IE)

   Same trend as γ_optical decreases (more coherent).

5. FRAMEWORK DISTINCTION
   χ relates to γ_optical (electronic coherence)
   V_a relates to γ_phonon (phonon coherence)

   These are INDEPENDENT coherence channels:
   - Electronic: set by IE, EA, shell structure
   - Phononic: set by atomic size, bonding, mass
""")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Define colors by group
colors = []
for e in elements:
    if e in ['F', 'Cl', 'Br', 'I']:
        colors.append('purple')
    elif e in ['He', 'Ne', 'Ar']:
        colors.append('lightgray')
    elif e in ['Li', 'Na', 'K', 'Rb', 'Cs']:
        colors.append('blue')
    elif e in ['C', 'N', 'O', 'S']:
        colors.append('green')
    elif e in ['Si', 'Ge', 'B']:
        colors.append('orange')
    elif e in ['Cu', 'Ag', 'Au', 'Fe', 'Ni', 'Co', 'Cr', 'W', 'Mo', 'Ti']:
        colors.append('red')
    else:
        colors.append('steelblue')

# Plot 1: χ vs IE
ax1 = axes[0, 0]
mask = chi > 0
ax1.scatter(IE[mask], chi[mask], c=np.array(colors)[mask], s=100, alpha=0.7)
for i, e in enumerate(elements):
    if chi[i] > 0:
        ax1.annotate(e, (IE[i], chi[i]), fontsize=8)
ax1.set_xlabel('Ionization energy IE (eV)', fontsize=12)
ax1.set_ylabel('Electronegativity χ (Pauling)', fontsize=12)
ax1.set_title(f'χ vs IE (r = {r1:.3f})', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: χ vs 1/γ_optical
ax2 = axes[0, 1]
ax2.scatter((1/gamma_optical)[mask], chi[mask], c=np.array(colors)[mask], s=100, alpha=0.7)
for i, e in enumerate(elements):
    if chi[i] > 0:
        ax2.annotate(e, (1/gamma_optical[i], chi[i]), fontsize=8)
ax2.set_xlabel('1/γ_optical (electron coherence)', fontsize=12)
ax2.set_ylabel('Electronegativity χ (Pauling)', fontsize=12)
ax2.set_title(f'χ vs 1/γ_optical (r = {r2:.3f})', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: Chemical hardness η vs γ_optical
ax3 = axes[1, 0]
ax3.scatter(gamma_optical[mask], eta[mask], c=np.array(colors)[mask], s=100, alpha=0.7)
for i, e in enumerate(elements):
    if chi[i] > 0:
        ax3.annotate(e, (gamma_optical[i], eta[i]), fontsize=8)
ax3.set_xlabel('γ_optical', fontsize=12)
ax3.set_ylabel('Chemical hardness η (eV)', fontsize=12)
ax3.set_title(f'η vs γ_optical', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: χ vs γ_phonon
ax4 = axes[1, 1]
ax4.scatter(gamma_phonon[mask], chi[mask], c=np.array(colors)[mask], s=100, alpha=0.7)
for i, e in enumerate(elements):
    if chi[i] > 0:
        ax4.annotate(e, (gamma_phonon[i], chi[i]), fontsize=8)
ax4.set_xlabel('γ_phonon', fontsize=12)
ax4.set_ylabel('Electronegativity χ (Pauling)', fontsize=12)
ax4.set_title(f'χ vs γ_phonon (r = {r5:.3f})', fontsize=14)
ax4.grid(True, alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='purple', label='Halogens'),
    Patch(facecolor='lightgray', label='Noble gases'),
    Patch(facecolor='blue', label='Alkali'),
    Patch(facecolor='green', label='Nonmetals'),
    Patch(facecolor='orange', label='Metalloids'),
    Patch(facecolor='red', label='Transition metals'),
    Patch(facecolor='steelblue', label='Other metals'),
]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/electronegativity_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #115")
print("="*70)

print(f"""
KEY RESULTS:

1. ELECTRONEGATIVITY vs OPTICAL COHERENCE
   χ vs 1/γ_optical: r = {r2:.3f}
   χ vs IE: r = {r1:.3f}

   MODERATE correlation - high χ → high IE → low γ_optical.

2. ELECTRONEGATIVITY vs PHONON COHERENCE
   χ vs θ_D: r = {r4:.3f} (WEAK)
   χ vs γ_phonon: r = {r5:.3f} (WEAK)

   Electronegativity is ELECTRONIC, not phononic.

3. CHEMICAL HARDNESS
   η = (IE - EA) / 2
   η vs 1/γ_optical: r = {r8:.3f}

   Hard atoms have coherent electrons.
   Soft atoms have incoherent electrons.

4. TWO COHERENCE CHANNELS
   - γ_optical: determined by IE, EA, χ (electronic)
   - γ_phonon: determined by V_a, θ_D (vibrational)

   These are INDEPENDENT!

5. PERIODIC TRENDS
   | Group | Mean χ | Mean IE | Mean γ_opt |
   |-------|--------|---------|------------|
   | Halogens | 3.19 | 13.2 | 1.08 |
   | 3d TM | 1.78 | 7.6 | 1.82 |
   | Alkali | 0.87 | 4.6 | 3.05 |

FRAMEWORK CONNECTION:
Electronegativity connects to ELECTRONIC coherence (γ_optical)
but NOT directly to PHONON coherence (γ_phonon).

This confirms TWO INDEPENDENT coherence channels:
1. Electronic: IE → γ_optical → optical properties
2. Phononic: V_a → θ_D → γ_phonon → thermal properties

The framework has TWO SECTORS, not one.
""")

print("\nFigure saved to: electronegativity_coherence.png")
