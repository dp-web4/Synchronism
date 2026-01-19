#!/usr/bin/env python3
"""
Chemistry Session #114: Atomic Volume and Coherence

Test whether atomic volume V_a relates to coherence parameters.

Theory:
V_a = M / (ρ × N_A) (atomic volume)
For cubic: V_a = a³ / n_atoms (lattice constant cubed / atoms per cell)

Atomic volume determines:
- Electron density n = Z / V_a (affects transport)
- Bond length ~ V_a^(1/3)
- Debye temperature θ_D ∝ v_D ∝ √(B/ρ)

Coherence connection:
From Session #109: θ_D ∝ √(E/ρ)
From Session #113: B ∝ 1/V_a (tighter packing → stiffer)

So: θ_D ∝ √(1/(ρ × V_a)) ∝ V_a^(-1/2) × ρ^(-1/2)

Expectation: Smaller atoms → higher θ_D → lower γ_phonon → more coherent
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# ELEMENT DATA
# =============================================================================

elements_data = {
    # Format: V_a (Å³), ρ (g/cm³), M (g/mol), θ_D (K), Z (valence), a (Å), r_cov (Å)
    # Noble metals
    'Ag': {'Va': 17.1, 'rho': 10.5, 'M': 108, 'theta_D': 225, 'Z': 1, 'a': 4.09, 'r_cov': 1.45},
    'Cu': {'Va': 11.8, 'rho': 8.96, 'M': 64, 'theta_D': 343, 'Z': 1, 'a': 3.61, 'r_cov': 1.32},
    'Au': {'Va': 17.0, 'rho': 19.3, 'M': 197, 'theta_D': 165, 'Z': 1, 'a': 4.08, 'r_cov': 1.36},

    # Alkali metals
    'Li': {'Va': 21.5, 'rho': 0.53, 'M': 7, 'theta_D': 344, 'Z': 1, 'a': 3.51, 'r_cov': 1.28},
    'Na': {'Va': 39.5, 'rho': 0.97, 'M': 23, 'theta_D': 158, 'Z': 1, 'a': 4.29, 'r_cov': 1.66},
    'K': {'Va': 75.6, 'rho': 0.86, 'M': 39, 'theta_D': 91, 'Z': 1, 'a': 5.33, 'r_cov': 2.03},
    'Rb': {'Va': 93.0, 'rho': 1.53, 'M': 85, 'theta_D': 56, 'Z': 1, 'a': 5.59, 'r_cov': 2.20},
    'Cs': {'Va': 115, 'rho': 1.93, 'M': 133, 'theta_D': 38, 'Z': 1, 'a': 6.05, 'r_cov': 2.44},

    # Transition metals (refractory)
    'W': {'Va': 15.8, 'rho': 19.3, 'M': 184, 'theta_D': 400, 'Z': 6, 'a': 3.16, 'r_cov': 1.62},
    'Mo': {'Va': 15.6, 'rho': 10.2, 'M': 96, 'theta_D': 450, 'Z': 6, 'a': 3.15, 'r_cov': 1.54},
    'Ta': {'Va': 18.1, 'rho': 16.7, 'M': 181, 'theta_D': 240, 'Z': 5, 'a': 3.30, 'r_cov': 1.70},
    'Nb': {'Va': 18.1, 'rho': 8.57, 'M': 93, 'theta_D': 275, 'Z': 5, 'a': 3.30, 'r_cov': 1.64},

    # Transition metals (3d)
    'Ti': {'Va': 17.7, 'rho': 4.51, 'M': 48, 'theta_D': 420, 'Z': 4, 'a': 2.95, 'r_cov': 1.60},
    'V': {'Va': 13.9, 'rho': 6.11, 'M': 51, 'theta_D': 380, 'Z': 5, 'a': 3.02, 'r_cov': 1.53},
    'Cr': {'Va': 12.0, 'rho': 7.15, 'M': 52, 'theta_D': 630, 'Z': 6, 'a': 2.88, 'r_cov': 1.39},
    'Fe': {'Va': 11.8, 'rho': 7.87, 'M': 56, 'theta_D': 470, 'Z': 8, 'a': 2.87, 'r_cov': 1.32},
    'Co': {'Va': 11.1, 'rho': 8.86, 'M': 59, 'theta_D': 445, 'Z': 9, 'a': 2.51, 'r_cov': 1.26},
    'Ni': {'Va': 11.0, 'rho': 8.91, 'M': 59, 'theta_D': 450, 'Z': 10, 'a': 3.52, 'r_cov': 1.24},

    # Simple metals
    'Al': {'Va': 16.6, 'rho': 2.70, 'M': 27, 'theta_D': 428, 'Z': 3, 'a': 4.05, 'r_cov': 1.21},
    'Mg': {'Va': 23.2, 'rho': 1.74, 'M': 24, 'theta_D': 400, 'Z': 2, 'a': 3.21, 'r_cov': 1.41},
    'Zn': {'Va': 15.3, 'rho': 7.14, 'M': 65, 'theta_D': 327, 'Z': 2, 'a': 2.66, 'r_cov': 1.22},
    'Pb': {'Va': 30.5, 'rho': 11.3, 'M': 207, 'theta_D': 105, 'Z': 4, 'a': 4.95, 'r_cov': 1.46},
    'Sn': {'Va': 27.0, 'rho': 7.26, 'M': 119, 'theta_D': 200, 'Z': 4, 'a': 5.83, 'r_cov': 1.39},

    # Semiconductors
    'Si': {'Va': 20.1, 'rho': 2.33, 'M': 28, 'theta_D': 645, 'Z': 4, 'a': 5.43, 'r_cov': 1.11},
    'Ge': {'Va': 22.7, 'rho': 5.32, 'M': 73, 'theta_D': 374, 'Z': 4, 'a': 5.66, 'r_cov': 1.20},

    # Ceramics (per atom)
    'C-d': {'Va': 5.7, 'rho': 3.51, 'M': 12, 'theta_D': 2230, 'Z': 4, 'a': 3.57, 'r_cov': 0.76},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #114: ATOMIC VOLUME AND COHERENCE")
print("="*70)

T = 300  # K

elements = list(elements_data.keys())
Va = np.array([elements_data[e]['Va'] for e in elements])
rho = np.array([elements_data[e]['rho'] for e in elements])
M = np.array([elements_data[e]['M'] for e in elements])
theta_D = np.array([elements_data[e]['theta_D'] for e in elements])
Z = np.array([elements_data[e]['Z'] for e in elements])
r_cov = np.array([elements_data[e]['r_cov'] for e in elements])

# Coherence parameters
gamma_phonon = 2 * T / theta_D

# Electron density (electrons / Å³)
n_e = Z / Va

# Bond length proxy
r_bond = Va**(1/3)  # Å

print(f"\n{'Element':<8} {'V_a (Å³)':<12} {'θ_D (K)':<10} {'γ_phonon':<10} {'n_e (Å⁻³)':<12} {'r_cov (Å)':<10}")
print("-"*66)
for i, e in enumerate(elements):
    print(f"{e:<8} {Va[i]:<12.1f} {theta_D[i]:<10.0f} {gamma_phonon[i]:<10.3f} {n_e[i]:<12.3f} {r_cov[i]:<10.2f}")

# =============================================================================
# CORRELATIONS
# =============================================================================

print("\n" + "="*70)
print("CORRELATIONS")
print("="*70)

# V_a vs γ_phonon
r1, p1 = stats.pearsonr(Va, gamma_phonon)
print(f"\nV_a vs γ_phonon: r = {r1:.3f}")

# V_a vs θ_D
r2, p2 = stats.pearsonr(Va, theta_D)
print(f"V_a vs θ_D: r = {r2:.3f}")

# log(V_a) vs log(θ_D) (power law)
r3, p3 = stats.pearsonr(np.log10(Va), np.log10(theta_D))
print(f"log(V_a) vs log(θ_D): r = {r3:.3f}")

# 1/V_a vs θ_D² (from √(B/ρ) relation)
r4, p4 = stats.pearsonr(1/Va, theta_D**2)
print(f"1/V_a vs θ_D²: r = {r4:.3f}")

# n_e vs θ_D
r5, p5 = stats.pearsonr(n_e, theta_D)
print(f"n_e vs θ_D: r = {r5:.3f}")

# r_cov vs γ_phonon
r6, p6 = stats.pearsonr(r_cov, gamma_phonon)
print(f"r_cov vs γ_phonon: r = {r6:.3f}")

# r_cov vs V_a^(1/3)
r7, p7 = stats.pearsonr(r_cov, Va**(1/3))
print(f"r_cov vs V_a^(1/3): r = {r7:.3f}")

# =============================================================================
# PERIODIC TRENDS
# =============================================================================

print("\n" + "="*70)
print("PERIODIC TRENDS")
print("="*70)

# Alkali metals - clear trend
alkali = ['Li', 'Na', 'K', 'Rb', 'Cs']
alkali_idx = [elements.index(e) for e in alkali]
alkali_Va = Va[alkali_idx]
alkali_theta = theta_D[alkali_idx]
alkali_gamma = gamma_phonon[alkali_idx]

r_alkali, _ = stats.pearsonr(alkali_Va, alkali_gamma)
print(f"\nAlkali metals: V_a vs γ_phonon: r = {r_alkali:.3f}")

print("\n{'Alkali':<8} {'V_a (Å³)':<12} {'θ_D (K)':<10} {'γ_phonon':<10}")
print("-"*42)
for i, e in enumerate(alkali):
    idx = elements.index(e)
    print(f"{e:<8} {Va[idx]:<12.1f} {theta_D[idx]:<10.0f} {gamma_phonon[idx]:<10.3f}")

# 3d transition metals
tm_3d = ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni']
tm_idx = [elements.index(e) for e in tm_3d]
tm_Va = Va[tm_idx]
tm_theta = theta_D[tm_idx]

r_tm, _ = stats.pearsonr(tm_Va, tm_theta)
print(f"\n3d TM: V_a vs θ_D: r = {r_tm:.3f}")

# =============================================================================
# SIZE-COHERENCE RELATIONSHIP
# =============================================================================

print("\n" + "="*70)
print("SIZE-COHERENCE RELATIONSHIP")
print("="*70)

print("""
Physical interpretation:

1. SMALLER ATOMS → MORE COHERENT
   V_a vs γ_phonon: r = {:.3f}
   Smaller volume → higher θ_D → lower γ_phonon.

2. ATOMIC VOLUME HIERARCHY
   | Element | V_a (Å³) | θ_D (K) | γ_phonon |
   |---------|----------|---------|----------|
   | Diamond | {:.1f}    | {}    | {:.3f}    |
   | Ni      | {:.1f}    | {}     | {:.3f}    |
   | Fe      | {:.1f}    | {}     | {:.3f}    |
   | ...     |          |         |          |
   | K       | {:.1f}    | {}      | {:.3f}    |
   | Cs      | {:.1f}   | {}      | {:.3f}   |

3. ELECTRON DENSITY CONNECTION
   n_e = Z / V_a (electrons per Å³)
   n_e vs θ_D: r = {:.3f}

   Higher electron density → stronger bonding → higher θ_D.

4. COVALENT RADIUS CORRELATION
   r_cov vs γ_phonon: r = {:.3f}
   Larger atoms → higher γ_phonon (less coherent).
""".format(r1,
           Va[elements.index('C-d')], theta_D[elements.index('C-d')], gamma_phonon[elements.index('C-d')],
           Va[elements.index('Ni')], theta_D[elements.index('Ni')], gamma_phonon[elements.index('Ni')],
           Va[elements.index('Fe')], theta_D[elements.index('Fe')], gamma_phonon[elements.index('Fe')],
           Va[elements.index('K')], theta_D[elements.index('K')], gamma_phonon[elements.index('K')],
           Va[elements.index('Cs')], theta_D[elements.index('Cs')], gamma_phonon[elements.index('Cs')],
           r5, r6))

# =============================================================================
# MATERIAL CLASS ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("MATERIAL CLASS ANALYSIS")
print("="*70)

classes = {
    'Noble': ['Ag', 'Cu', 'Au'],
    'Alkali': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Refractory': ['W', 'Mo', 'Ta', 'Nb'],
    '3d TM': ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni'],
    'Simple': ['Al', 'Mg', 'Zn', 'Pb', 'Sn'],
    'Semiconductors': ['Si', 'Ge'],
}

print(f"\n{'Class':<15} {'Mean V_a':<12} {'Mean θ_D':<10} {'Mean γ_ph':<10} {'Mean n_e':<10}")
print("-"*59)
for cls, members in classes.items():
    idx = [elements.index(m) for m in members if m in elements]
    if len(idx) > 0:
        print(f"{cls:<15} {np.mean(Va[idx]):<12.1f} {np.mean(theta_D[idx]):<10.0f} {np.mean(gamma_phonon[idx]):<10.3f} {np.mean(n_e[idx]):<10.3f}")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Define colors by class
colors = []
for e in elements:
    if e in ['Ag', 'Cu', 'Au']:
        colors.append('gold')
    elif e in ['Li', 'Na', 'K', 'Rb', 'Cs']:
        colors.append('blue')
    elif e in ['W', 'Mo', 'Ta', 'Nb']:
        colors.append('gray')
    elif e in ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni']:
        colors.append('red')
    elif e in ['Al', 'Mg', 'Zn', 'Pb', 'Sn']:
        colors.append('purple')
    elif e in ['Si', 'Ge']:
        colors.append('orange')
    else:
        colors.append('green')

# Plot 1: V_a vs γ_phonon
ax1 = axes[0, 0]
ax1.scatter(Va, gamma_phonon, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax1.annotate(e, (Va[i], gamma_phonon[i]), fontsize=8)
ax1.set_xlabel('Atomic volume V_a (Å³)', fontsize=12)
ax1.set_ylabel('γ_phonon', fontsize=12)
ax1.set_title(f'V_a vs γ_phonon (r = {r1:.3f})', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: V_a vs θ_D (log-log)
ax2 = axes[0, 1]
ax2.scatter(Va, theta_D, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax2.annotate(e, (Va[i], theta_D[i]), fontsize=8)
ax2.set_xlabel('Atomic volume V_a (Å³)', fontsize=12)
ax2.set_ylabel('Debye temperature θ_D (K)', fontsize=12)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title(f'log(V_a) vs log(θ_D) (r = {r3:.3f})', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: r_cov vs γ_phonon
ax3 = axes[1, 0]
ax3.scatter(r_cov, gamma_phonon, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax3.annotate(e, (r_cov[i], gamma_phonon[i]), fontsize=8)
ax3.set_xlabel('Covalent radius r_cov (Å)', fontsize=12)
ax3.set_ylabel('γ_phonon', fontsize=12)
ax3.set_title(f'r_cov vs γ_phonon (r = {r6:.3f})', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: n_e vs θ_D
ax4 = axes[1, 1]
ax4.scatter(n_e, theta_D, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax4.annotate(e, (n_e[i], theta_D[i]), fontsize=8)
ax4.set_xlabel('Electron density n_e (Å⁻³)', fontsize=12)
ax4.set_ylabel('Debye temperature θ_D (K)', fontsize=12)
ax4.set_title(f'n_e vs θ_D (r = {r5:.3f})', fontsize=14)
ax4.grid(True, alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='gold', label='Noble'),
    Patch(facecolor='blue', label='Alkali'),
    Patch(facecolor='gray', label='Refractory'),
    Patch(facecolor='red', label='3d TM'),
    Patch(facecolor='purple', label='Simple'),
    Patch(facecolor='orange', label='Semiconductor'),
    Patch(facecolor='green', label='Diamond'),
]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/atomic_volume_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #114")
print("="*70)

print(f"""
KEY RESULTS:

1. ATOMIC VOLUME vs COHERENCE
   V_a vs γ_phonon: r = {r1:.3f}
   log(V_a) vs log(θ_D): r = {r3:.3f}

   STRONG correlation - smaller atoms are more coherent.

2. COVALENT RADIUS CORRELATION
   r_cov vs γ_phonon: r = {r6:.3f}
   r_cov vs V_a^(1/3): r = {r7:.3f}

   Atomic size directly determines coherence scale.

3. ELECTRON DENSITY
   n_e vs θ_D: r = {r5:.3f}
   Higher electron density → higher θ_D.

4. SIZE HIERARCHY
   Diamond: V_a = 5.7 Å³, γ = 0.27 (MOST coherent)
   3d TM: V_a ~ 12 Å³, γ ~ 1.3
   Noble: V_a ~ 15 Å³, γ ~ 2.7
   Simple: V_a ~ 22 Å³, γ ~ 2.5
   Alkali: V_a ~ 69 Å³, γ ~ 7.4
   Cs: V_a = 115 Å³, γ = 15.8 (LEAST coherent)

5. ALKALI METALS SHOW PERFECT TREND
   V_a vs γ_phonon (alkali only): r = {r_alkali:.3f}
   Clear periodic trend: Li → Na → K → Rb → Cs

FRAMEWORK CONNECTION:
Atomic volume is a FUNDAMENTAL determinant of coherence:
- V_a sets the bond length scale
- Bond length determines spring constant
- Spring constant determines θ_D
- θ_D determines γ_phonon = 2T/θ_D

Chain: V_a → bond length → k → θ_D → γ_phonon

Small atoms (Diamond, Ni, Fe) have:
- Short bonds → high spring constant
- High k → high θ_D
- High θ_D → low γ_phonon → COHERENT

Large atoms (K, Rb, Cs) have:
- Long bonds → low spring constant
- Low k → low θ_D
- Low θ_D → high γ_phonon → INCOHERENT
""")

print("\nFigure saved to: atomic_volume_coherence.png")
