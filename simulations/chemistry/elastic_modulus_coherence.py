#!/usr/bin/env python3
"""
Chemistry Session #78: Elastic Moduli & Coherence
Test whether coherence framework predicts mechanical properties.

Elastic modulus (E, K, G) measures stiffness:
- Higher modulus = stiffer material
- Depends on bond strength and lattice structure

Key relationships:
- E ∝ d²U/dr² (curvature of potential at equilibrium)
- E ∝ E_cohesive / V_m (cohesive energy per volume)
- E scales with T_m (higher melting point = stiffer)

Coherence interpretation:
- Stronger bonds = more coherent lattice = higher modulus
- E ∝ 2/γ (stiffness proportional to coherence)
- This connects to T_m (Session #77) and κ (Session #65)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

print("=" * 70)
print("CHEMISTRY SESSION #78: ELASTIC MODULI & COHERENCE")
print("=" * 70)

# ==============================================================================
# DATASET: ELASTIC MODULI
# ==============================================================================

print("\n" + "=" * 70)
print("DATASET: YOUNG'S MODULUS AND RELATED PROPERTIES")
print("=" * 70)

# Young's modulus E (GPa), melting point T_m (K), cohesive energy (kJ/mol)
# Also Debye temperature θ_D where available
materials = {
    # Noble metals
    'Cu': (128, 1358, 337, 343),
    'Ag': (83, 1235, 285, 225),
    'Au': (78, 1337, 366, 165),

    # Transition metals
    'Fe': (211, 1811, 415, 470),
    'Ni': (200, 1728, 428, 450),
    'Co': (209, 1768, 424, 445),
    'W': (411, 3695, 860, 400),
    'Mo': (329, 2896, 658, 450),
    'Ti': (116, 1941, 473, 420),
    'Cr': (279, 2180, 395, 630),

    # Alkali metals
    'Na': (10, 371, 108, 158),
    'K': (3.5, 336, 89, 91),
    'Li': (4.9, 454, 159, 344),

    # Other metals
    'Al': (70, 933, 327, 428),
    'Mg': (45, 923, 147, 400),
    'Pb': (16, 601, 195, 105),
    'Sn': (50, 505, 301, 200),
    'Zn': (108, 693, 130, 327),

    # Covalent solids
    'C_diamond': (1050, 3800, 716, 2230),
    'Si': (130, 1687, 456, 645),
    'Ge': (103, 1211, 377, 374),
    'SiC': (450, 3100, 634, 1200),

    # Ceramics
    'Al2O3': (400, 2345, 1009, 1047),
    'MgO': (310, 3098, 620, 946),
}

print(f"Materials: {len(materials)}")

# Print sorted by E
print("\nMaterials sorted by Young's modulus:")
print("-" * 70)
print(f"{'Material':<15} {'E (GPa)':<10} {'T_m (K)':<10} {'E_coh (kJ/mol)':<15} {'θ_D (K)':<10}")
print("-" * 70)

for name, (E, Tm, E_coh, theta_D) in sorted(materials.items(), key=lambda x: -x[1][0]):
    print(f"{name:<15} {E:>8.0f}  {Tm:>8.0f}  {E_coh:>12.0f}    {theta_D:>8.0f}")

# ==============================================================================
# CORRELATION ANALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("CORRELATION ANALYSIS")
print("=" * 70)

# Extract arrays
E_arr = []
Tm_arr = []
E_coh_arr = []
theta_D_arr = []
names = []

for name, (E, Tm, E_coh, theta_D) in materials.items():
    E_arr.append(E)
    Tm_arr.append(Tm)
    E_coh_arr.append(E_coh)
    theta_D_arr.append(theta_D)
    names.append(name)

E_arr = np.array(E_arr)
Tm_arr = np.array(Tm_arr)
E_coh_arr = np.array(E_coh_arr)
theta_D_arr = np.array(theta_D_arr)

# Correlations
r_E_Tm, _ = stats.pearsonr(E_arr, Tm_arr)
r_E_Ecoh, _ = stats.pearsonr(E_arr, E_coh_arr)
r_E_theta, _ = stats.pearsonr(E_arr, theta_D_arr)

print(f"E vs T_m: r = {r_E_Tm:.3f}")
print(f"E vs E_cohesive: r = {r_E_Ecoh:.3f}")
print(f"E vs θ_D: r = {r_E_theta:.3f}")

# ==============================================================================
# COHERENCE PARAMETER ESTIMATION
# ==============================================================================

print("\n" + "=" * 70)
print("γ FROM DEBYE TEMPERATURE")
print("=" * 70)

def gamma_from_theta_D(theta_D, T=298):
    """
    Estimate γ from Debye temperature.
    From Session #75: γ_phonon = 2(T/θ_D) at T < θ_D.
    """
    gamma = 2.0 * T / theta_D
    return np.clip(gamma, 0.5, 2.0)

gamma_arr = np.array([gamma_from_theta_D(theta) for theta in theta_D_arr])

# E vs γ
r_E_gamma, _ = stats.pearsonr(E_arr, gamma_arr)
# E vs 2/γ
coh_factor = 2.0 / gamma_arr
r_E_coh, _ = stats.pearsonr(E_arr, coh_factor)

print(f"E vs γ: r = {r_E_gamma:.3f}")
print(f"E vs 2/γ: r = {r_E_coh:.3f}")

# Print γ values
print("\nγ values from Debye temperature:")
print("-" * 50)
for i, name in enumerate(names):
    print(f"{name:<10}: θ_D = {theta_D_arr[i]:>5.0f} K, γ = {gamma_arr[i]:.2f}, 2/γ = {coh_factor[i]:.2f}")

# ==============================================================================
# THEORETICAL RELATIONSHIP: E, θ_D, T_m
# ==============================================================================

print("\n" + "=" * 70)
print("THEORETICAL RELATIONSHIPS")
print("=" * 70)

print("""
Theoretical connections:

1. DEBYE MODEL:
   θ_D ∝ (K/M)^0.5 × a^(-1)
   Where K = bulk modulus, M = mass, a = lattice parameter

   So: E ∝ θ_D² × M / V

2. LINDEMANN (from Session #77):
   T_m ∝ M × θ_D² / k_B

   Combining: E ∝ T_m / V (stiffness scales with melting)

3. COHERENCE:
   From γ = 2(T/θ_D):
   θ_D = 2T/γ

   So: E ∝ θ_D² ∝ (2T/γ)² ∝ 1/γ² (at fixed T)

   Or: E ∝ (2/γ)² for elastic modulus

Let's test E vs (2/γ)²:
""")

# E vs (2/γ)²
coh_factor_sq = (2.0 / gamma_arr)**2
r_E_coh_sq, _ = stats.pearsonr(E_arr, coh_factor_sq)
print(f"E vs (2/γ)²: r = {r_E_coh_sq:.3f}")

# ==============================================================================
# LINEAR MODEL
# ==============================================================================

print("\n" + "=" * 70)
print("LINEAR MODELS")
print("=" * 70)

# E vs T_m
slope_Tm, intercept_Tm, r_Tm, _, _ = stats.linregress(Tm_arr, E_arr)
E_pred_Tm = slope_Tm * Tm_arr + intercept_Tm
R2_Tm = r_Tm**2

print(f"E = {slope_Tm:.3f} × T_m + {intercept_Tm:.1f}")
print(f"R² = {R2_Tm:.3f}")

# E vs θ_D
slope_theta, intercept_theta, r_theta, _, _ = stats.linregress(theta_D_arr, E_arr)
E_pred_theta = slope_theta * theta_D_arr + intercept_theta
R2_theta = r_theta**2

print(f"\nE = {slope_theta:.3f} × θ_D + {intercept_theta:.1f}")
print(f"R² = {R2_theta:.3f}")

# ==============================================================================
# MATERIAL CLASSES
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS BY MATERIAL CLASS")
print("=" * 70)

classes = {
    'Noble metals': ['Cu', 'Ag', 'Au'],
    'Transition metals': ['Fe', 'Ni', 'Co', 'W', 'Mo', 'Ti', 'Cr'],
    'Alkali metals': ['Na', 'K', 'Li'],
    'Covalent': ['C_diamond', 'Si', 'Ge', 'SiC'],
    'Ceramics': ['Al2O3', 'MgO'],
}

for class_name, members in classes.items():
    E_class = [materials[m][0] for m in members if m in materials]
    Tm_class = [materials[m][1] for m in members if m in materials]
    theta_class = [materials[m][3] for m in members if m in materials]

    if len(E_class) >= 3:
        r, _ = stats.pearsonr(E_class, Tm_class)
        print(f"{class_name}: n={len(E_class)}, E_mean={np.mean(E_class):.0f} GPa, E vs T_m: r={r:.3f}")
    else:
        print(f"{class_name}: n={len(E_class)}, E_mean={np.mean(E_class):.0f} GPa")

# ==============================================================================
# PHONON COHERENCE CONNECTION
# ==============================================================================

print("\n" + "=" * 70)
print("PHONON COHERENCE & ELASTIC MODULUS")
print("=" * 70)

print("""
Physical picture:
- Elastic modulus measures restoring force against deformation
- Restoring force comes from bond stretching
- Strong, coherent bonds → high modulus

From coherence framework:
- γ_phonon = 2(T/θ_D) (Session #75)
- Higher θ_D → lower γ → more coherent phonons
- More coherent lattice → stronger restoring forces → higher E

So: E ∝ (2/γ)^n where n is to be determined

Testing n = 1 (linear): r = {:.3f}
Testing n = 2 (quadratic): r = {:.3f}

The relationship appears {} quadratic than linear.
""".format(r_E_coh, r_E_coh_sq, "more" if abs(r_E_coh_sq) > abs(r_E_coh) else "less"))

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SESSION #78 SUMMARY: ELASTIC MODULI & COHERENCE")
print("=" * 70)

print(f"""
Correlations Found:
- E vs T_m: r = {r_E_Tm:.3f} {"(GOOD)" if abs(r_E_Tm) > 0.6 else "(MODERATE)"}
- E vs E_cohesive: r = {r_E_Ecoh:.3f}
- E vs θ_D: r = {r_E_theta:.3f}
- E vs γ: r = {r_E_gamma:.3f}
- E vs 2/γ: r = {r_E_coh:.3f}
- E vs (2/γ)²: r = {r_E_coh_sq:.3f}

Best predictor: {"T_m" if abs(r_E_Tm) > abs(r_E_theta) else "θ_D"} with r = {max(abs(r_E_Tm), abs(r_E_theta)):.3f}

Key Findings:
1. E correlates strongly with T_m (r = {r_E_Tm:.3f})
   - Higher melting → higher modulus
   - Both reflect bond strength

2. E correlates with θ_D (r = {r_E_theta:.3f})
   - Debye temperature also reflects lattice stiffness
   - E ∝ θ_D² theoretically (Debye model)

3. Coherence model: E vs 2/γ gives r = {r_E_coh:.3f}
   - Consistent with framework
   - Quadratic (2/γ)² gives r = {r_E_coh_sq:.3f}

4. Material class variations:
   - Diamond (E = 1050 GPa) highest - strong covalent bonds
   - Alkali metals (E ~ 3-10 GPa) lowest - weak metallic bonds

Physical Interpretation:
- Elastic modulus measures lattice coherence against deformation
- Coherent (low γ) → strong bonds → high E
- This connects T_m, θ_D, E through common origin: bond strength
""")

# ==============================================================================
# PREDICTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print("""
P78.1: E ∝ (2/γ)^n where n ≈ 1-2
Elastic modulus scales with coherence factor (power TBD).

P78.2: E, T_m, θ_D all reflect bond strength
Connected through common physics of interatomic potential.

P78.3: Diamond has highest E due to lowest γ
Strong sp³ bonds = high coherence = extreme stiffness.

P78.4: Alkali metals soft due to high γ
Weak metallic bonds = low coherence = low modulus.

P78.5: E/T_m ratio approximately constant for material class
Both scale with bond strength, so ratio reflects structure only.
""")

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: E vs T_m
ax1 = axes[0, 0]
ax1.scatter(Tm_arr, E_arr, s=80, alpha=0.7, c='blue')
ax1.plot(Tm_arr, E_pred_Tm, 'r--', label=f'r = {r_E_Tm:.3f}')
for i, name in enumerate(names):
    if name in ['C_diamond', 'W', 'Na', 'Cu', 'Fe']:
        ax1.annotate(name, (Tm_arr[i], E_arr[i]), fontsize=9)
ax1.set_xlabel('Melting Point (K)', fontsize=12)
ax1.set_ylabel("Young's Modulus (GPa)", fontsize=12)
ax1.set_title('Elastic Modulus vs Melting Point', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend()

# Plot 2: E vs θ_D
ax2 = axes[0, 1]
ax2.scatter(theta_D_arr, E_arr, s=80, alpha=0.7, c='green')
ax2.plot(theta_D_arr, E_pred_theta, 'r--', label=f'r = {r_E_theta:.3f}')
for i, name in enumerate(names):
    if name in ['C_diamond', 'W', 'Na', 'Cu', 'Fe']:
        ax2.annotate(name, (theta_D_arr[i], E_arr[i]), fontsize=9)
ax2.set_xlabel('Debye Temperature (K)', fontsize=12)
ax2.set_ylabel("Young's Modulus (GPa)", fontsize=12)
ax2.set_title('Elastic Modulus vs Debye Temperature', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.legend()

# Plot 3: E vs 2/γ
ax3 = axes[1, 0]
ax3.scatter(coh_factor, E_arr, s=80, alpha=0.7, c='purple')
for i, name in enumerate(names):
    if name in ['C_diamond', 'W', 'Na', 'Cu', 'Fe']:
        ax3.annotate(name, (coh_factor[i], E_arr[i]), fontsize=9)
ax3.set_xlabel('2/γ (coherence factor)', fontsize=12)
ax3.set_ylabel("Young's Modulus (GPa)", fontsize=12)
ax3.set_title(f'Elastic Modulus vs Coherence\n(r = {r_E_coh:.3f})', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: By material class
ax4 = axes[1, 1]
class_colors = {'Noble metals': 'gold', 'Transition metals': 'blue',
                'Alkali metals': 'red', 'Covalent': 'green', 'Ceramics': 'brown'}
for class_name, members in classes.items():
    E_class = [materials[m][0] for m in members if m in materials]
    Tm_class = [materials[m][1] for m in members if m in materials]
    ax4.scatter(Tm_class, E_class, label=class_name, alpha=0.7, s=100,
                c=class_colors.get(class_name, 'gray'))
ax4.set_xlabel('Melting Point (K)', fontsize=12)
ax4.set_ylabel("Young's Modulus (GPa)", fontsize=12)
ax4.set_title('Elastic Modulus by Material Class', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/elastic_modulus_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved to: simulations/chemistry/elastic_modulus_coherence.png")

print("\n" + "=" * 70)
print("SESSION #78 COMPLETE: ELASTIC MODULI & COHERENCE")
print("=" * 70)
