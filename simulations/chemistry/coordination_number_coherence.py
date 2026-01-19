#!/usr/bin/env python3
"""
Chemistry Session #123: Coordination Number and Coherence

Test whether coordination number Z relates to coherence parameters.

Z = number of nearest neighbors in crystal structure
FCC: Z = 12, BCC: Z = 8, HCP: Z = 12, Diamond: Z = 4, SC: Z = 6

Hypothesis: Higher Z means more bonds sharing each atom's electrons,
potentially affecting both phonon and electronic coherence.

Physical model:
- More neighbors → distributed bonding → weaker individual bonds
- More neighbors → more vibrational modes → different phonon coherence
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Coordination number data
# Z values by crystal structure
materials = {
    # FCC (Z = 12)
    'Cu':  {'Z': 12, 'structure': 'FCC', 'theta_D': 343, 'E': 130, 'T_m': 1358, 'E_coh': 3.49},
    'Ag':  {'Z': 12, 'structure': 'FCC', 'theta_D': 225, 'E': 83, 'T_m': 1235, 'E_coh': 2.95},
    'Au':  {'Z': 12, 'structure': 'FCC', 'theta_D': 165, 'E': 78, 'T_m': 1337, 'E_coh': 3.81},
    'Al':  {'Z': 12, 'structure': 'FCC', 'theta_D': 428, 'E': 70, 'T_m': 933, 'E_coh': 3.39},
    'Ni':  {'Z': 12, 'structure': 'FCC', 'theta_D': 450, 'E': 200, 'T_m': 1728, 'E_coh': 4.44},
    'Pt':  {'Z': 12, 'structure': 'FCC', 'theta_D': 240, 'E': 168, 'T_m': 2041, 'E_coh': 5.84},
    'Pb':  {'Z': 12, 'structure': 'FCC', 'theta_D': 105, 'E': 16, 'T_m': 601, 'E_coh': 2.03},
    'Ca':  {'Z': 12, 'structure': 'FCC', 'theta_D': 230, 'E': 20, 'T_m': 1115, 'E_coh': 1.84},

    # BCC (Z = 8)
    'Li':  {'Z': 8, 'structure': 'BCC', 'theta_D': 344, 'E': 5, 'T_m': 454, 'E_coh': 1.63},
    'Na':  {'Z': 8, 'structure': 'BCC', 'theta_D': 158, 'E': 10, 'T_m': 371, 'E_coh': 1.11},
    'K':   {'Z': 8, 'structure': 'BCC', 'theta_D': 91, 'E': 3, 'T_m': 337, 'E_coh': 0.93},
    'Fe':  {'Z': 8, 'structure': 'BCC', 'theta_D': 470, 'E': 211, 'T_m': 1811, 'E_coh': 4.28},
    'Cr':  {'Z': 8, 'structure': 'BCC', 'theta_D': 630, 'E': 279, 'T_m': 2180, 'E_coh': 4.10},
    'V':   {'Z': 8, 'structure': 'BCC', 'theta_D': 380, 'E': 128, 'T_m': 2183, 'E_coh': 5.31},
    'Mo':  {'Z': 8, 'structure': 'BCC', 'theta_D': 450, 'E': 329, 'T_m': 2896, 'E_coh': 6.82},
    'W':   {'Z': 8, 'structure': 'BCC', 'theta_D': 400, 'E': 411, 'T_m': 3695, 'E_coh': 8.90},
    'Ta':  {'Z': 8, 'structure': 'BCC', 'theta_D': 240, 'E': 186, 'T_m': 3290, 'E_coh': 8.10},
    'Nb':  {'Z': 8, 'structure': 'BCC', 'theta_D': 275, 'E': 105, 'T_m': 2750, 'E_coh': 7.57},

    # HCP (Z = 12)
    'Mg':  {'Z': 12, 'structure': 'HCP', 'theta_D': 400, 'E': 45, 'T_m': 923, 'E_coh': 1.51},
    'Zn':  {'Z': 12, 'structure': 'HCP', 'theta_D': 327, 'E': 108, 'T_m': 693, 'E_coh': 1.35},
    'Ti':  {'Z': 12, 'structure': 'HCP', 'theta_D': 420, 'E': 116, 'T_m': 1941, 'E_coh': 4.85},
    'Co':  {'Z': 12, 'structure': 'HCP', 'theta_D': 445, 'E': 209, 'T_m': 1768, 'E_coh': 4.39},
    'Be':  {'Z': 12, 'structure': 'HCP', 'theta_D': 1440, 'E': 287, 'T_m': 1560, 'E_coh': 3.32},

    # Diamond (Z = 4)
    'C':   {'Z': 4, 'structure': 'Diamond', 'theta_D': 2230, 'E': 1050, 'T_m': 4000, 'E_coh': 7.37},
    'Si':  {'Z': 4, 'structure': 'Diamond', 'theta_D': 645, 'E': 130, 'T_m': 1687, 'E_coh': 4.63},
    'Ge':  {'Z': 4, 'structure': 'Diamond', 'theta_D': 374, 'E': 103, 'T_m': 1211, 'E_coh': 3.85},
    'Sn':  {'Z': 4, 'structure': 'Diamond', 'theta_D': 200, 'E': 50, 'T_m': 505, 'E_coh': 3.14},  # gray Sn
}

# Calculate coherence parameters
T = 300

for mat in materials:
    data = materials[mat]
    data['gamma_phonon'] = 2 * T / data['theta_D']
    data['E_per_bond'] = data['E_coh'] / (data['Z'] / 2)  # Each bond shared between 2 atoms

# Extract arrays
names = list(materials.keys())
Z = np.array([materials[m]['Z'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
E_modulus = np.array([materials[m]['E'] for m in names])
T_m = np.array([materials[m]['T_m'] for m in names])
E_coh = np.array([materials[m]['E_coh'] for m in names])
E_per_bond = np.array([materials[m]['E_per_bond'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #123: COORDINATION NUMBER AND COHERENCE")
print("=" * 70)

print(f"\nDataset: {len(names)} materials")
print(f"Z range: {np.min(Z)} - {np.max(Z)}")

# Test 1: Z vs γ_phonon
r1, p1 = stats.pearsonr(Z, gamma_phonon)
print(f"\n1. Z vs γ_phonon: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.5:
    print("   Coordination affects phonon coherence")
else:
    print("   Coordination has WEAK effect on phonon coherence")

# Test 2: Z vs θ_D
r2, p2 = stats.pearsonr(Z, theta_D)
print(f"2. Z vs θ_D: r = {r2:.3f}")

# Test 3: Z vs E (Young's modulus)
r3, p3 = stats.pearsonr(Z, E_modulus)
print(f"3. Z vs E: r = {r3:.3f}")

# Test 4: E_coh/Z (bond energy) vs γ_phonon
r4, p4 = stats.pearsonr(E_per_bond, gamma_phonon)
print(f"\n4. E_coh/Z (bond energy) vs γ_phonon: r = {r4:.3f}")
if abs(r4) > 0.7:
    print("   Bond energy determines coherence")

# Test 5: Z × E_coh (total bonding) vs θ_D
Z_E_coh = Z * E_coh
r5, p5 = stats.pearsonr(Z_E_coh, theta_D)
print(f"5. Z × E_coh vs θ_D: r = {r5:.3f}")

# Structure analysis
print("\n" + "=" * 70)
print("STRUCTURE ANALYSIS")
print("=" * 70)

structures = {
    'Diamond (Z=4)': [m for m in names if materials[m]['Z'] == 4],
    'BCC (Z=8)': [m for m in names if materials[m]['Z'] == 8],
    'FCC/HCP (Z=12)': [m for m in names if materials[m]['Z'] == 12],
}

print(f"\n{'Structure':<20} {'Mean θ_D':<10} {'Mean γ_ph':<10} {'Mean E':<10} {'Mean E_coh/Z':<12}")
print("-" * 65)

for struct_name, members in structures.items():
    if len(members) >= 2:
        struct_theta = [materials[m]['theta_D'] for m in members]
        struct_gamma = [materials[m]['gamma_phonon'] for m in members]
        struct_E = [materials[m]['E'] for m in members]
        struct_EZ = [materials[m]['E_per_bond'] for m in members]
        print(f"{struct_name:<20} {np.mean(struct_theta):<10.0f} {np.mean(struct_gamma):<10.2f} "
              f"{np.mean(struct_E):<10.0f} {np.mean(struct_EZ):<12.2f}")

# Compare close-packed (Z=12) vs open (Z=4, Z=8)
print("\n" + "-" * 70)
print("CLOSE-PACKED (Z=12) vs OPEN STRUCTURES (Z=4, Z=8)")
print("-" * 70)

close_packed = [m for m in names if materials[m]['Z'] == 12]
open_struct = [m for m in names if materials[m]['Z'] in [4, 8]]

cp_gamma = [materials[m]['gamma_phonon'] for m in close_packed]
op_gamma = [materials[m]['gamma_phonon'] for m in open_struct]

from scipy.stats import mannwhitneyu
stat, p_mw = mannwhitneyu(cp_gamma, op_gamma, alternative='two-sided')
print(f"\nClose-packed mean γ: {np.mean(cp_gamma):.2f} ± {np.std(cp_gamma):.2f}")
print(f"Open structure mean γ: {np.mean(op_gamma):.2f} ± {np.std(op_gamma):.2f}")
print(f"Mann-Whitney p = {p_mw:.3f}")

if p_mw < 0.05:
    print("Significant difference between close-packed and open structures")
else:
    print("NO significant difference - Z has weak effect on coherence")

# E_coh/Z analysis
print("\n" + "=" * 70)
print("BOND ENERGY ANALYSIS")
print("=" * 70)

print(f"\n{'Element':<5} {'Z':<3} {'E_coh':<8} {'E_coh/Z':<10} {'θ_D':<8} {'γ_phonon':<10}")
print("-" * 55)

# Sort by E_per_bond
sorted_mats = sorted(names, key=lambda m: materials[m]['E_per_bond'], reverse=True)
for m in sorted_mats[:10]:
    print(f"{m:<5} {materials[m]['Z']:<3} {materials[m]['E_coh']:<8.2f} "
          f"{materials[m]['E_per_bond']:<10.2f} {materials[m]['theta_D']:<8} "
          f"{materials[m]['gamma_phonon']:<10.2f}")

# Summary
print("\n" + "=" * 70)
print("SESSION #123 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- Z vs γ_phonon: r = {r1:.3f}
- Z vs θ_D: r = {r2:.3f}
- E_coh/Z vs γ_phonon: r = {r4:.3f}
- Close-packed vs Open: p = {p_mw:.3f}

Physical Interpretation:
""")

if abs(r1) < 0.3:
    print("Coordination number has WEAK effect on phonon coherence")
    print("\nWhy Z doesn't determine γ_phonon:")
    print("  - Diamond (Z=4) has LOWEST γ despite fewest neighbors")
    print("  - This is because bond STRENGTH matters more than bond COUNT")
    print("  - Diamond's sp³ bonds are extremely stiff")
    print("")
    print("The key insight:")
    print("  - γ_phonon is determined by BOND STIFFNESS (k)")
    print("  - k depends on bond length and bond order, not Z")
    print("  - Higher Z → distributed electrons → weaker individual bonds")
    print("  - But the AVERAGE stiffness can still be high/low")
else:
    print("Coordination shows correlation with coherence")

if abs(r4) > 0.5:
    print(f"\nBond energy (E_coh/Z) vs γ_phonon: r = {r4:.3f}")
    print("Bond energy PER BOND is a better coherence predictor than total E_coh")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Z vs γ_phonon
ax1 = axes[0, 0]
colors = {'Diamond': 'purple', 'BCC': 'red', 'FCC': 'blue', 'HCP': 'green'}

for m in names:
    struct = materials[m]['structure']
    ax1.scatter(materials[m]['Z'], materials[m]['gamma_phonon'],
                c=colors[struct], s=80, alpha=0.7, label=struct)
    ax1.annotate(m, (materials[m]['Z'], materials[m]['gamma_phonon']), fontsize=7)

# Remove duplicate legend entries
handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc='best')

ax1.set_xlabel('Coordination Number Z')
ax1.set_ylabel('γ_phonon = 2T/θ_D')
ax1.set_title(f'Z vs γ_phonon\nr = {r1:.3f}')
ax1.grid(True, alpha=0.3)

# Plot 2: E_coh/Z vs γ_phonon
ax2 = axes[0, 1]
for m in names:
    struct = materials[m]['structure']
    ax2.scatter(materials[m]['E_per_bond'], materials[m]['gamma_phonon'],
                c=colors[struct], s=80, alpha=0.7)
    ax2.annotate(m, (materials[m]['E_per_bond'], materials[m]['gamma_phonon']), fontsize=7)

ax2.set_xlabel('E_coh/Z (eV/bond)')
ax2.set_ylabel('γ_phonon = 2T/θ_D')
ax2.set_title(f'Bond Energy vs γ_phonon\nr = {r4:.3f}')
ax2.grid(True, alpha=0.3)

# Plot 3: Z vs θ_D
ax3 = axes[1, 0]
for m in names:
    struct = materials[m]['structure']
    ax3.scatter(materials[m]['Z'], materials[m]['theta_D'],
                c=colors[struct], s=80, alpha=0.7)
    ax3.annotate(m, (materials[m]['Z'], materials[m]['theta_D']), fontsize=6, alpha=0.8)

ax3.set_xlabel('Coordination Number Z')
ax3.set_ylabel('Debye Temperature θ_D (K)')
ax3.set_title(f'Z vs θ_D\nr = {r2:.3f}')
ax3.grid(True, alpha=0.3)

# Plot 4: Structure comparison boxplot
ax4 = axes[1, 1]
data_for_box = []
labels_for_box = []

for struct_name, members in structures.items():
    gammas = [materials[m]['gamma_phonon'] for m in members]
    data_for_box.append(gammas)
    labels_for_box.append(struct_name)

bp = ax4.boxplot(data_for_box, labels=labels_for_box, patch_artist=True)
colors_box = ['purple', 'red', 'blue']
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.5)

ax4.set_ylabel('γ_phonon')
ax4.set_title('γ_phonon by Crystal Structure')
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/coordination_number_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: coordination_number_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r1) < 0.3 and abs(r4) > 0.5:
    print("\n○ COORDINATION is NOT a coherence determinant")
    print(f"  Z vs γ_phonon: r = {r1:.3f} (WEAK)")
    print("\n✓ BOND ENERGY is better coherence proxy")
    print(f"  E_coh/Z vs γ_phonon: r = {r4:.3f}")
    print("\n  Key insight: Bond STRENGTH matters, not bond COUNT")
elif abs(r1) > 0.5:
    print("\n✓ COORDINATION affects coherence")
    print(f"  Z vs γ_phonon: r = {r1:.3f}")
else:
    print("\n○ Mixed results - structure effects are complex")
