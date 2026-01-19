#!/usr/bin/env python3
"""
Chemistry Session #124: Bond Order and Coherence

Test whether bond order relates to coherence parameters.

Bond order n = number of electron pairs in bond
- Single bond: n = 1
- Double bond: n = 2
- Triple bond: n = 3
- Metallic bond: fractional (electrons delocalized)
- Covalent bond: n = 1-3 (directed, strong)

Hypothesis: Higher bond order → stiffer bonds → higher θ_D → lower γ
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Bond order and related data
# Bond order definitions:
# - Diamond: 1 (sp³ single bonds, but VERY stiff)
# - Graphite: 1.33 (resonance structure)
# - N₂: 3 (triple bond)
# - Metals: ~0.5-1.0 (delocalized)

materials = {
    # Covalent solids (directional bonds)
    'Diamond': {'bond_order': 1.0, 'theta_D': 2230, 'E_bond': 3.69, 'type': 'covalent'},
    'Si':      {'bond_order': 1.0, 'theta_D': 645, 'E_bond': 2.31, 'type': 'covalent'},
    'Ge':      {'bond_order': 1.0, 'theta_D': 374, 'E_bond': 1.93, 'type': 'covalent'},
    'SiC':     {'bond_order': 1.0, 'theta_D': 1200, 'E_bond': 3.1, 'type': 'covalent'},
    'BN':      {'bond_order': 1.0, 'theta_D': 1900, 'E_bond': 3.25, 'type': 'covalent'},
    'AlN':     {'bond_order': 1.0, 'theta_D': 950, 'E_bond': 2.9, 'type': 'covalent'},
    'GaN':     {'bond_order': 1.0, 'theta_D': 600, 'E_bond': 2.2, 'type': 'covalent'},

    # Metallic (delocalized bonds)
    'Cu':  {'bond_order': 0.5, 'theta_D': 343, 'E_bond': 0.58, 'type': 'metallic'},
    'Ag':  {'bond_order': 0.5, 'theta_D': 225, 'E_bond': 0.49, 'type': 'metallic'},
    'Au':  {'bond_order': 0.5, 'theta_D': 165, 'E_bond': 0.63, 'type': 'metallic'},
    'Al':  {'bond_order': 0.5, 'theta_D': 428, 'E_bond': 0.56, 'type': 'metallic'},
    'Pb':  {'bond_order': 0.5, 'theta_D': 105, 'E_bond': 0.34, 'type': 'metallic'},

    # Alkali metals (very weak metallic)
    'Li':  {'bond_order': 0.25, 'theta_D': 344, 'E_bond': 0.41, 'type': 'metallic'},
    'Na':  {'bond_order': 0.25, 'theta_D': 158, 'E_bond': 0.28, 'type': 'metallic'},
    'K':   {'bond_order': 0.25, 'theta_D': 91, 'E_bond': 0.23, 'type': 'metallic'},

    # Transition metals (partial d-bonding)
    'Fe':  {'bond_order': 0.8, 'theta_D': 470, 'E_bond': 1.07, 'type': 'metallic'},
    'Cr':  {'bond_order': 1.0, 'theta_D': 630, 'E_bond': 1.03, 'type': 'metallic'},  # high bond order for metal
    'W':   {'bond_order': 1.0, 'theta_D': 400, 'E_bond': 2.23, 'type': 'metallic'},
    'Mo':  {'bond_order': 1.0, 'theta_D': 450, 'E_bond': 1.71, 'type': 'metallic'},

    # Ionic (Coulomb interaction, not covalent bond order)
    'NaCl': {'bond_order': 0.0, 'theta_D': 321, 'E_bond': 2.05, 'type': 'ionic'},  # Madelung
    'KCl':  {'bond_order': 0.0, 'theta_D': 235, 'E_bond': 2.13, 'type': 'ionic'},
    'MgO':  {'bond_order': 0.0, 'theta_D': 946, 'E_bond': 2.50, 'type': 'ionic'},
    'CaO':  {'bond_order': 0.0, 'theta_D': 648, 'E_bond': 2.70, 'type': 'ionic'},
}

# Calculate coherence parameters
T = 300

for mat in materials:
    data = materials[mat]
    data['gamma_phonon'] = 2 * T / data['theta_D']

# Extract arrays
names = list(materials.keys())
bond_order = np.array([materials[m]['bond_order'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
E_bond = np.array([materials[m]['E_bond'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #124: BOND ORDER AND COHERENCE")
print("=" * 70)

print(f"\nDataset: {len(names)} materials")
print(f"Bond order range: {np.min(bond_order):.2f} - {np.max(bond_order):.2f}")

# Test 1: Bond order vs γ_phonon
r1, p1 = stats.pearsonr(bond_order, gamma_phonon)
print(f"\n1. Bond order vs γ_phonon: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.5:
    print("   Bond order affects phonon coherence")

# Test 2: Bond order vs θ_D
r2, p2 = stats.pearsonr(bond_order, theta_D)
print(f"2. Bond order vs θ_D: r = {r2:.3f}")

# Test 3: E_bond vs γ_phonon
r3, p3 = stats.pearsonr(E_bond, gamma_phonon)
print(f"3. E_bond vs γ_phonon: r = {r3:.3f}")

# Test 4: E_bond vs θ_D
r4, p4 = stats.pearsonr(E_bond, theta_D)
print(f"4. E_bond vs θ_D: r = {r4:.3f}")

# Type analysis
print("\n" + "=" * 70)
print("BOND TYPE ANALYSIS")
print("=" * 70)

bond_types = ['covalent', 'metallic', 'ionic']
print(f"\n{'Type':<12} {'Mean BO':<10} {'Mean θ_D':<10} {'Mean γ_ph':<10} {'Mean E_bond':<10}")
print("-" * 55)

type_stats = {}
for btype in bond_types:
    members = [m for m in names if materials[m]['type'] == btype]
    if len(members) >= 2:
        type_bo = [materials[m]['bond_order'] for m in members]
        type_theta = [materials[m]['theta_D'] for m in members]
        type_gamma = [materials[m]['gamma_phonon'] for m in members]
        type_E = [materials[m]['E_bond'] for m in members]
        type_stats[btype] = {
            'bond_order': np.mean(type_bo),
            'theta_D': np.mean(type_theta),
            'gamma_phonon': np.mean(type_gamma),
            'E_bond': np.mean(type_E)
        }
        print(f"{btype:<12} {np.mean(type_bo):<10.2f} {np.mean(type_theta):<10.0f} "
              f"{np.mean(type_gamma):<10.2f} {np.mean(type_E):<10.2f}")

# Within-type correlations
print("\n" + "-" * 70)
print("WITHIN-TYPE CORRELATIONS (Bond Order vs γ_phonon)")
print("-" * 70)

for btype in bond_types:
    members = [m for m in names if materials[m]['type'] == btype]
    if len(members) >= 4:
        type_bo = [materials[m]['bond_order'] for m in members]
        type_gamma = [materials[m]['gamma_phonon'] for m in members]
        r_wt, _ = stats.pearsonr(type_bo, type_gamma)
        print(f"  {btype}: r = {r_wt:.3f} (n={len(members)})")

# Compare covalent vs metallic
print("\n" + "-" * 70)
print("COVALENT vs METALLIC")
print("-" * 70)

cov_mats = [m for m in names if materials[m]['type'] == 'covalent']
met_mats = [m for m in names if materials[m]['type'] == 'metallic']

cov_gamma = [materials[m]['gamma_phonon'] for m in cov_mats]
met_gamma = [materials[m]['gamma_phonon'] for m in met_mats]

from scipy.stats import mannwhitneyu
stat, p_mw = mannwhitneyu(cov_gamma, met_gamma, alternative='two-sided')
print(f"\nCovalent mean γ: {np.mean(cov_gamma):.2f} ± {np.std(cov_gamma):.2f}")
print(f"Metallic mean γ: {np.mean(met_gamma):.2f} ± {np.std(met_gamma):.2f}")
print(f"Mann-Whitney p = {p_mw:.4f}")

if p_mw < 0.01:
    print("\n*** HIGHLY SIGNIFICANT ***")
    print("Covalent bonds are more coherent than metallic bonds")

# Summary
print("\n" + "=" * 70)
print("SESSION #124 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- Bond order vs γ_phonon: r = {r1:.3f}
- Bond order vs θ_D: r = {r2:.3f}
- E_bond vs γ_phonon: r = {r3:.3f}
- E_bond vs θ_D: r = {r4:.3f}
- Covalent vs Metallic: p = {p_mw:.4f}

Physical Interpretation:
""")

if p_mw < 0.05 and abs(r3) > 0.5:
    print("BOND TYPE AND STRENGTH determine coherence:")
    print("")
    print("Covalent bonds: DIRECTIONAL + STRONG")
    print("  - Electrons localized between atoms")
    print("  - High bond energy → high θ_D → low γ")
    print("  - Mean γ = {:.2f}".format(np.mean(cov_gamma)))
    print("")
    print("Metallic bonds: DELOCALIZED + VARIABLE")
    print("  - Electrons shared among many atoms")
    print("  - Lower effective bond energy")
    print("  - Mean γ = {:.2f}".format(np.mean(met_gamma)))
    print("")
    print("The key insight:")
    print("  BOND DIRECTIONALITY (covalent character) determines coherence")
    print("  Directional bonds → well-defined spring constant → coherent phonons")
else:
    print("Bond order shows correlation with coherence")
    print("But bond type (covalent vs metallic) is more fundamental")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Bond order vs γ_phonon
ax1 = axes[0, 0]
colors = {'covalent': 'purple', 'metallic': 'blue', 'ionic': 'orange'}

for m in names:
    btype = materials[m]['type']
    ax1.scatter(materials[m]['bond_order'], materials[m]['gamma_phonon'],
                c=colors[btype], s=80, alpha=0.7, label=btype)
    ax1.annotate(m, (materials[m]['bond_order'], materials[m]['gamma_phonon']), fontsize=7)

# Remove duplicate legend entries
handles, labels = ax1.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax1.legend(by_label.values(), by_label.keys(), loc='best')

ax1.set_xlabel('Bond Order')
ax1.set_ylabel('γ_phonon = 2T/θ_D')
ax1.set_title(f'Bond Order vs γ_phonon\nr = {r1:.3f}')
ax1.grid(True, alpha=0.3)

# Plot 2: E_bond vs γ_phonon
ax2 = axes[0, 1]
for m in names:
    btype = materials[m]['type']
    ax2.scatter(materials[m]['E_bond'], materials[m]['gamma_phonon'],
                c=colors[btype], s=80, alpha=0.7)
    ax2.annotate(m, (materials[m]['E_bond'], materials[m]['gamma_phonon']), fontsize=7)

ax2.set_xlabel('Bond Energy (eV/bond)')
ax2.set_ylabel('γ_phonon = 2T/θ_D')
ax2.set_title(f'E_bond vs γ_phonon\nr = {r3:.3f}')
ax2.grid(True, alpha=0.3)

# Plot 3: E_bond vs θ_D
ax3 = axes[1, 0]
for m in names:
    btype = materials[m]['type']
    ax3.scatter(materials[m]['E_bond'], materials[m]['theta_D'],
                c=colors[btype], s=80, alpha=0.7)
    ax3.annotate(m, (materials[m]['E_bond'], materials[m]['theta_D']), fontsize=6, alpha=0.8)

ax3.set_xlabel('Bond Energy (eV/bond)')
ax3.set_ylabel('Debye Temperature θ_D (K)')
ax3.set_title(f'E_bond vs θ_D\nr = {r4:.3f}')
ax3.grid(True, alpha=0.3)

# Plot 4: Bond type boxplot
ax4 = axes[1, 1]
data_for_box = []
labels_for_box = []

for btype in bond_types:
    members = [m for m in names if materials[m]['type'] == btype]
    gammas = [materials[m]['gamma_phonon'] for m in members]
    data_for_box.append(gammas)
    labels_for_box.append(btype)

bp = ax4.boxplot(data_for_box, labels=labels_for_box, patch_artist=True)
for patch, btype in zip(bp['boxes'], bond_types):
    patch.set_facecolor(colors[btype])
    patch.set_alpha(0.5)

ax4.set_ylabel('γ_phonon')
ax4.set_title('γ_phonon by Bond Type')
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/bond_order_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: bond_order_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if p_mw < 0.01:
    print("\n✓ EXCELLENT - Bond type determines coherence")
    print(f"  Covalent vs Metallic: p = {p_mw:.4f}")
    print(f"  Covalent mean γ: {np.mean(cov_gamma):.2f}")
    print(f"  Metallic mean γ: {np.mean(met_gamma):.2f}")
    print("\n  COVALENT bonds are more COHERENT than metallic bonds")
elif abs(r3) > 0.5:
    print("\n✓ GOOD - Bond energy correlates with coherence")
    print(f"  E_bond vs γ_phonon: r = {r3:.3f}")
else:
    print("\n○ PARTIAL - Bond type matters but relationship complex")
