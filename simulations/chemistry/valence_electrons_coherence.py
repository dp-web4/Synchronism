#!/usr/bin/env python3
"""
Chemistry Session #125: Valence Electrons and Coherence

Test whether valence electron count relates to coherence parameters.

Valence electrons n_v = electrons in outermost shell available for bonding
- Group 1 (alkali): n_v = 1
- Group 2 (alkaline earth): n_v = 2
- Group 4 (Si, Ge): n_v = 4
- Transition metals: n_v = varies (d-electrons)
- Group 14 (C, Si): n_v = 4

Hypothesis: More valence electrons → more bonding options
This may affect electronic coherence more than phononic coherence.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Valence electron data
# n_v = nominal valence electrons
materials = {
    # Alkali (n_v = 1)
    'Li':  {'n_v': 1, 'theta_D': 344, 'IE': 5.39, 'EN': 0.98, 'sigma': 1.1e7},
    'Na':  {'n_v': 1, 'theta_D': 158, 'IE': 5.14, 'EN': 0.93, 'sigma': 2.1e7},
    'K':   {'n_v': 1, 'theta_D': 91, 'IE': 4.34, 'EN': 0.82, 'sigma': 1.4e7},

    # Alkaline earth (n_v = 2)
    'Be':  {'n_v': 2, 'theta_D': 1440, 'IE': 9.32, 'EN': 1.57, 'sigma': 2.5e7},
    'Mg':  {'n_v': 2, 'theta_D': 400, 'IE': 7.65, 'EN': 1.31, 'sigma': 2.3e7},
    'Ca':  {'n_v': 2, 'theta_D': 230, 'IE': 6.11, 'EN': 1.00, 'sigma': 2.9e7},

    # Group 13 (n_v = 3)
    'Al':  {'n_v': 3, 'theta_D': 428, 'IE': 5.99, 'EN': 1.61, 'sigma': 3.8e7},

    # Group 14 (n_v = 4)
    'C':   {'n_v': 4, 'theta_D': 2230, 'IE': 11.26, 'EN': 2.55, 'sigma': 1e-18},  # insulator
    'Si':  {'n_v': 4, 'theta_D': 645, 'IE': 8.15, 'EN': 1.90, 'sigma': 4e-4},
    'Ge':  {'n_v': 4, 'theta_D': 374, 'IE': 7.90, 'EN': 2.01, 'sigma': 2.2},
    'Sn':  {'n_v': 4, 'theta_D': 200, 'IE': 7.34, 'EN': 1.96, 'sigma': 9.2e6},  # white Sn
    'Pb':  {'n_v': 4, 'theta_D': 105, 'IE': 7.42, 'EN': 2.33, 'sigma': 4.8e6},

    # Transition metals (varied n_v)
    'Ti':  {'n_v': 4, 'theta_D': 420, 'IE': 6.83, 'EN': 1.54, 'sigma': 2.4e6},
    'V':   {'n_v': 5, 'theta_D': 380, 'IE': 6.75, 'EN': 1.63, 'sigma': 5.0e6},
    'Cr':  {'n_v': 6, 'theta_D': 630, 'IE': 6.77, 'EN': 1.66, 'sigma': 7.7e6},
    'Fe':  {'n_v': 8, 'theta_D': 470, 'IE': 7.90, 'EN': 1.83, 'sigma': 1.0e7},
    'Co':  {'n_v': 9, 'theta_D': 445, 'IE': 7.88, 'EN': 1.88, 'sigma': 1.7e7},
    'Ni':  {'n_v': 10, 'theta_D': 450, 'IE': 7.64, 'EN': 1.91, 'sigma': 1.4e7},
    'Cu':  {'n_v': 11, 'theta_D': 343, 'IE': 7.73, 'EN': 1.90, 'sigma': 5.9e7},
    'Zn':  {'n_v': 12, 'theta_D': 327, 'IE': 9.39, 'EN': 1.65, 'sigma': 1.7e7},

    # Noble metals
    'Ag':  {'n_v': 11, 'theta_D': 225, 'IE': 7.58, 'EN': 1.93, 'sigma': 6.3e7},
    'Au':  {'n_v': 11, 'theta_D': 165, 'IE': 9.22, 'EN': 2.54, 'sigma': 4.1e7},

    # 5d transition
    'W':   {'n_v': 6, 'theta_D': 400, 'IE': 7.98, 'EN': 2.36, 'sigma': 1.8e7},
    'Mo':  {'n_v': 6, 'theta_D': 450, 'IE': 7.09, 'EN': 2.16, 'sigma': 2.0e7},
}

# Calculate coherence parameters
T = 300
IE_ref = 13.6

for mat in materials:
    data = materials[mat]
    data['gamma_phonon'] = 2 * T / data['theta_D']
    data['gamma_optical'] = IE_ref / data['IE']

# Extract arrays
names = list(materials.keys())
n_v = np.array([materials[m]['n_v'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
gamma_optical = np.array([materials[m]['gamma_optical'] for m in names])
IE = np.array([materials[m]['IE'] for m in names])
EN = np.array([materials[m]['EN'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #125: VALENCE ELECTRONS AND COHERENCE")
print("=" * 70)

print(f"\nDataset: {len(names)} materials")
print(f"n_v range: {np.min(n_v)} - {np.max(n_v)}")

# Test 1: n_v vs γ_phonon
r1, p1 = stats.pearsonr(n_v, gamma_phonon)
print(f"\n1. n_v vs γ_phonon: r = {r1:.3f}, p = {p1:.2e}")

# Test 2: n_v vs γ_optical
r2, p2 = stats.pearsonr(n_v, gamma_optical)
print(f"2. n_v vs γ_optical: r = {r2:.3f}")

# Test 3: n_v vs θ_D
r3, p3 = stats.pearsonr(n_v, theta_D)
print(f"3. n_v vs θ_D: r = {r3:.3f}")

# Test 4: n_v vs IE
r4, p4 = stats.pearsonr(n_v, IE)
print(f"4. n_v vs IE: r = {r4:.3f}")

# Test 5: n_v vs EN
r5, p5 = stats.pearsonr(n_v, EN)
print(f"5. n_v vs EN: r = {r5:.3f}")

# Group analysis
print("\n" + "=" * 70)
print("GROUP ANALYSIS")
print("=" * 70)

groups = {
    'Alkali (n_v=1)': [m for m in names if materials[m]['n_v'] == 1],
    'Alkaline Earth (n_v=2)': [m for m in names if materials[m]['n_v'] == 2],
    'Group 14 (n_v=4)': [m for m in names if m in ['C', 'Si', 'Ge', 'Sn', 'Pb']],
    'Late TM (n_v>6)': [m for m in names if materials[m]['n_v'] > 6 and m not in ['Sn', 'Pb']],
}

print(f"\n{'Group':<25} {'Mean n_v':<10} {'Mean θ_D':<10} {'Mean γ_ph':<10} {'Mean γ_opt':<10}")
print("-" * 70)

for group_name, members in groups.items():
    if len(members) >= 2:
        g_nv = [materials[m]['n_v'] for m in members]
        g_theta = [materials[m]['theta_D'] for m in members]
        g_gamma_ph = [materials[m]['gamma_phonon'] for m in members]
        g_gamma_opt = [materials[m]['gamma_optical'] for m in members]
        print(f"{group_name:<25} {np.mean(g_nv):<10.1f} {np.mean(g_theta):<10.0f} "
              f"{np.mean(g_gamma_ph):<10.2f} {np.mean(g_gamma_opt):<10.2f}")

# Separate s-block vs d-block
print("\n" + "-" * 70)
print("s-BLOCK vs d-BLOCK COMPARISON")
print("-" * 70)

s_block = ['Li', 'Na', 'K', 'Be', 'Mg', 'Ca']
d_block = ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ag', 'Au', 'W', 'Mo']

s_gamma_ph = [materials[m]['gamma_phonon'] for m in s_block]
d_gamma_ph = [materials[m]['gamma_phonon'] for m in d_block]

s_gamma_opt = [materials[m]['gamma_optical'] for m in s_block]
d_gamma_opt = [materials[m]['gamma_optical'] for m in d_block]

print(f"\ns-block mean γ_phonon: {np.mean(s_gamma_ph):.2f} ± {np.std(s_gamma_ph):.2f}")
print(f"d-block mean γ_phonon: {np.mean(d_gamma_ph):.2f} ± {np.std(d_gamma_ph):.2f}")

print(f"\ns-block mean γ_optical: {np.mean(s_gamma_opt):.2f} ± {np.std(s_gamma_opt):.2f}")
print(f"d-block mean γ_optical: {np.mean(d_gamma_opt):.2f} ± {np.std(d_gamma_opt):.2f}")

from scipy.stats import mannwhitneyu
stat_ph, p_ph = mannwhitneyu(s_gamma_ph, d_gamma_ph)
stat_opt, p_opt = mannwhitneyu(s_gamma_opt, d_gamma_opt)
print(f"\nMann-Whitney (γ_phonon): p = {p_ph:.3f}")
print(f"Mann-Whitney (γ_optical): p = {p_opt:.3f}")

# Summary
print("\n" + "=" * 70)
print("SESSION #125 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- n_v vs γ_phonon: r = {r1:.3f}
- n_v vs γ_optical: r = {r2:.3f}
- n_v vs IE: r = {r4:.3f}
- n_v vs EN: r = {r5:.3f}

Physical Interpretation:
""")

if abs(r1) < 0.3 and abs(r2) < 0.3:
    print("Valence electron count has WEAK effect on coherence")
    print("\nWhy n_v doesn't determine γ:")
    print("  - n_v determines bonding CAPACITY, not bonding STRENGTH")
    print("  - Diamond (n_v=4) and Cu (n_v=11) have very different γ")
    print("  - The difference comes from BOND TYPE, not electron count")
    print("")
    print("Compare:")
    print("  C (n_v=4, covalent): γ_phonon = {:.2f}".format(materials['C']['gamma_phonon']))
    print("  Cu (n_v=11, metallic): γ_phonon = {:.2f}".format(materials['Cu']['gamma_phonon']))
    print("  More electrons ≠ more coherent")
elif abs(r4) > 0.5 or abs(r5) > 0.5:
    print("n_v correlates with IE/EN, which affects electronic coherence")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: n_v vs γ_phonon
ax1 = axes[0, 0]
colors_group = {'Alkali': 'green', 'Alkaline Earth': 'lime', 'Group 14': 'purple',
                'TM': 'blue', 'Noble': 'gold'}

def get_color(mat):
    if mat in ['Li', 'Na', 'K']:
        return 'green'
    elif mat in ['Be', 'Mg', 'Ca']:
        return 'lime'
    elif mat in ['C', 'Si', 'Ge', 'Sn', 'Pb']:
        return 'purple'
    elif mat in ['Cu', 'Ag', 'Au']:
        return 'gold'
    else:
        return 'blue'

for m in names:
    ax1.scatter(materials[m]['n_v'], materials[m]['gamma_phonon'],
                c=get_color(m), s=80, alpha=0.7)
    ax1.annotate(m, (materials[m]['n_v'], materials[m]['gamma_phonon']), fontsize=7)

ax1.set_xlabel('Valence Electrons n_v')
ax1.set_ylabel('γ_phonon = 2T/θ_D')
ax1.set_title(f'n_v vs γ_phonon\nr = {r1:.3f}')
ax1.grid(True, alpha=0.3)

# Plot 2: n_v vs γ_optical
ax2 = axes[0, 1]
for m in names:
    ax2.scatter(materials[m]['n_v'], materials[m]['gamma_optical'],
                c=get_color(m), s=80, alpha=0.7)
    ax2.annotate(m, (materials[m]['n_v'], materials[m]['gamma_optical']), fontsize=7)

ax2.set_xlabel('Valence Electrons n_v')
ax2.set_ylabel('γ_optical = IE_ref/IE')
ax2.set_title(f'n_v vs γ_optical\nr = {r2:.3f}')
ax2.grid(True, alpha=0.3)

# Plot 3: n_v vs IE
ax3 = axes[1, 0]
for m in names:
    ax3.scatter(materials[m]['n_v'], materials[m]['IE'],
                c=get_color(m), s=80, alpha=0.7)
    ax3.annotate(m, (materials[m]['n_v'], materials[m]['IE']), fontsize=6, alpha=0.8)

ax3.set_xlabel('Valence Electrons n_v')
ax3.set_ylabel('Ionization Energy (eV)')
ax3.set_title(f'n_v vs IE\nr = {r4:.3f}')
ax3.grid(True, alpha=0.3)

# Plot 4: s-block vs d-block comparison
ax4 = axes[1, 1]
data_box = [s_gamma_ph, d_gamma_ph]
bp = ax4.boxplot(data_box, labels=['s-block', 'd-block'], patch_artist=True)
bp['boxes'][0].set_facecolor('green')
bp['boxes'][1].set_facecolor('blue')
for box in bp['boxes']:
    box.set_alpha(0.5)

ax4.set_ylabel('γ_phonon')
ax4.set_title(f's-block vs d-block γ_phonon\np = {p_ph:.3f}')
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/valence_electrons_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: valence_electrons_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r1) < 0.3 and abs(r2) < 0.3:
    print("\n○ VALENCE COUNT is NOT a coherence determinant")
    print(f"  n_v vs γ_phonon: r = {r1:.3f} (WEAK)")
    print(f"  n_v vs γ_optical: r = {r2:.3f} (WEAK)")
    print("\n  Valence count determines bonding CAPACITY, not QUALITY")
    print("  Bond TYPE (#124) matters more than electron COUNT")
elif abs(r4) > 0.5:
    print("\n✓ n_v correlates with electronic properties")
    print(f"  n_v vs IE: r = {r4:.3f}")
    print("  But this is periodic table effect, not coherence")
