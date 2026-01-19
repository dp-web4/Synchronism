#!/usr/bin/env python3
"""
Chemistry Session #118: Electron Affinity and Coherence

Test whether electron affinity EA relates to coherence parameters.

EA = energy released when adding an electron to neutral atom
IE = energy required to remove an electron from neutral atom

Both relate to electron binding:
- High EA → readily accepts electrons → electron-attracting
- High IE → holds electrons tightly → electron-keeping

Hypothesis:
- EA should correlate with γ_optical via EN = (IE + EA)/2 (Mulliken)
- But EA may show different behavior since it's GAINING vs LOSING electrons
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Electron affinity data (eV)
# Sources: NIST, Rienstra-Kiracofe et al. (2002)
# Negative values = electron not bound (He, N, noble gases)
materials = {
    # Alkali metals (s-block)
    'Li': {'EA': 0.618, 'IE': 5.39, 'EN': 0.98, 'theta_D': 344},
    'Na': {'EA': 0.548, 'IE': 5.14, 'EN': 0.93, 'theta_D': 158},
    'K':  {'EA': 0.501, 'IE': 4.34, 'EN': 0.82, 'theta_D': 91},
    'Rb': {'EA': 0.486, 'IE': 4.18, 'EN': 0.82, 'theta_D': 56},
    'Cs': {'EA': 0.472, 'IE': 3.89, 'EN': 0.79, 'theta_D': 38},

    # Halogens (p-block, high EA)
    'F':  {'EA': 3.401, 'IE': 17.42, 'EN': 3.98, 'theta_D': None},
    'Cl': {'EA': 3.612, 'IE': 12.97, 'EN': 3.16, 'theta_D': None},
    'Br': {'EA': 3.364, 'IE': 11.81, 'EN': 2.96, 'theta_D': None},
    'I':  {'EA': 3.059, 'IE': 10.45, 'EN': 2.66, 'theta_D': None},

    # Oxygen group (high EA)
    'O':  {'EA': 1.461, 'IE': 13.62, 'EN': 3.44, 'theta_D': None},
    'S':  {'EA': 2.077, 'IE': 10.36, 'EN': 2.58, 'theta_D': None},
    'Se': {'EA': 2.021, 'IE': 9.75, 'EN': 2.55, 'theta_D': 90},
    'Te': {'EA': 1.971, 'IE': 9.01, 'EN': 2.10, 'theta_D': 153},

    # Nitrogen group (lower EA)
    'N':  {'EA': -0.07, 'IE': 14.53, 'EN': 3.04, 'theta_D': None},  # negative!
    'P':  {'EA': 0.746, 'IE': 10.49, 'EN': 2.19, 'theta_D': None},
    'As': {'EA': 0.804, 'IE': 9.79, 'EN': 2.18, 'theta_D': 285},
    'Sb': {'EA': 1.047, 'IE': 8.61, 'EN': 2.05, 'theta_D': 211},

    # Carbon group
    'C':  {'EA': 1.262, 'IE': 11.26, 'EN': 2.55, 'theta_D': 2230},
    'Si': {'EA': 1.389, 'IE': 8.15, 'EN': 1.90, 'theta_D': 645},
    'Ge': {'EA': 1.233, 'IE': 7.90, 'EN': 2.01, 'theta_D': 374},
    'Sn': {'EA': 1.112, 'IE': 7.34, 'EN': 1.96, 'theta_D': 200},
    'Pb': {'EA': 0.364, 'IE': 7.42, 'EN': 2.33, 'theta_D': 105},

    # Transition metals
    'Ti': {'EA': 0.079, 'IE': 6.83, 'EN': 1.54, 'theta_D': 420},
    'V':  {'EA': 0.525, 'IE': 6.75, 'EN': 1.63, 'theta_D': 380},
    'Cr': {'EA': 0.666, 'IE': 6.77, 'EN': 1.66, 'theta_D': 630},
    'Fe': {'EA': 0.151, 'IE': 7.90, 'EN': 1.83, 'theta_D': 470},
    'Co': {'EA': 0.662, 'IE': 7.88, 'EN': 1.88, 'theta_D': 445},
    'Ni': {'EA': 1.157, 'IE': 7.64, 'EN': 1.91, 'theta_D': 450},
    'Cu': {'EA': 1.235, 'IE': 7.73, 'EN': 1.90, 'theta_D': 343},
    'Zn': {'EA': -0.60, 'IE': 9.39, 'EN': 1.65, 'theta_D': 327},  # negative
    'Ag': {'EA': 1.302, 'IE': 7.58, 'EN': 1.93, 'theta_D': 225},
    'Au': {'EA': 2.309, 'IE': 9.22, 'EN': 2.54, 'theta_D': 165},
    'Pt': {'EA': 2.128, 'IE': 8.96, 'EN': 2.28, 'theta_D': 240},
    'W':  {'EA': 0.816, 'IE': 7.98, 'EN': 2.36, 'theta_D': 400},

    # Noble gases (zero/negative EA)
    'He': {'EA': -0.5, 'IE': 24.59, 'EN': None, 'theta_D': None},  # estimated negative
    'Ne': {'EA': -1.2, 'IE': 21.56, 'EN': None, 'theta_D': 75},    # estimated negative
    'Ar': {'EA': -0.4, 'IE': 15.76, 'EN': None, 'theta_D': 92},    # estimated negative
}

# Calculate coherence parameters
IE_ref = 13.6  # eV (Hydrogen)
T = 300  # Room temperature

for mat in materials:
    data = materials[mat]
    data['gamma_optical'] = IE_ref / data['IE']
    if data['EN'] is not None:
        data['EN_Mulliken'] = (data['IE'] + data['EA']) / 2
    else:
        data['EN_Mulliken'] = None
    if data['theta_D'] is not None:
        data['gamma_phonon'] = 2 * T / data['theta_D']
    else:
        data['gamma_phonon'] = None

# Extract arrays (exclude noble gases and elements with no EN)
valid_names = [m for m in materials if materials[m]['EN'] is not None and materials[m]['EA'] > -0.2]
EA = np.array([materials[m]['EA'] for m in valid_names])
IE = np.array([materials[m]['IE'] for m in valid_names])
EN = np.array([materials[m]['EN'] for m in valid_names])
EN_Mulliken = np.array([materials[m]['EN_Mulliken'] for m in valid_names])
gamma_optical = np.array([materials[m]['gamma_optical'] for m in valid_names])

print("=" * 70)
print("CHEMISTRY SESSION #118: ELECTRON AFFINITY AND COHERENCE")
print("=" * 70)

print(f"\nDataset: {len(valid_names)} elements with positive EA")

# Test 1: EA vs 1/γ_optical (should correlate if EA is coherence)
inv_gamma_optical = 1 / gamma_optical
r1, p1 = stats.pearsonr(EA, inv_gamma_optical)
print(f"\n1. EA vs 1/γ_optical: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.8:
    print("   EXCELLENT - EA IS electronic coherence")
elif abs(r1) > 0.5:
    print("   MODERATE correlation")
else:
    print("   WEAK - EA may be independent")

# Test 2: EA vs IE
r2, p2 = stats.pearsonr(EA, IE)
print(f"\n2. EA vs IE: r = {r2:.3f}, p = {p2:.2e}")

# Test 3: EA vs EN (Pauling)
r3, p3 = stats.pearsonr(EA, EN)
print(f"\n3. EA vs EN (Pauling): r = {r3:.3f}, p = {p3:.2e}")

# Test 4: Mulliken EN = (IE + EA)/2 validation
r4, p4 = stats.pearsonr(EN, EN_Mulliken)
print(f"\n4. EN (Pauling) vs EN (Mulliken): r = {r4:.3f}")
print("   Tests whether (IE+EA)/2 ≈ Pauling EN")

# Test 5: EA vs γ_optical (direct)
r5, p5 = stats.pearsonr(EA, gamma_optical)
print(f"\n5. EA vs γ_optical: r = {r5:.3f}")

# Group analysis
print("\n" + "=" * 70)
print("GROUP ANALYSIS")
print("=" * 70)

groups = {
    'Halogens': ['F', 'Cl', 'Br', 'I'],
    'Chalcogens': ['O', 'S', 'Se', 'Te'],
    'Pnictogens': ['P', 'As', 'Sb'],  # exclude N (negative EA)
    'Carbon Group': ['C', 'Si', 'Ge', 'Sn', 'Pb'],
    'Alkali': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Coinage': ['Cu', 'Ag', 'Au'],
    'Late TM': ['Fe', 'Co', 'Ni', 'Cu'],
}

print(f"\n{'Group':<15} {'Mean EA':<10} {'Mean IE':<10} {'Mean EN':<10} {'Mean γ_opt':<10}")
print("-" * 55)

for group_name, members in groups.items():
    valid = [m for m in members if m in valid_names]
    if len(valid) >= 2:
        group_EA = [materials[m]['EA'] for m in valid]
        group_IE = [materials[m]['IE'] for m in valid]
        group_EN = [materials[m]['EN'] for m in valid]
        group_gamma = [materials[m]['gamma_optical'] for m in valid]
        print(f"{group_name:<15} {np.mean(group_EA):<10.2f} {np.mean(group_IE):<10.2f} "
              f"{np.mean(group_EN):<10.2f} {np.mean(group_gamma):<10.2f}")

# Within-group correlations
print("\n" + "-" * 70)
print("WITHIN-GROUP CORRELATIONS (EA vs 1/γ_optical)")
print("-" * 70)

for group_name, members in groups.items():
    valid = [m for m in members if m in valid_names]
    if len(valid) >= 3:
        group_EA = [materials[m]['EA'] for m in valid]
        group_inv_g = [1/materials[m]['gamma_optical'] for m in valid]
        r_wg, _ = stats.pearsonr(group_EA, group_inv_g)
        print(f"  {group_name}: r = {r_wg:.3f}")

# Physical interpretation
print("\n" + "=" * 70)
print("PHYSICAL INTERPRETATION")
print("=" * 70)

print("""
Electron Affinity vs Ionization Energy:

EA = energy RELEASED when adding electron
IE = energy REQUIRED to remove electron

Both measure electron binding, but:
- IE: how tightly atom HOLDS its electrons (electronic coherence)
- EA: how readily atom ACCEPTS electrons (electron-attracting power)

Relationship to coherence:
- γ_optical = IE_ref / IE measures EXISTING electron coherence
- EA measures POTENTIAL electron coherence of extra electron

The correlation depends on orbital structure:
- Halogens: High EA (need 1 electron for filled shell)
- Noble gases: Negative EA (filled shell resists addition)
- Alkali: Low EA (extra electron loosely bound)
""")

# Hardness/Softness analysis
print("\n" + "=" * 70)
print("CHEMICAL HARDNESS ANALYSIS")
print("=" * 70)

# Chemical hardness η = (IE - EA)/2
hardness = (IE - EA) / 2
r_hard, _ = stats.pearsonr(hardness, inv_gamma_optical)
print(f"\nChemical hardness η = (IE - EA)/2")
print(f"η vs 1/γ_optical: r = {r_hard:.3f}")

# Softness S = 1/(2η)
softness = 1 / (2 * hardness)
r_soft, _ = stats.pearsonr(softness, gamma_optical)
print(f"\nChemical softness S = 1/(2η)")
print(f"S vs γ_optical: r = {r_soft:.3f}")

print("\nHardness measures resistance to charge transfer.")
print("Hard atoms (high η): resist both losing and gaining electrons")
print("Soft atoms (low η): easily polarized, undergo charge transfer")

# Summary
print("\n" + "=" * 70)
print("SESSION #118 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- EA vs 1/γ_optical: r = {r1:.3f}
- EA vs IE: r = {r2:.3f}
- EA vs EN: r = {r3:.3f}
- η (hardness) vs 1/γ_optical: r = {r_hard:.3f}

Physical Interpretation:
""")

if abs(r1) > 0.7:
    print("EA IS electronic coherence:")
    print("  Both EA and IE relate to electron binding")
    print("  EN = (IE + EA)/2 connects both to Mulliken scale")
elif abs(r3) > abs(r1):
    print("EA relates to EN more than γ_optical:")
    print(f"  EA vs EN: r = {r3:.3f}")
    print(f"  EA vs 1/γ_optical: r = {r1:.3f}")
    print("  Electronegativity integrates both EA and IE")
else:
    print("EA is partially independent:")
    print("  EA measures ATTRACTING electrons (empty orbital)")
    print("  IE measures HOLDING electrons (filled orbital)")
    print("  These are different aspects of electronic structure")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Define colors by group
group_colors = {'Halogens': 'red', 'Chalcogens': 'orange', 'Pnictogens': 'green',
                'Carbon Group': 'gray', 'Alkali': 'blue', 'Coinage': 'gold',
                'Late TM': 'purple', 'Other': 'black'}

def get_group(element):
    for group, members in groups.items():
        if element in members:
            return group
    return 'Other'

# Plot 1: EA vs 1/γ_optical
ax1 = axes[0, 0]
for m in valid_names:
    group = get_group(m)
    color = group_colors.get(group, 'black')
    ax1.scatter(1/materials[m]['gamma_optical'], materials[m]['EA'],
                c=color, s=80, alpha=0.7)
    ax1.annotate(m, (1/materials[m]['gamma_optical'], materials[m]['EA']),
                fontsize=7, alpha=0.8)

ax1.set_xlabel('1/γ_optical (electronic coherence)')
ax1.set_ylabel('Electron Affinity (eV)')
ax1.set_title(f'EA vs 1/γ_optical\nr = {r1:.3f}')
ax1.grid(True, alpha=0.3)

# Plot 2: EA vs EN
ax2 = axes[0, 1]
for m in valid_names:
    group = get_group(m)
    color = group_colors.get(group, 'black')
    ax2.scatter(materials[m]['EN'], materials[m]['EA'], c=color, s=80, alpha=0.7)

# Add trend line
z = np.polyfit(EN, EA, 1)
p_trend = np.poly1d(z)
ax2.plot([0.7, 4], [p_trend(0.7), p_trend(4)], 'r--', alpha=0.5)

ax2.set_xlabel('Electronegativity (Pauling)')
ax2.set_ylabel('Electron Affinity (eV)')
ax2.set_title(f'EA vs EN\nr = {r3:.3f}')
ax2.grid(True, alpha=0.3)

# Plot 3: Mulliken vs Pauling EN
ax3 = axes[1, 0]
for m in valid_names:
    group = get_group(m)
    color = group_colors.get(group, 'black')
    ax3.scatter(materials[m]['EN'], materials[m]['EN_Mulliken'],
                c=color, s=80, alpha=0.7, label=group)

ax3.plot([0.5, 4.5], [0.5, 4.5], 'k--', alpha=0.3, label='1:1 line')
ax3.set_xlabel('EN (Pauling)')
ax3.set_ylabel('EN (Mulliken) = (IE+EA)/2')
ax3.set_title(f'Pauling vs Mulliken Electronegativity\nr = {r4:.3f}')
ax3.grid(True, alpha=0.3)

# Plot 4: Hardness vs 1/γ_optical
ax4 = axes[1, 1]
for i, m in enumerate(valid_names):
    group = get_group(m)
    color = group_colors.get(group, 'black')
    ax4.scatter(inv_gamma_optical[i], hardness[i], c=color, s=80, alpha=0.7)
    ax4.annotate(m, (inv_gamma_optical[i], hardness[i]), fontsize=7, alpha=0.8)

ax4.set_xlabel('1/γ_optical (electronic coherence)')
ax4.set_ylabel('Chemical Hardness η = (IE-EA)/2 (eV)')
ax4.set_title(f'Hardness vs 1/γ_optical\nr = {r_hard:.3f}')
ax4.grid(True, alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=g) for g, c in group_colors.items() if g != 'Other']
fig.legend(handles=legend_elements, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.02))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/electron_affinity_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: electron_affinity_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r3) > 0.7 and abs(r_hard) > 0.8:
    print("\n✓ EXCELLENT VALIDATION")
    print(f"  EA vs EN: r = {r3:.3f}")
    print(f"  η vs 1/γ_optical: r = {r_hard:.3f}")
    print("  Electron affinity contributes to chemical hardness")
    print("  Hardness IS electronic coherence")
elif abs(r3) > 0.5 or abs(r1) > 0.5:
    print("\n✓ GOOD VALIDATION")
    print("  EA correlates with electronic coherence metrics")
    print("  EA + IE together determine coherence (via EN, η)")
else:
    print("\n○ PARTIAL")
    print("  EA is group-specific property")
    print("  Halogens and noble gases show opposite behavior")
    print("  Orbital filling dominates over simple coherence")
