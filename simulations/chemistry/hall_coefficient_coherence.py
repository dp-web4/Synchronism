#!/usr/bin/env python3
"""
Chemistry Session #102: Hall Coefficient and Coherence
=======================================================

Test whether Hall coefficient R_H relates to coherence parameters.

Physical reasoning:
- R_H = 1/(n×e) for free electrons (simple metals)
- Sign: negative for electron-like, positive for hole-like
- Magnitude: inversely proportional to carrier density
- Deviations from free-electron indicate d-band effects

Connection to coherence:
- Hall effect measures carrier response to magnetic field
- Multi-band effects can cause anomalous Hall
- Does R_H correlate with γ_electron or E_F?

Hypothesis:
- |R_H| should correlate with 1/n (carrier density)
- Deviations from free electron indicate coherence effects
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Dataset: Hall coefficients at room temperature
# R_H: Hall coefficient (10^-11 m³/C) - negative for electrons, positive for holes
# n: Nominal electron density (10^28 m^-3)
# θ_D: Debye temperature (K)
# λ_ep: Electron-phonon coupling

materials = {
    # Element, R_H (10^-11 m³/C), n (10^28/m³), θ_D (K), λ_ep, valence
    # Alkali metals (nearly free electron)
    'Li': (-17.0, 4.70, 344, 0.45, 1),
    'Na': (-24.8, 2.65, 158, 0.16, 1),
    'K': (-44.5, 1.40, 91, 0.13, 1),
    'Rb': (-58.9, 1.15, 56, 0.12, 1),
    'Cs': (-73.3, 0.91, 38, 0.13, 1),

    # Noble metals (s-d hybridization)
    'Cu': (-5.5, 8.47, 343, 0.13, 1),
    'Ag': (-9.0, 5.86, 225, 0.12, 1),
    'Au': (-7.2, 5.90, 165, 0.16, 1),

    # Divalent (two s-electrons)
    'Mg': (-9.4, 8.61, 400, 0.30, 2),
    'Zn': (+4.1, 13.2, 327, 0.45, 2),  # POSITIVE - holes!
    'Cd': (+6.0, 9.28, 209, 0.38, 2),  # POSITIVE - holes!
    'Ca': (-17.8, 4.61, 230, 0.35, 2),
    'Ba': (-23.7, 3.15, 110, 0.40, 2),

    # Trivalent
    'Al': (-3.5, 18.1, 428, 0.43, 3),
    'Ga': (-8.0, 15.3, 325, 0.40, 3),
    'In': (-2.4, 11.5, 108, 0.80, 3),

    # Tetravalent
    'Pb': (-0.09, 13.2, 105, 1.55, 4),  # Nearly zero - compensated
    'Sn': (-0.2, 14.8, 200, 0.72, 4),

    # Transition metals (complex band structure)
    'Fe': (+100.0, 17.0, 470, 0.32, 2),  # LARGE POSITIVE - anomalous!
    'Co': (+36.0, 18.0, 445, 0.45, 2),   # POSITIVE
    'Ni': (-6.0, 18.1, 450, 0.49, 2),    # Negative
    'Cr': (+36.0, 8.5, 630, 0.35, 1),    # POSITIVE
    'Pd': (-69.0, 8.5, 274, 0.42, 0),    # Large negative
    'Pt': (-21.0, 6.6, 240, 0.66, 0),    # Negative
    'W': (+8.0, 12.1, 400, 0.28, 2),     # Positive
    'Mo': (+18.0, 8.5, 450, 0.41, 1),    # Positive
    'Ti': (-24.0, 5.5, 420, 0.38, 2),
    'V': (+7.6, 10.0, 380, 0.60, 2),

    # Heavy elements
    'Bi': (-50000.0, 0.141, 119, 0.20, 5),  # SEMIMETAL - very large!
    'As': (+50000.0, 0.20, 282, 0.25, 5),   # SEMIMETAL - large positive
}

# Extract data
names = list(materials.keys())
R_H = np.array([materials[m][0] for m in names])  # 10^-11 m³/C
n = np.array([materials[m][1] for m in names])  # 10^28/m³
theta_D = np.array([materials[m][2] for m in names])  # K
lambda_ep = np.array([materials[m][3] for m in names])
valence = np.array([materials[m][4] for m in names])

# Exclude semimetals for main analysis (extreme outliers)
semimetal_mask = np.abs(R_H) > 1000
normal_mask = ~semimetal_mask

# Calculate derived quantities
T = 300  # K
gamma_phonon = 2 * T / theta_D
gamma_electron = 2 * lambda_ep / (1 + lambda_ep)

# Free electron prediction: R_H = 1/(n×e) = 6.24×10^-29 / n (in SI)
# In our units (10^-11 m³/C with n in 10^28 m^-3):
# R_H_free = -1/(n_SI × e) = -1/(n × 10^28 × 1.6×10^-19) = -6.24×10^-10 / n
# Converting: R_H (10^-11) = -62.4 / n (for n in 10^28)
R_H_free = -62.4 / n  # Free electron prediction (negative for electrons)

print("="*70)
print("CHEMISTRY SESSION #102: HALL COEFFICIENT AND COHERENCE")
print("="*70)
print(f"\nDataset: {len(names)} metals ({np.sum(normal_mask)} normal, {np.sum(semimetal_mask)} semimetals)")
print(f"R_H range: {R_H[normal_mask].min():.1f} to {R_H[normal_mask].max():.1f} × 10^-11 m³/C")
print(f"Carrier density range: {n[normal_mask].min():.2f} - {n[normal_mask].max():.2f} × 10^28 m^-3")

# Analysis 1: Sign of R_H
print("\n" + "="*70)
print("ANALYSIS 1: CARRIER SIGN")
print("="*70)

n_negative = np.sum(R_H[normal_mask] < 0)
n_positive = np.sum(R_H[normal_mask] > 0)
print(f"\nNormal metals:")
print(f"  Electron-like (R_H < 0): {n_negative}")
print(f"  Hole-like (R_H > 0): {n_positive}")

print("\nHole-like metals:")
for name in names:
    if materials[name][0] > 0 and abs(materials[name][0]) < 1000:
        print(f"  {name}: R_H = +{materials[name][0]:.1f} × 10^-11 m³/C")

# Analysis 2: Free electron comparison
print("\n" + "="*70)
print("ANALYSIS 2: FREE ELECTRON MODEL")
print("="*70)

# Ratio of observed to free electron
ratio = R_H / R_H_free

print("\n|R_H|_obs / |R_H|_free ratios (normal metals):")
for name in names:
    if not semimetal_mask[names.index(name)]:
        i = names.index(name)
        r = abs(R_H[i] / R_H_free[i])
        sign = '+' if R_H[i] > 0 else '-'
        if r > 2 or r < 0.5 or R_H[i] > 0:
            print(f"  {name}: {r:.2f}× ({sign})")

# For monovalent alkalis (best free electron)
alkali_mask = np.array([n in ['Li', 'Na', 'K', 'Rb', 'Cs'] for n in names])
r_alkali, p_alkali = stats.pearsonr(R_H[alkali_mask], R_H_free[alkali_mask])
print(f"\nAlkali metals: R_H_obs vs R_H_free, r = {r_alkali:.3f}")

# Analysis 3: Hall coefficient vs carrier density
print("\n" + "="*70)
print("ANALYSIS 3: |R_H| vs 1/n")
print("="*70)

# Use only negative R_H metals (electron-like)
electron_mask = (R_H < 0) & normal_mask
r_density, p_density = stats.pearsonr(np.abs(R_H[electron_mask]), 1/n[electron_mask])
print(f"\nElectron-like metals: |R_H| vs 1/n, r = {r_density:.3f}")

# Analysis 4: Deviations and coherence
print("\n" + "="*70)
print("ANALYSIS 4: DEVIATIONS FROM FREE ELECTRON")
print("="*70)

# Define Hall enhancement factor
hall_factor = np.abs(R_H / R_H_free)  # For negative R_H, this should be ~1

print("\nHall factor |R_H_obs|/|R_H_free| analysis:")

# Correlation with λ_ep (normal electron-like metals)
normal_electron_mask = electron_mask
r_lambda, p_lambda = stats.pearsonr(hall_factor[normal_electron_mask], lambda_ep[normal_electron_mask])
print(f"Hall factor vs λ_ep: r = {r_lambda:.3f}")

r_gamma_e, p_gamma_e = stats.pearsonr(hall_factor[normal_electron_mask], gamma_electron[normal_electron_mask])
print(f"Hall factor vs γ_electron: r = {r_gamma_e:.3f}")

r_gamma_p, p_gamma_p = stats.pearsonr(hall_factor[normal_electron_mask], gamma_phonon[normal_electron_mask])
print(f"Hall factor vs γ_phonon: r = {r_gamma_p:.3f}")

# Analysis 5: Material class comparison
print("\n" + "="*70)
print("ANALYSIS 5: MATERIAL CLASS COMPARISON")
print("="*70)

alkali = ['Li', 'Na', 'K', 'Rb', 'Cs']
noble = ['Cu', 'Ag', 'Au']
divalent = ['Mg', 'Zn', 'Cd', 'Ca', 'Ba']
trivalent = ['Al', 'Ga', 'In']
tetravalent = ['Pb', 'Sn']
transition = ['Fe', 'Co', 'Ni', 'Cr', 'Pd', 'Pt', 'W', 'Mo', 'Ti', 'V']

print("\n| Class | Mean |R_H| | Mean n | Hall Factor | % Positive |")
print("|-------|----------|--------|-------------|------------|")

for class_name, class_list in [('Alkali', alkali), ('Noble', noble),
                                ('Divalent', divalent), ('Trivalent', trivalent),
                                ('Transition', transition)]:
    mask = np.array([name in class_list for name in names]) & normal_mask
    if np.sum(mask) >= 2:
        mean_RH = np.mean(np.abs(R_H[mask]))
        mean_n = np.mean(n[mask])
        mean_hf = np.mean(hall_factor[mask])
        pct_pos = 100 * np.sum(R_H[mask] > 0) / np.sum(mask)
        print(f"| {class_name:10} | {mean_RH:8.1f} | {mean_n:6.2f} | {mean_hf:11.2f} | {pct_pos:10.0f}% |")

# Analysis 6: Anomalous Hall effect
print("\n" + "="*70)
print("ANALYSIS 6: ANOMALOUS HALL EFFECT")
print("="*70)

print("""
ANOMALOUS HALL EFFECT occurs when:
1. R_H has WRONG SIGN (positive for electron metals)
2. R_H has WRONG MAGNITUDE (>>1 or <<1 × free electron)

Physical origins:
- Multi-band effects: different bands have different signs
- Magnetic scattering: spin-orbit gives extra contribution
- Band structure: hole pockets in Fermi surface

Examples from this dataset:
""")

for name in ['Fe', 'Co', 'Zn', 'Cd', 'Mo', 'Cr']:
    i = names.index(name)
    print(f"  {name}: R_H = {R_H[i]:+.1f}, ratio = {hall_factor[i]:.2f}×")

# Analysis 7: Connection to coherence interpretation
print("\n" + "="*70)
print("ANALYSIS 7: COHERENCE INTERPRETATION")
print("="*70)

print("""
Hall coefficient is primarily determined by BAND STRUCTURE:

1. CARRIER TYPE (sign)
   - Simple metals: electron-like (R_H < 0)
   - Divalent/transition: often hole-like (R_H > 0)

2. CARRIER DENSITY (magnitude)
   - |R_H| ∝ 1/n (inverse carrier density)
   - This is EXTENSIVE property like E_F (Session #100)

3. SCATTERING EFFECTS
   - Hall mobility μ_H = R_H × σ
   - If σ ∝ 1/γ_electron, then μ_H should correlate with γ

Let's test: Hall mobility μ_H = R_H × σ (not exactly, but related)
""")

# Create visualization
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Plot 1: R_H vs 1/n
ax1 = axes[0, 0]
colors = {'alkali': 'red', 'noble': 'gold', 'divalent': 'blue',
          'trivalent': 'green', 'transition': 'purple'}
class_map = {n: 'alkali' for n in alkali}
class_map.update({n: 'noble' for n in noble})
class_map.update({n: 'divalent' for n in divalent})
class_map.update({n: 'trivalent' for n in trivalent})
class_map.update({n: 'transition' for n in transition})

for name in names:
    if not semimetal_mask[names.index(name)]:
        i = names.index(name)
        mat_class = class_map.get(name, 'other')
        c = colors.get(mat_class, 'gray')
        ax1.scatter(1/n[i], R_H[i], c=c, s=80, edgecolors='black', linewidth=0.5)
        if R_H[i] > 20 or R_H[i] < -50 or abs(R_H[i] - R_H_free[i]) > 30:
            ax1.annotate(name, (1/n[i], R_H[i]), fontsize=8, xytext=(3, 3),
                        textcoords='offset points')

# Free electron line
x_range = np.linspace(0.05, 1.2, 100)
ax1.plot(x_range, -62.4 * x_range, 'k--', linewidth=2, label='Free electron')
ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax1.set_xlabel('1/n (10^-28 m³)', fontsize=11)
ax1.set_ylabel('R_H (10^-11 m³/C)', fontsize=11)
ax1.set_title('Hall Coefficient vs Inverse Carrier Density', fontsize=13)
ax1.legend()

# Plot 2: Observed vs Free electron
ax2 = axes[0, 1]
for name in names:
    if not semimetal_mask[names.index(name)]:
        i = names.index(name)
        mat_class = class_map.get(name, 'other')
        c = colors.get(mat_class, 'gray')
        ax2.scatter(R_H_free[i], R_H[i], c=c, s=80, edgecolors='black', linewidth=0.5)
ax2.plot([-80, 0], [-80, 0], 'k--', linewidth=2, label='1:1')
ax2.set_xlabel('R_H (free electron)', fontsize=11)
ax2.set_ylabel('R_H (observed)', fontsize=11)
ax2.set_title(f'Free Electron Model Test (alkali r = {r_alkali:.3f})', fontsize=13)
ax2.legend()

# Plot 3: Hall factor distribution
ax3 = axes[0, 2]
hf_electron = hall_factor[electron_mask]
ax3.hist(hf_electron, bins=15, color='steelblue', edgecolor='black', alpha=0.7)
ax3.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='Free electron')
ax3.set_xlabel('|R_H_obs| / |R_H_free|', fontsize=11)
ax3.set_ylabel('Count', fontsize=11)
ax3.set_title('Hall Factor Distribution (electron-like)', fontsize=13)
ax3.legend()

# Plot 4: Hall factor vs γ_electron
ax4 = axes[1, 0]
for name in names:
    i = names.index(name)
    if normal_mask[i] and R_H[i] < 0:
        mat_class = class_map.get(name, 'other')
        c = colors.get(mat_class, 'gray')
        ax4.scatter(gamma_electron[i], hall_factor[i], c=c, s=80,
                   edgecolors='black', linewidth=0.5)
ax4.set_xlabel('γ_electron = 2λ/(1+λ)', fontsize=11)
ax4.set_ylabel('Hall Factor |R_H_obs|/|R_H_free|', fontsize=11)
ax4.set_title(f'Hall Factor vs γ_electron (r = {r_gamma_e:.3f})', fontsize=13)
ax4.axhline(y=1.0, color='red', linestyle='--', alpha=0.5)

# Plot 5: Sign classification
ax5 = axes[1, 1]
for name in names:
    if not semimetal_mask[names.index(name)]:
        i = names.index(name)
        mat_class = class_map.get(name, 'other')
        c = colors.get(mat_class, 'gray')
        marker = 'o' if R_H[i] < 0 else '^'
        ax5.scatter(n[i], np.abs(R_H[i]), c=c, s=100, marker=marker,
                   edgecolors='black', linewidth=0.5)
# Custom legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                          markersize=10, label='Electron-like (R_H < 0)'),
                   Line2D([0], [0], marker='^', color='w', markerfacecolor='gray',
                          markersize=10, label='Hole-like (R_H > 0)')]
ax5.legend(handles=legend_elements)
ax5.set_xlabel('Carrier density n (10^28 m^-3)', fontsize=11)
ax5.set_ylabel('|R_H| (10^-11 m³/C)', fontsize=11)
ax5.set_title('Carrier Sign Classification', fontsize=13)
ax5.set_yscale('log')

# Plot 6: Class averages
ax6 = axes[1, 2]
class_data = []
for class_name, class_list in [('Alkali', alkali), ('Noble', noble),
                                ('Divalent', divalent), ('Trivalent', trivalent),
                                ('Transition', transition)]:
    mask = np.array([name in class_list for name in names]) & normal_mask
    if np.sum(mask) >= 2:
        class_data.append((class_name, np.mean(hall_factor[mask]),
                          np.mean(gamma_electron[mask]), colors.get(class_name.lower(), 'gray')))

class_names_plot = [d[0] for d in class_data]
class_hf = [d[1] for d in class_data]
class_ge = [d[2] for d in class_data]
class_colors_plot = [d[3] for d in class_data]

for i, (name, hf, ge, c) in enumerate(class_data):
    ax6.scatter(ge, hf, s=200, c=c, edgecolors='black', linewidth=1.5)
    ax6.annotate(name, (ge, hf), textcoords="offset points", xytext=(10, 5),
                fontsize=11, fontweight='bold')
ax6.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='Free electron')
ax6.set_xlabel('Mean γ_electron', fontsize=11)
ax6.set_ylabel('Mean Hall Factor', fontsize=11)
ax6.set_title('Material Class Averages', fontsize=13)
ax6.legend()

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/hall_coefficient_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "="*70)
print("SESSION #102 SUMMARY")
print("="*70)

print(f"""
HALL COEFFICIENT AND COHERENCE

Key Results:
1. Alkali metals: R_H_obs vs R_H_free, r = {r_alkali:.3f} (FREE ELECTRON WORKS)
2. Hall factor vs γ_electron: r = {r_gamma_e:.3f} (WEAK)
3. Hall factor vs λ_ep: r = {r_lambda:.3f} (WEAK)
4. Many transition metals have POSITIVE R_H (hole-like)

Physical Interpretation:

1. **Hall Effect is BAND STRUCTURE Dominated**
   R_H determined by Fermi surface topology:
   - Electron pockets → negative R_H
   - Hole pockets → positive R_H
   - Multi-band → mixed signs, anomalous R_H

2. **Alkali Metals are FREE ELECTRON**
   R_H_obs/R_H_free ≈ 1 for Li, Na, K, Rb, Cs
   Single s-band, spherical Fermi surface.

3. **Transition Metals are ANOMALOUS**
   Fe: R_H = +100 (large positive!)
   Co: R_H = +36
   These have hole-dominated transport despite electron count.

4. **Coherence is SECONDARY**
   Hall factor shows WEAK correlation with γ_electron (r = {r_gamma_e:.3f})
   Band structure dominates over scattering effects.

Framework Classification:

Hall coefficient R_H is BAND STRUCTURE DOMINATED:
- Primary: Fermi surface topology (extensive)
- Secondary: Carrier scattering (γ_electron)

Compare to:
- Conductivity σ: γ_electron HELPS (Session #86)
- Sommerfeld γ_S: λ_ep HELPS within class (Session #101)
- Hall R_H: Band structure DOMINATES

This is similar to:
- Thermionic emission (#98): φ dominates, weak γ
- Magnetostriction (#94): SOC dominates, weak γ

Coherence matters for TRANSPORT EFFICIENCY, not CARRIER TYPE.

The Hall coefficient tells you WHAT carries current.
The coherence parameter tells you HOW WELL they carry it.
""")

print("\n[Plot saved to hall_coefficient_coherence.png]")
