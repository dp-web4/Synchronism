#!/usr/bin/env python3
"""
Chemistry Session #120: Bulk Modulus and Coherence (Comprehensive)

Follow-up to:
- Session #110: Elastic moduli (G, E, B vs coherence)
- Session #113: Compressibility (κ_T = 1/B vs γ_phonon)

Comprehensive analysis of bulk modulus B as coherence parameter.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Comprehensive bulk modulus data (GPa)
# Sources: Kittel, ASM Handbook, various experimental compilations
materials = {
    # Diamond and carbides (ultra-high B)
    'Diamond': {'B': 443, 'theta_D': 2230, 'V_a': 5.7, 'E_coh': 7.37, 'IE': 11.26},
    'SiC':     {'B': 225, 'theta_D': 1200, 'V_a': 10.5, 'E_coh': 6.2, 'IE': None},
    'BN':      {'B': 400, 'theta_D': 1900, 'V_a': 5.9, 'E_coh': 6.5, 'IE': None},

    # Refractory metals
    'W':  {'B': 310, 'theta_D': 400, 'V_a': 15.8, 'E_coh': 8.90, 'IE': 7.98},
    'Mo': {'B': 230, 'theta_D': 450, 'V_a': 15.6, 'E_coh': 6.82, 'IE': 7.09},
    'Ta': {'B': 200, 'theta_D': 240, 'V_a': 18.0, 'E_coh': 8.10, 'IE': 7.55},
    'Nb': {'B': 170, 'theta_D': 275, 'V_a': 18.0, 'E_coh': 7.57, 'IE': 6.76},

    # 3d transition metals
    'Fe': {'B': 170, 'theta_D': 470, 'V_a': 11.8, 'E_coh': 4.28, 'IE': 7.90},
    'Co': {'B': 180, 'theta_D': 445, 'V_a': 11.1, 'E_coh': 4.39, 'IE': 7.88},
    'Ni': {'B': 180, 'theta_D': 450, 'V_a': 11.0, 'E_coh': 4.44, 'IE': 7.64},
    'Cr': {'B': 160, 'theta_D': 630, 'V_a': 12.0, 'E_coh': 4.10, 'IE': 6.77},
    'V':  {'B': 162, 'theta_D': 380, 'V_a': 13.8, 'E_coh': 5.31, 'IE': 6.75},
    'Ti': {'B': 110, 'theta_D': 420, 'V_a': 17.7, 'E_coh': 4.85, 'IE': 6.83},
    'Cu': {'B': 140, 'theta_D': 343, 'V_a': 11.8, 'E_coh': 3.49, 'IE': 7.73},
    'Zn': {'B': 72, 'theta_D': 327, 'V_a': 15.2, 'E_coh': 1.35, 'IE': 9.39},

    # Noble metals
    'Ag': {'B': 100, 'theta_D': 225, 'V_a': 17.1, 'E_coh': 2.95, 'IE': 7.58},
    'Au': {'B': 180, 'theta_D': 165, 'V_a': 17.0, 'E_coh': 3.81, 'IE': 9.22},
    'Pt': {'B': 230, 'theta_D': 240, 'V_a': 15.1, 'E_coh': 5.84, 'IE': 8.96},

    # Alkali metals
    'Li': {'B': 11.0, 'theta_D': 344, 'V_a': 21.3, 'E_coh': 1.63, 'IE': 5.39},
    'Na': {'B': 6.3, 'theta_D': 158, 'V_a': 39.5, 'E_coh': 1.11, 'IE': 5.14},
    'K':  {'B': 3.1, 'theta_D': 91, 'V_a': 75.6, 'E_coh': 0.93, 'IE': 4.34},
    'Rb': {'B': 2.5, 'theta_D': 56, 'V_a': 93.0, 'E_coh': 0.85, 'IE': 4.18},
    'Cs': {'B': 1.6, 'theta_D': 38, 'V_a': 115, 'E_coh': 0.80, 'IE': 3.89},

    # Alkaline earth
    'Be': {'B': 130, 'theta_D': 1440, 'V_a': 8.1, 'E_coh': 3.32, 'IE': 9.32},
    'Mg': {'B': 45, 'theta_D': 400, 'V_a': 23.2, 'E_coh': 1.51, 'IE': 7.65},
    'Ca': {'B': 15, 'theta_D': 230, 'V_a': 43.2, 'E_coh': 1.84, 'IE': 6.11},

    # Simple metals
    'Al': {'B': 76, 'theta_D': 428, 'V_a': 16.6, 'E_coh': 3.39, 'IE': 5.99},
    'Pb': {'B': 46, 'theta_D': 105, 'V_a': 30.3, 'E_coh': 2.03, 'IE': 7.42},
    'Sn': {'B': 58, 'theta_D': 200, 'V_a': 27.0, 'E_coh': 3.14, 'IE': 7.34},

    # Semiconductors
    'Si': {'B': 98, 'theta_D': 645, 'V_a': 20.0, 'E_coh': 4.63, 'IE': 8.15},
    'Ge': {'B': 75, 'theta_D': 374, 'V_a': 22.6, 'E_coh': 3.85, 'IE': 7.90},
    'GaAs': {'B': 75, 'theta_D': 344, 'V_a': 22.7, 'E_coh': 3.3, 'IE': None},

    # Ionic crystals
    'NaCl': {'B': 25, 'theta_D': 321, 'V_a': 22.3, 'E_coh': 3.28, 'IE': None},
    'KCl':  {'B': 18, 'theta_D': 235, 'V_a': 31.5, 'E_coh': 3.40, 'IE': None},
    'MgO':  {'B': 160, 'theta_D': 946, 'V_a': 9.3, 'E_coh': 5.0, 'IE': None},
}

# Calculate coherence parameters
T = 300  # Room temperature
IE_ref = 13.6  # eV

for mat in materials:
    data = materials[mat]
    data['gamma_phonon'] = 2 * T / data['theta_D']
    data['inv_gamma_phonon'] = 1 / data['gamma_phonon']
    data['kappa_T'] = 1 / data['B']  # Compressibility
    data['log_B'] = np.log10(data['B'])
    if data['IE'] is not None:
        data['gamma_optical'] = IE_ref / data['IE']
    else:
        data['gamma_optical'] = None

# Extract arrays
names = list(materials.keys())
B = np.array([materials[m]['B'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
inv_gamma_phonon = np.array([materials[m]['inv_gamma_phonon'] for m in names])
V_a = np.array([materials[m]['V_a'] for m in names])
E_coh = np.array([materials[m]['E_coh'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #120: BULK MODULUS AND COHERENCE (COMPREHENSIVE)")
print("=" * 70)

print(f"\nDataset: {len(names)} materials")
print(f"B range: {np.min(B):.1f} - {np.max(B):.0f} GPa ({np.max(B)/np.min(B):.0f}× range)")

# Test 1: B vs 1/γ_phonon (main coherence test)
r1, p1 = stats.pearsonr(B, inv_gamma_phonon)
print(f"\n1. B vs 1/γ_phonon: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.8:
    print("   EXCELLENT - Bulk modulus IS phonon coherence")
elif abs(r1) > 0.5:
    print("   MODERATE correlation")

# Test 2: B vs θ_D² (Debye model predicts B ∝ θ_D²)
theta_D_sq = theta_D**2
r2, p2 = stats.pearsonr(B, theta_D_sq)
print(f"\n2. B vs θ_D²: r = {r2:.3f}")

# Test 3: B vs 1/V_a (should correlate if B ~ E_coh/V)
inv_V_a = 1 / V_a
r3, p3 = stats.pearsonr(B, inv_V_a)
print(f"\n3. B vs 1/V_a: r = {r3:.3f}")

# Test 4: B vs E_coh/V_a (cohesive energy density)
E_coh_density = E_coh / V_a
r4, p4 = stats.pearsonr(B, E_coh_density)
print(f"\n4. B vs E_coh/V_a: r = {r4:.3f}")
if abs(r4) > 0.8:
    print("   EXCELLENT - B IS cohesive energy density")

# Test 5: B × V_a (energy) vs E_coh
BV = B * V_a / 160.2  # Convert GPa·Å³ to eV (1 GPa·Å³ = 0.00624 eV)
r5, p5 = stats.pearsonr(BV, E_coh)
print(f"\n5. B×V_a vs E_coh: r = {r5:.3f}")

# Test 6: log-log relationship
log_B = np.log10(B)
log_theta = np.log10(theta_D)
r6, p6 = stats.pearsonr(log_B, log_theta)
slope6, intercept6, _, _, _ = stats.linregress(log_theta, log_B)
print(f"\n6. log(B) vs log(θ_D): r = {r6:.3f}")
print(f"   Power law: B ∝ θ_D^{slope6:.2f}")

# Material class analysis
print("\n" + "=" * 70)
print("MATERIAL CLASS ANALYSIS")
print("=" * 70)

classes = {
    'Ceramics': ['Diamond', 'SiC', 'BN', 'MgO'],
    'Refractory': ['W', 'Mo', 'Ta', 'Nb'],
    '3d TM': ['Fe', 'Co', 'Ni', 'Cr', 'V', 'Ti', 'Cu', 'Zn'],
    'Noble': ['Ag', 'Au', 'Pt'],
    'Alkali': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Simple': ['Al', 'Pb', 'Sn', 'Mg'],
    'Semiconductors': ['Si', 'Ge', 'GaAs'],
    'Ionic': ['NaCl', 'KCl'],
}

print(f"\n{'Class':<15} {'Mean B':<10} {'Mean θ_D':<10} {'Mean γ_ph':<10} {'Mean V_a':<10}")
print("-" * 60)

class_stats = {}
for class_name, members in classes.items():
    valid = [m for m in members if m in names]
    class_B = [materials[m]['B'] for m in valid]
    class_theta = [materials[m]['theta_D'] for m in valid]
    class_gamma = [materials[m]['gamma_phonon'] for m in valid]
    class_V = [materials[m]['V_a'] for m in valid]
    class_stats[class_name] = {
        'B': np.mean(class_B),
        'theta_D': np.mean(class_theta),
        'gamma_phonon': np.mean(class_gamma),
        'V_a': np.mean(class_V)
    }
    print(f"{class_name:<15} {np.mean(class_B):<10.0f} {np.mean(class_theta):<10.0f} "
          f"{np.mean(class_gamma):<10.2f} {np.mean(class_V):<10.1f}")

# Within-class correlations
print("\n" + "-" * 70)
print("WITHIN-CLASS CORRELATIONS (B vs 1/γ_phonon)")
print("-" * 70)

for class_name, members in classes.items():
    valid = [m for m in members if m in names]
    if len(valid) >= 3:
        class_B = [materials[m]['B'] for m in valid]
        class_inv_g = [1/materials[m]['gamma_phonon'] for m in valid]
        r_wc, _ = stats.pearsonr(class_B, class_inv_g)
        print(f"  {class_name}: r = {r_wc:.3f}")

# Physical model comparison
print("\n" + "=" * 70)
print("PHYSICAL MODEL COMPARISON")
print("=" * 70)

# Model 1: B ∝ 1/γ_phonon (simple coherence)
model1 = inv_gamma_phonon
r_m1, _ = stats.pearsonr(B, model1)
print(f"\n1. B ∝ 1/γ_phonon: r = {r_m1:.3f}")

# Model 2: B ∝ θ_D²/V_a (Debye + volume)
model2 = theta_D**2 / V_a
r_m2, _ = stats.pearsonr(B, model2)
print(f"2. B ∝ θ_D²/V_a: r = {r_m2:.3f}")

# Model 3: B ∝ E_coh/V_a (cohesive energy density)
model3 = E_coh / V_a
r_m3, _ = stats.pearsonr(B, model3)
print(f"3. B ∝ E_coh/V_a: r = {r_m3:.3f}")

# Model 4: B ∝ (1/γ)²/V_a (coherence² × density)
model4 = inv_gamma_phonon**2 / V_a
r_m4, _ = stats.pearsonr(B, model4)
print(f"4. B ∝ (1/γ)²/V_a: r = {r_m4:.3f}")

# Fermi gas model: B ∝ n^(5/3) for free electrons
n = 1 / V_a  # Simple electron density proxy
model5 = n**(5/3)
r_m5, _ = stats.pearsonr(B, model5)
print(f"5. B ∝ n^(5/3) (Fermi gas): r = {r_m5:.3f}")

# Summary
print("\n" + "=" * 70)
print("SESSION #120 SUMMARY")
print("=" * 70)

best_model = max([(r_m1, '1/γ'), (r_m2, 'θ_D²/V_a'), (r_m3, 'E_coh/V_a'),
                  (r_m4, '(1/γ)²/V_a'), (r_m5, 'n^(5/3)')], key=lambda x: abs(x[0]))

print(f"""
Key Results:
- B vs 1/γ_phonon: r = {r1:.3f}
- B vs θ_D²: r = {r2:.3f}
- B vs E_coh/V_a: r = {r4:.3f}
- B ∝ θ_D^{slope6:.2f}: r = {r6:.3f}

Best model: B ∝ {best_model[1]}: r = {best_model[0]:.3f}

Physical Interpretation:
""")

if abs(r4) >= abs(r1) and abs(r4) > 0.8:
    print("B IS cohesive energy density:")
    print(f"  B vs E_coh/V_a: r = {r4:.3f}")
    print("  B = -V × dP/dV ∝ E_coh/V")
    print("  Resistance to compression = bonding energy per volume")
    print("\n  Coherence enters via E_coh:")
    print("    E_coh ∝ k × r² (bond energy)")
    print("    k ∝ θ_D² (spring constant)")
    print("    θ_D ∝ 1/γ_phonon (coherence)")
elif abs(r1) > 0.7:
    print("B IS coherence-dependent:")
    print(f"  B vs 1/γ_phonon: r = {r1:.3f}")
    print("  Stiffer lattices (low γ) have higher B")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: B vs 1/γ_phonon
ax1 = axes[0, 0]
colors = {'Ceramics': 'purple', 'Refractory': 'red', '3d TM': 'blue',
          'Noble': 'gold', 'Alkali': 'green', 'Simple': 'orange',
          'Semiconductors': 'gray', 'Ionic': 'cyan'}

for class_name, members in classes.items():
    valid = [m for m in members if m in names]
    class_B = [materials[m]['B'] for m in valid]
    class_inv_g = [1/materials[m]['gamma_phonon'] for m in valid]
    ax1.scatter(class_inv_g, class_B, c=colors[class_name],
                label=class_name, s=80, alpha=0.7)

ax1.set_xlabel('1/γ_phonon (coherence)')
ax1.set_ylabel('Bulk Modulus B (GPa)')
ax1.set_title(f'B vs 1/γ_phonon\nr = {r1:.3f}')
ax1.legend(loc='best', fontsize=7)
ax1.grid(True, alpha=0.3)

# Plot 2: B vs E_coh/V_a
ax2 = axes[0, 1]
for class_name, members in classes.items():
    valid = [m for m in members if m in names]
    class_B = [materials[m]['B'] for m in valid]
    class_EV = [materials[m]['E_coh']/materials[m]['V_a'] for m in valid]
    ax2.scatter(class_EV, class_B, c=colors[class_name],
                label=class_name, s=80, alpha=0.7)

# Add trend line
z = np.polyfit(E_coh_density, B, 1)
p = np.poly1d(z)
ax2.plot([0, 1.5], [p(0), p(1.5)], 'r--', alpha=0.5)

ax2.set_xlabel('E_coh/V_a (eV/Å³)')
ax2.set_ylabel('Bulk Modulus B (GPa)')
ax2.set_title(f'B vs E_coh/V_a\nr = {r4:.3f}')
ax2.legend(loc='best', fontsize=7)
ax2.grid(True, alpha=0.3)

# Plot 3: log(B) vs log(θ_D)
ax3 = axes[1, 0]
for class_name, members in classes.items():
    valid = [m for m in members if m in names]
    class_log_B = [np.log10(materials[m]['B']) for m in valid]
    class_log_theta = [np.log10(materials[m]['theta_D']) for m in valid]
    ax3.scatter(class_log_theta, class_log_B, c=colors[class_name],
                label=class_name, s=80, alpha=0.7)
    for i, m in enumerate(valid):
        ax3.annotate(m, (class_log_theta[i], class_log_B[i]), fontsize=6, alpha=0.8)

# Add power law fit
ax3.plot([1.5, 3.5], [slope6*1.5+intercept6, slope6*3.5+intercept6], 'r--',
         alpha=0.5, label=f'B ∝ θ_D^{slope6:.1f}')

ax3.set_xlabel('log₁₀(θ_D)')
ax3.set_ylabel('log₁₀(B)')
ax3.set_title(f'Power Law: B ∝ θ_D^{slope6:.2f}\nr = {r6:.3f}')
ax3.legend(loc='best', fontsize=7)
ax3.grid(True, alpha=0.3)

# Plot 4: Model comparison
ax4 = axes[1, 1]
model_names = ['1/γ', 'θ_D²/V_a', 'E_coh/V_a', '(1/γ)²/V_a', 'n^(5/3)']
model_r = [r_m1, r_m2, r_m3, r_m4, r_m5]
bars = ax4.barh(model_names, model_r)
for i, (bar, r) in enumerate(zip(bars, model_r)):
    bar.set_color('green' if r > 0.85 else 'blue' if r > 0.7 else 'gray')
ax4.axvline(x=0.8, color='red', linestyle='--', alpha=0.5, label='r=0.8')
ax4.set_xlabel('Correlation r')
ax4.set_title('Model Comparison: B ∝ ?')
ax4.legend()
ax4.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/bulk_modulus_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: bulk_modulus_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r4) > 0.9:
    print("\n✓ EXCELLENT VALIDATION")
    print(f"  B vs E_coh/V_a: r = {r4:.3f}")
    print("  Bulk modulus = cohesive energy density")
    print("  Coherence enters via E_coh-θ_D relationship")
elif abs(r1) > 0.7 or abs(r2) > 0.7:
    print("\n✓ GOOD VALIDATION")
    print(f"  B vs 1/γ_phonon: r = {r1:.3f}")
    print(f"  B vs θ_D²: r = {r2:.3f}")
    print("  Bulk modulus shows coherence dependence")
else:
    print("\n○ PARTIAL")
    print("  Multiple factors contribute to B")
    print(f"  Best model: B ∝ {best_model[1]}: r = {best_model[0]:.3f}")
