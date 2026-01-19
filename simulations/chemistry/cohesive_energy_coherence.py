#!/usr/bin/env python3
"""
Chemistry Session #116: Cohesive Energy and Coherence

Test whether cohesive energy E_coh relates to coherence parameters.

Hypothesis: E_coh should correlate with 1/γ_phonon since both relate
to bond stiffness:
- E_coh = energy to break all bonds
- θ_D = phonon cutoff from bond stiffness
- γ_phonon = 2T/θ_D

Expected: E_coh ∝ 1/γ_phonon (stronger bonds = more coherent)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Cohesive energy data (eV/atom)
# Sources: Kittel, CRC Handbook, various DFT compilations
materials = {
    # Refractory metals (high E_coh)
    'W':  {'E_coh': 8.90, 'theta_D': 400, 'V_a': 15.8, 'T_m': 3695},
    'Mo': {'E_coh': 6.82, 'theta_D': 450, 'V_a': 15.6, 'T_m': 2896},
    'Ta': {'E_coh': 8.10, 'theta_D': 240, 'V_a': 18.0, 'T_m': 3290},
    'Nb': {'E_coh': 7.57, 'theta_D': 275, 'V_a': 18.0, 'T_m': 2750},

    # 3d transition metals
    'Fe': {'E_coh': 4.28, 'theta_D': 470, 'V_a': 11.8, 'T_m': 1811},
    'Co': {'E_coh': 4.39, 'theta_D': 445, 'V_a': 11.1, 'T_m': 1768},
    'Ni': {'E_coh': 4.44, 'theta_D': 450, 'V_a': 11.0, 'T_m': 1728},
    'Cr': {'E_coh': 4.10, 'theta_D': 630, 'V_a': 12.0, 'T_m': 2180},
    'V':  {'E_coh': 5.31, 'theta_D': 380, 'V_a': 13.8, 'T_m': 2183},
    'Ti': {'E_coh': 4.85, 'theta_D': 420, 'V_a': 17.7, 'T_m': 1941},
    'Mn': {'E_coh': 2.92, 'theta_D': 410, 'V_a': 12.6, 'T_m': 1519},

    # Noble metals
    'Cu': {'E_coh': 3.49, 'theta_D': 343, 'V_a': 11.8, 'T_m': 1358},
    'Ag': {'E_coh': 2.95, 'theta_D': 225, 'V_a': 17.1, 'T_m': 1235},
    'Au': {'E_coh': 3.81, 'theta_D': 165, 'V_a': 17.0, 'T_m': 1337},

    # Alkali metals
    'Li': {'E_coh': 1.63, 'theta_D': 344, 'V_a': 21.3, 'T_m': 454},
    'Na': {'E_coh': 1.11, 'theta_D': 158, 'V_a': 39.5, 'T_m': 371},
    'K':  {'E_coh': 0.93, 'theta_D': 91, 'V_a': 75.6, 'T_m': 337},
    'Rb': {'E_coh': 0.85, 'theta_D': 56, 'V_a': 93.0, 'T_m': 312},
    'Cs': {'E_coh': 0.80, 'theta_D': 38, 'V_a': 115, 'T_m': 302},

    # Simple metals
    'Al': {'E_coh': 3.39, 'theta_D': 428, 'V_a': 16.6, 'T_m': 933},
    'Pb': {'E_coh': 2.03, 'theta_D': 105, 'V_a': 30.3, 'T_m': 601},
    'Sn': {'E_coh': 3.14, 'theta_D': 200, 'V_a': 27.0, 'T_m': 505},
    'Zn': {'E_coh': 1.35, 'theta_D': 327, 'V_a': 15.2, 'T_m': 693},
    'Mg': {'E_coh': 1.51, 'theta_D': 400, 'V_a': 23.2, 'T_m': 923},

    # Semiconductors/Covalent
    'Si': {'E_coh': 4.63, 'theta_D': 645, 'V_a': 20.0, 'T_m': 1687},
    'Ge': {'E_coh': 3.85, 'theta_D': 374, 'V_a': 22.6, 'T_m': 1211},
    'Diamond': {'E_coh': 7.37, 'theta_D': 2230, 'V_a': 5.7, 'T_m': 4000},  # sublimes
}

# Calculate coherence parameters
T = 300  # Room temperature

for mat in materials:
    data = materials[mat]
    data['gamma_phonon'] = 2 * T / data['theta_D']
    data['inv_gamma'] = 1 / data['gamma_phonon']

# Extract arrays
names = list(materials.keys())
E_coh = np.array([materials[m]['E_coh'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
inv_gamma = np.array([materials[m]['inv_gamma'] for m in names])
V_a = np.array([materials[m]['V_a'] for m in names])
T_m = np.array([materials[m]['T_m'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #116: COHESIVE ENERGY AND COHERENCE")
print("=" * 70)

# Test 1: E_coh vs 1/γ_phonon
r1, p1 = stats.pearsonr(E_coh, inv_gamma)
print(f"\n1. E_coh vs 1/γ_phonon: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.8:
    print("   EXCELLENT correlation - E_coh IS coherence-related")
elif abs(r1) > 0.5:
    print("   MODERATE correlation")
else:
    print("   WEAK correlation - may need different model")

# Test 2: E_coh vs θ_D
r2, p2 = stats.pearsonr(E_coh, theta_D)
print(f"\n2. E_coh vs θ_D: r = {r2:.3f}, p = {p2:.2e}")

# Test 3: E_coh vs T_m (should be excellent - both measure bond strength)
r3, p3 = stats.pearsonr(E_coh, T_m)
print(f"\n3. E_coh vs T_m: r = {r3:.3f}, p = {p3:.2e}")
if abs(r3) > 0.9:
    print("   EXCELLENT - validates E_coh as bond strength measure")

# Test 4: E_coh vs V_a (atomic volume)
r4, p4 = stats.pearsonr(E_coh, V_a)
print(f"\n4. E_coh vs V_a: r = {r4:.3f}, p = {p4:.2e}")

# Test 5: E_coh/V_a (cohesive energy density) vs γ_phonon
E_coh_density = E_coh / V_a
r5, p5 = stats.pearsonr(E_coh_density, inv_gamma)
print(f"\n5. E_coh/V_a vs 1/γ_phonon: r = {r5:.3f}, p = {p5:.2e}")
if abs(r5) > abs(r1):
    print("   Energy DENSITY is better coherence correlate")

# Test 6: Log-log relationship
log_E_coh = np.log(E_coh)
log_theta_D = np.log(theta_D)
r6, p6 = stats.pearsonr(log_E_coh, log_theta_D)
print(f"\n6. log(E_coh) vs log(θ_D): r = {r6:.3f}")
if abs(r6) > 0.8:
    slope6, intercept6, _, _, _ = stats.linregress(log_theta_D, log_E_coh)
    print(f"   Power law: E_coh ∝ θ_D^{slope6:.2f}")

# Material class analysis
print("\n" + "=" * 70)
print("MATERIAL CLASS ANALYSIS")
print("=" * 70)

classes = {
    'Refractory': ['W', 'Mo', 'Ta', 'Nb'],
    '3d TM': ['Fe', 'Co', 'Ni', 'Cr', 'V', 'Ti', 'Mn'],
    'Noble': ['Cu', 'Ag', 'Au'],
    'Alkali': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Simple': ['Al', 'Pb', 'Sn', 'Zn', 'Mg'],
    'Covalent': ['Si', 'Ge', 'Diamond'],
}

print(f"\n{'Class':<12} {'Mean E_coh':<12} {'Mean θ_D':<10} {'Mean γ_ph':<10}")
print("-" * 50)

class_stats = {}
for class_name, members in classes.items():
    class_E_coh = [materials[m]['E_coh'] for m in members]
    class_theta = [materials[m]['theta_D'] for m in members]
    class_gamma = [materials[m]['gamma_phonon'] for m in members]
    class_stats[class_name] = {
        'E_coh': np.mean(class_E_coh),
        'theta_D': np.mean(class_theta),
        'gamma_phonon': np.mean(class_gamma)
    }
    print(f"{class_name:<12} {np.mean(class_E_coh):<12.2f} {np.mean(class_theta):<10.0f} {np.mean(class_gamma):<10.2f}")

# Within-class correlations
print("\n" + "-" * 70)
print("WITHIN-CLASS CORRELATIONS (E_coh vs 1/γ_phonon)")
print("-" * 70)

within_class_r = {}
for class_name, members in classes.items():
    if len(members) >= 3:
        class_E = [materials[m]['E_coh'] for m in members]
        class_inv_g = [1/materials[m]['gamma_phonon'] for m in members]
        r_wc, _ = stats.pearsonr(class_E, class_inv_g)
        within_class_r[class_name] = r_wc
        print(f"{class_name}: r = {r_wc:.3f}")

# Model: E_coh ∝ k × r_bond² / V_a (spring constant × displacement squared / volume)
# Since k ∝ θ_D² × M and V_a ∝ r_bond³
# E_coh ∝ θ_D² / V_a^(1/3)

print("\n" + "=" * 70)
print("PHYSICAL MODEL TEST")
print("=" * 70)

# Test: E_coh vs θ_D² / V_a^(1/3)
model_param = theta_D**2 / V_a**(1/3)
r_model, p_model = stats.pearsonr(E_coh, model_param)
print(f"\nE_coh vs θ_D²/V_a^(1/3): r = {r_model:.3f}")
if abs(r_model) > abs(r1):
    print("   Combined model is BETTER than simple 1/γ")

# Test: E_coh vs θ_D² × V_a^(-2/3) (energy per bond ~ force × distance)
model_param2 = theta_D**2 * V_a**(-2/3)
r_model2, p_model2 = stats.pearsonr(E_coh, model_param2)
print(f"E_coh vs θ_D²/V_a^(2/3): r = {r_model2:.3f}")

# Best model: Include coherence factor
# E_coh ∝ (2/γ)^α × θ_D^β
best_alpha = 0.5
best_r = 0
for alpha in np.arange(0, 2.1, 0.1):
    model_test = inv_gamma**alpha * theta_D
    r_test, _ = stats.pearsonr(E_coh, model_test)
    if abs(r_test) > abs(best_r):
        best_r = r_test
        best_alpha = alpha

print(f"\nBest fit: E_coh ∝ (1/γ)^{best_alpha:.1f} × θ_D: r = {best_r:.3f}")

# Print summary
print("\n" + "=" * 70)
print("SESSION #116 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- E_coh vs 1/γ_phonon: r = {r1:.3f}
- E_coh vs θ_D: r = {r2:.3f}
- E_coh vs T_m: r = {r3:.3f}
- E_coh vs V_a: r = {r4:.3f}
- log(E_coh) vs log(θ_D): r = {r6:.3f}

Physical Interpretation:
""")

if abs(r1) > 0.8:
    print("E_coh IS coherence-related:")
    print("  - Stronger bonds → higher θ_D → lower γ_phonon")
    print("  - E_coh ∝ 1/γ_phonon validates cohesion = coherence")
    print("  - Both measure lattice phase-locking strength")
elif abs(r3) > abs(r1):
    print("E_coh relates to T_m more than to coherence:")
    print("  - Both are THERMODYNAMIC bond strength measures")
    print("  - Coherence (γ) modulates but doesn't determine E_coh")
    print("  - E_coh and T_m are DUAL: one sets the other")

print("\nCausal chain:")
print("  Atomic structure → Bond stiffness k → E_coh ~ k×r²")
print("  Same k → θ_D ~ √(k/M) → γ_phonon = 2T/θ_D")
print("  Therefore E_coh and 1/γ_phonon share COMMON CAUSE (k)")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: E_coh vs 1/γ_phonon
ax1 = axes[0, 0]
colors = {'Refractory': 'red', '3d TM': 'blue', 'Noble': 'gold',
          'Alkali': 'green', 'Simple': 'orange', 'Covalent': 'purple'}
for class_name, members in classes.items():
    class_E = [materials[m]['E_coh'] for m in members]
    class_inv_g = [1/materials[m]['gamma_phonon'] for m in members]
    ax1.scatter(class_inv_g, class_E, c=colors[class_name],
                label=class_name, s=100, alpha=0.7)
    for i, m in enumerate(members):
        ax1.annotate(m, (class_inv_g[i], class_E[i]), fontsize=8)

ax1.set_xlabel('1/γ_phonon (coherence)')
ax1.set_ylabel('Cohesive Energy (eV/atom)')
ax1.set_title(f'E_coh vs 1/γ_phonon\nr = {r1:.3f}')
ax1.legend(loc='best', fontsize=8)
ax1.grid(True, alpha=0.3)

# Plot 2: E_coh vs θ_D
ax2 = axes[0, 1]
for class_name, members in classes.items():
    class_E = [materials[m]['E_coh'] for m in members]
    class_theta = [materials[m]['theta_D'] for m in members]
    ax2.scatter(class_theta, class_E, c=colors[class_name],
                label=class_name, s=100, alpha=0.7)

ax2.set_xlabel('Debye Temperature θ_D (K)')
ax2.set_ylabel('Cohesive Energy (eV/atom)')
ax2.set_title(f'E_coh vs θ_D\nr = {r2:.3f}')
ax2.legend(loc='best', fontsize=8)
ax2.grid(True, alpha=0.3)

# Plot 3: E_coh vs T_m
ax3 = axes[1, 0]
for class_name, members in classes.items():
    class_E = [materials[m]['E_coh'] for m in members]
    class_Tm = [materials[m]['T_m'] for m in members]
    ax3.scatter(class_Tm, class_E, c=colors[class_name],
                label=class_name, s=100, alpha=0.7)

ax3.set_xlabel('Melting Point T_m (K)')
ax3.set_ylabel('Cohesive Energy (eV/atom)')
ax3.set_title(f'E_coh vs T_m\nr = {r3:.3f}')
ax3.legend(loc='best', fontsize=8)
ax3.grid(True, alpha=0.3)

# Plot 4: Material class comparison
ax4 = axes[1, 1]
class_names = list(class_stats.keys())
class_E_means = [class_stats[c]['E_coh'] for c in class_names]
class_gamma_means = [class_stats[c]['gamma_phonon'] for c in class_names]

ax4.bar(range(len(class_names)), class_E_means, color=[colors[c] for c in class_names])
ax4.set_xticks(range(len(class_names)))
ax4.set_xticklabels(class_names, rotation=45, ha='right')
ax4.set_ylabel('Mean Cohesive Energy (eV/atom)')
ax4.set_title('E_coh by Material Class')
ax4.grid(True, alpha=0.3, axis='y')

# Add γ_phonon annotations
for i, (e, g) in enumerate(zip(class_E_means, class_gamma_means)):
    ax4.annotate(f'γ={g:.1f}', (i, e + 0.2), ha='center', fontsize=9)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/cohesive_energy_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: cohesive_energy_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r1) > 0.85 and abs(r3) > 0.9:
    print("\n✓ EXCELLENT VALIDATION")
    print(f"  E_coh vs 1/γ_phonon: r = {r1:.3f}")
    print(f"  E_coh vs T_m: r = {r3:.3f}")
    print("  Cohesive energy IS coherence-related via bond stiffness")
elif abs(r1) > 0.7:
    print("\n✓ GOOD VALIDATION")
    print(f"  E_coh vs 1/γ_phonon: r = {r1:.3f}")
    print("  Cohesive energy correlates with coherence")
elif abs(r3) > abs(r1) + 0.2:
    print("\n◐ PARTIAL - THERMODYNAMIC PROPERTY")
    print(f"  E_coh vs T_m: r = {r3:.3f} (stronger)")
    print(f"  E_coh vs 1/γ_phonon: r = {r1:.3f} (weaker)")
    print("  E_coh and γ share common cause but aren't directly related")
else:
    print("\n○ MIXED RESULTS")
    print("  Need further analysis")
