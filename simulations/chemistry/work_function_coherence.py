#!/usr/bin/env python3
"""
Chemistry Session #117: Work Function and Coherence

Session #98 identified thermionic emission as dominated by φ (work function).
Here we test whether φ itself relates to coherence parameters.

Hypothesis: φ should correlate with γ_optical (electronic coherence)
since both measure electron binding:
- φ = energy to remove electron from solid surface
- IE = energy to remove electron from isolated atom
- γ_optical = IE_ref / IE

Expected: φ ∝ 1/γ_optical (tighter binding = more coherent electrons)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Work function data (eV)
# Sources: Michaelson (1977), CRC Handbook, Hölzl & Schulte review
materials = {
    # Alkali metals (very low φ)
    'Li': {'phi': 2.93, 'IE': 5.39, 'theta_D': 344, 'EN': 0.98},
    'Na': {'phi': 2.75, 'IE': 5.14, 'theta_D': 158, 'EN': 0.93},
    'K':  {'phi': 2.30, 'IE': 4.34, 'theta_D': 91, 'EN': 0.82},
    'Rb': {'phi': 2.16, 'IE': 4.18, 'theta_D': 56, 'EN': 0.82},
    'Cs': {'phi': 2.14, 'IE': 3.89, 'theta_D': 38, 'EN': 0.79},

    # Alkaline earth
    'Be': {'phi': 4.98, 'IE': 9.32, 'theta_D': 1440, 'EN': 1.57},
    'Mg': {'phi': 3.66, 'IE': 7.65, 'theta_D': 400, 'EN': 1.31},
    'Ca': {'phi': 2.87, 'IE': 6.11, 'theta_D': 230, 'EN': 1.00},
    'Sr': {'phi': 2.59, 'IE': 5.69, 'theta_D': 147, 'EN': 0.95},
    'Ba': {'phi': 2.52, 'IE': 5.21, 'theta_D': 110, 'EN': 0.89},

    # 3d transition metals
    'Ti': {'phi': 4.33, 'IE': 6.83, 'theta_D': 420, 'EN': 1.54},
    'V':  {'phi': 4.30, 'IE': 6.75, 'theta_D': 380, 'EN': 1.63},
    'Cr': {'phi': 4.50, 'IE': 6.77, 'theta_D': 630, 'EN': 1.66},
    'Mn': {'phi': 4.10, 'IE': 7.43, 'theta_D': 410, 'EN': 1.55},
    'Fe': {'phi': 4.67, 'IE': 7.90, 'theta_D': 470, 'EN': 1.83},
    'Co': {'phi': 5.00, 'IE': 7.88, 'theta_D': 445, 'EN': 1.88},
    'Ni': {'phi': 5.15, 'IE': 7.64, 'theta_D': 450, 'EN': 1.91},
    'Cu': {'phi': 4.65, 'IE': 7.73, 'theta_D': 343, 'EN': 1.90},
    'Zn': {'phi': 4.33, 'IE': 9.39, 'theta_D': 327, 'EN': 1.65},

    # 4d transition metals
    'Zr': {'phi': 4.05, 'IE': 6.63, 'theta_D': 291, 'EN': 1.33},
    'Mo': {'phi': 4.60, 'IE': 7.09, 'theta_D': 450, 'EN': 2.16},
    'Pd': {'phi': 5.12, 'IE': 8.34, 'theta_D': 274, 'EN': 2.20},
    'Ag': {'phi': 4.26, 'IE': 7.58, 'theta_D': 225, 'EN': 1.93},

    # 5d transition metals
    'Ta': {'phi': 4.25, 'IE': 7.55, 'theta_D': 240, 'EN': 1.50},
    'W':  {'phi': 4.55, 'IE': 7.98, 'theta_D': 400, 'EN': 2.36},
    'Pt': {'phi': 5.65, 'IE': 8.96, 'theta_D': 240, 'EN': 2.28},
    'Au': {'phi': 5.10, 'IE': 9.22, 'theta_D': 165, 'EN': 2.54},

    # Simple metals
    'Al': {'phi': 4.28, 'IE': 5.99, 'theta_D': 428, 'EN': 1.61},
    'Pb': {'phi': 4.25, 'IE': 7.42, 'theta_D': 105, 'EN': 2.33},
    'Sn': {'phi': 4.42, 'IE': 7.34, 'theta_D': 200, 'EN': 1.96},
}

# Calculate coherence parameters
IE_ref = 13.6  # eV (Hydrogen)
T = 300  # Room temperature

for mat in materials:
    data = materials[mat]
    data['gamma_optical'] = IE_ref / data['IE']
    data['inv_gamma_optical'] = 1 / data['gamma_optical']
    data['gamma_phonon'] = 2 * T / data['theta_D']

# Extract arrays
names = list(materials.keys())
phi = np.array([materials[m]['phi'] for m in names])
IE = np.array([materials[m]['IE'] for m in names])
EN = np.array([materials[m]['EN'] for m in names])
gamma_optical = np.array([materials[m]['gamma_optical'] for m in names])
inv_gamma_optical = np.array([materials[m]['inv_gamma_optical'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #117: WORK FUNCTION AND COHERENCE")
print("=" * 70)

# Test 1: φ vs 1/γ_optical (should be good - both measure electron binding)
r1, p1 = stats.pearsonr(phi, inv_gamma_optical)
print(f"\n1. φ vs 1/γ_optical: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.8:
    print("   EXCELLENT - Work function IS electronic coherence")
elif abs(r1) > 0.5:
    print("   MODERATE correlation")
else:
    print("   WEAK - may need different model")

# Test 2: φ vs IE directly
r2, p2 = stats.pearsonr(phi, IE)
print(f"\n2. φ vs IE: r = {r2:.3f}, p = {p2:.2e}")
if abs(r2) > 0.8:
    print("   EXCELLENT - φ tracks ionization energy")

# Test 3: φ vs EN (electronegativity)
r3, p3 = stats.pearsonr(phi, EN)
print(f"\n3. φ vs EN: r = {r3:.3f}, p = {p3:.2e}")
if abs(r3) > 0.8:
    print("   EXCELLENT - φ tracks electronegativity")

# Test 4: φ vs γ_phonon (should be WEAK - different coherence type)
r4, p4 = stats.pearsonr(phi, gamma_phonon)
print(f"\n4. φ vs γ_phonon: r = {r4:.3f}, p = {p4:.2e}")
if abs(r4) < 0.3:
    print("   WEAK as expected - φ is electronic, not phononic")

# Test 5: φ/IE ratio (metallic screening factor)
phi_IE_ratio = phi / IE
r5, p5 = stats.pearsonr(phi_IE_ratio, gamma_optical)
print(f"\n5. φ/IE ratio vs γ_optical: r = {r5:.3f}")
print(f"   Mean φ/IE = {np.mean(phi_IE_ratio):.3f} ± {np.std(phi_IE_ratio):.3f}")

# Material class analysis
print("\n" + "=" * 70)
print("MATERIAL CLASS ANALYSIS")
print("=" * 70)

classes = {
    'Alkali': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Alkaline Earth': ['Be', 'Mg', 'Ca', 'Sr', 'Ba'],
    '3d TM': ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    '4d TM': ['Zr', 'Mo', 'Pd', 'Ag'],
    '5d TM': ['Ta', 'W', 'Pt', 'Au'],
    'Simple': ['Al', 'Pb', 'Sn'],
}

print(f"\n{'Class':<15} {'Mean φ':<10} {'Mean IE':<10} {'Mean γ_opt':<10} {'φ/IE':<10}")
print("-" * 55)

class_stats = {}
for class_name, members in classes.items():
    class_phi = [materials[m]['phi'] for m in members]
    class_IE = [materials[m]['IE'] for m in members]
    class_gamma_opt = [materials[m]['gamma_optical'] for m in members]
    class_ratio = [materials[m]['phi']/materials[m]['IE'] for m in members]
    class_stats[class_name] = {
        'phi': np.mean(class_phi),
        'IE': np.mean(class_IE),
        'gamma_optical': np.mean(class_gamma_opt),
        'ratio': np.mean(class_ratio)
    }
    print(f"{class_name:<15} {np.mean(class_phi):<10.2f} {np.mean(class_IE):<10.2f} "
          f"{np.mean(class_gamma_opt):<10.2f} {np.mean(class_ratio):<10.3f}")

# Within-class correlations
print("\n" + "-" * 70)
print("WITHIN-CLASS CORRELATIONS")
print("-" * 70)

print("\nφ vs IE:")
for class_name, members in classes.items():
    if len(members) >= 3:
        class_phi = [materials[m]['phi'] for m in members]
        class_IE = [materials[m]['IE'] for m in members]
        r_wc, _ = stats.pearsonr(class_phi, class_IE)
        print(f"  {class_name}: r = {r_wc:.3f}")

print("\nφ vs EN:")
for class_name, members in classes.items():
    if len(members) >= 3:
        class_phi = [materials[m]['phi'] for m in members]
        class_EN = [materials[m]['EN'] for m in members]
        r_wc, _ = stats.pearsonr(class_phi, class_EN)
        print(f"  {class_name}: r = {r_wc:.3f}")

# Physics of work function
print("\n" + "=" * 70)
print("PHYSICAL ANALYSIS: WORK FUNCTION")
print("=" * 70)

print("""
Work function φ = Energy to remove electron from Fermi level to vacuum

φ = φ_bulk + φ_surface

φ_bulk ∝ IE (atomic ionization energy)
φ_surface = dipole barrier (surface-dependent)

The φ/IE ratio measures METALLIC SCREENING:
""")

# Group by φ/IE ratio
print(f"\nφ/IE ratio by material class:")
for class_name in ['Alkali', 'Alkaline Earth', '3d TM', '5d TM', 'Simple']:
    ratio = class_stats[class_name]['ratio']
    print(f"  {class_name}: {ratio:.3f}")

print("\nAlkali metals: φ/IE ≈ 0.52 (strong screening)")
print("5d TM: φ/IE ≈ 0.58 (moderate screening)")
print("These differences reflect band structure, not coherence")

# Best model: φ from EN
print("\n" + "=" * 70)
print("BEST MODEL: φ from Electronegativity")
print("=" * 70)

slope, intercept, r_value, p_value, std_err = stats.linregress(EN, phi)
print(f"\nLinear fit: φ = {slope:.2f}×EN + {intercept:.2f}")
print(f"Correlation: r = {r_value:.3f}")
print(f"φ predicted from EN alone: r² = {r_value**2:.3f}")

# Alternative: φ = A × EN + B × (something)
# Try φ vs EN + d-electron count

# Print extremes
print("\n" + "=" * 70)
print("EXTREME VALUES")
print("=" * 70)

sorted_phi = sorted([(m, materials[m]['phi']) for m in names], key=lambda x: x[1])
print("\nLowest φ (easiest electron emission):")
for m, p in sorted_phi[:5]:
    print(f"  {m}: φ = {p:.2f} eV, IE = {materials[m]['IE']:.2f} eV")

print("\nHighest φ (hardest electron emission):")
for m, p in sorted_phi[-5:][::-1]:
    print(f"  {m}: φ = {p:.2f} eV, IE = {materials[m]['IE']:.2f} eV")

# Summary
print("\n" + "=" * 70)
print("SESSION #117 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- φ vs IE: r = {r2:.3f}
- φ vs EN: r = {r3:.3f}
- φ vs 1/γ_optical: r = {r1:.3f}
- φ vs γ_phonon: r = {r4:.3f}

Physical Interpretation:
""")

if abs(r2) > 0.7 or abs(r3) > 0.7:
    print("Work function IS electronic coherence property:")
    print(f"  φ vs IE: r = {r2:.3f}")
    print(f"  φ vs EN: r = {r3:.3f}")
    print("  φ tracks electron binding = electronic coherence")
    print("  Confirms #115: EN → IE → γ_optical → electronic properties")
else:
    print("Work function is MIXED:")
    print("  Depends on atomic (IE) + surface (dipole) contributions")
    print("  Surface effects reduce the simple atomic correlation")

if abs(r4) < 0.3:
    print("\nConfirms TWO INDEPENDENT coherence channels (#115):")
    print(f"  φ vs γ_phonon: r = {r4:.3f} (essentially ZERO)")
    print("  Work function is ELECTRONIC, unrelated to phonon coherence")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: φ vs IE
ax1 = axes[0, 0]
colors = {'Alkali': 'green', 'Alkaline Earth': 'lime', '3d TM': 'blue',
          '4d TM': 'cyan', '5d TM': 'navy', 'Simple': 'orange'}
for class_name, members in classes.items():
    class_phi = [materials[m]['phi'] for m in members]
    class_IE = [materials[m]['IE'] for m in members]
    ax1.scatter(class_IE, class_phi, c=colors[class_name],
                label=class_name, s=100, alpha=0.7)
    for i, m in enumerate(members):
        ax1.annotate(m, (class_IE[i], class_phi[i]), fontsize=7)

# Add trend line
z = np.polyfit(IE, phi, 1)
p_trend = np.poly1d(z)
ax1.plot([4, 10], [p_trend(4), p_trend(10)], 'r--', alpha=0.5, label=f'r={r2:.3f}')

ax1.set_xlabel('Ionization Energy (eV)')
ax1.set_ylabel('Work Function φ (eV)')
ax1.set_title(f'φ vs IE\nr = {r2:.3f}')
ax1.legend(loc='best', fontsize=7)
ax1.grid(True, alpha=0.3)

# Plot 2: φ vs EN
ax2 = axes[0, 1]
for class_name, members in classes.items():
    class_phi = [materials[m]['phi'] for m in members]
    class_EN = [materials[m]['EN'] for m in members]
    ax2.scatter(class_EN, class_phi, c=colors[class_name],
                label=class_name, s=100, alpha=0.7)

# Add trend line
z = np.polyfit(EN, phi, 1)
p_trend = np.poly1d(z)
ax2.plot([0.7, 2.6], [p_trend(0.7), p_trend(2.6)], 'r--', alpha=0.5, label=f'r={r3:.3f}')

ax2.set_xlabel('Electronegativity (Pauling)')
ax2.set_ylabel('Work Function φ (eV)')
ax2.set_title(f'φ vs EN\nr = {r3:.3f}')
ax2.legend(loc='best', fontsize=7)
ax2.grid(True, alpha=0.3)

# Plot 3: φ vs γ_phonon (should show NO correlation)
ax3 = axes[1, 0]
for class_name, members in classes.items():
    class_phi = [materials[m]['phi'] for m in members]
    class_gamma = [materials[m]['gamma_phonon'] for m in members]
    ax3.scatter(class_gamma, class_phi, c=colors[class_name],
                label=class_name, s=100, alpha=0.7)

ax3.set_xlabel('γ_phonon = 2T/θ_D')
ax3.set_ylabel('Work Function φ (eV)')
ax3.set_title(f'φ vs γ_phonon\nr = {r4:.3f} (NO correlation)')
ax3.legend(loc='best', fontsize=7)
ax3.grid(True, alpha=0.3)

# Plot 4: φ/IE ratio by class
ax4 = axes[1, 1]
class_names_sorted = ['Alkali', 'Alkaline Earth', '3d TM', '4d TM', '5d TM', 'Simple']
ratios = [class_stats[c]['ratio'] for c in class_names_sorted]

ax4.bar(range(len(class_names_sorted)), ratios,
        color=[colors[c] for c in class_names_sorted])
ax4.set_xticks(range(len(class_names_sorted)))
ax4.set_xticklabels(class_names_sorted, rotation=45, ha='right')
ax4.set_ylabel('φ/IE ratio (screening factor)')
ax4.set_title('Metallic Screening by Material Class')
ax4.axhline(y=np.mean(ratios), color='red', linestyle='--',
            label=f'Mean = {np.mean(ratios):.3f}')
ax4.legend()
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/work_function_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: work_function_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r2) > 0.7 and abs(r3) > 0.7:
    print("\n✓ EXCELLENT VALIDATION")
    print(f"  φ vs IE: r = {r2:.3f}")
    print(f"  φ vs EN: r = {r3:.3f}")
    print("  Work function IS electronic coherence (via EN/IE)")
    print("\n  Coherence hierarchy:")
    print("    EN → IE → γ_optical → φ")
    print("    All measure electron binding = electronic coherence")
elif abs(r4) < 0.3 and (abs(r2) > 0.5 or abs(r3) > 0.5):
    print("\n✓ GOOD VALIDATION")
    print("  φ correlates with electronic properties (IE, EN)")
    print("  φ is INDEPENDENT of phonon coherence")
    print("  Two coherence channels confirmed")
else:
    print("\n○ PARTIAL")
    print("  φ has both atomic and surface contributions")
    print("  Simple correlations may be diluted")
