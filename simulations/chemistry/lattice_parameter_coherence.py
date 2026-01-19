#!/usr/bin/env python3
"""
Chemistry Session #121: Lattice Parameter and Coherence

Test lattice parameter a₀ as coherence determinant.
Related to V_a (#114) and r_a (#119) but may show different behavior.

Hypothesis: a₀ → bond length → spring constant → θ_D → γ_phonon
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Lattice parameter data (Å)
# For FCC: V_a = a³/4, BCC: V_a = a³/2, SC: V_a = a³
# Sources: Wyckoff, Kittel, various crystallographic databases
materials = {
    # FCC metals
    'Cu':  {'a': 3.61, 'structure': 'FCC', 'theta_D': 343, 'mass': 63.5, 'Z': 4},
    'Ag':  {'a': 4.09, 'structure': 'FCC', 'theta_D': 225, 'mass': 107.9, 'Z': 4},
    'Au':  {'a': 4.08, 'structure': 'FCC', 'theta_D': 165, 'mass': 197.0, 'Z': 4},
    'Al':  {'a': 4.05, 'structure': 'FCC', 'theta_D': 428, 'mass': 27.0, 'Z': 4},
    'Ni':  {'a': 3.52, 'structure': 'FCC', 'theta_D': 450, 'mass': 58.7, 'Z': 4},
    'Pt':  {'a': 3.92, 'structure': 'FCC', 'theta_D': 240, 'mass': 195.1, 'Z': 4},
    'Pb':  {'a': 4.95, 'structure': 'FCC', 'theta_D': 105, 'mass': 207.2, 'Z': 4},
    'Ca':  {'a': 5.59, 'structure': 'FCC', 'theta_D': 230, 'mass': 40.1, 'Z': 4},
    'Sr':  {'a': 6.08, 'structure': 'FCC', 'theta_D': 147, 'mass': 87.6, 'Z': 4},

    # BCC metals
    'Li':  {'a': 3.49, 'structure': 'BCC', 'theta_D': 344, 'mass': 6.94, 'Z': 2},
    'Na':  {'a': 4.29, 'structure': 'BCC', 'theta_D': 158, 'mass': 23.0, 'Z': 2},
    'K':   {'a': 5.33, 'structure': 'BCC', 'theta_D': 91, 'mass': 39.1, 'Z': 2},
    'Rb':  {'a': 5.71, 'structure': 'BCC', 'theta_D': 56, 'mass': 85.5, 'Z': 2},
    'Cs':  {'a': 6.14, 'structure': 'BCC', 'theta_D': 38, 'mass': 132.9, 'Z': 2},
    'Fe':  {'a': 2.87, 'structure': 'BCC', 'theta_D': 470, 'mass': 55.8, 'Z': 2},
    'Cr':  {'a': 2.88, 'structure': 'BCC', 'theta_D': 630, 'mass': 52.0, 'Z': 2},
    'V':   {'a': 3.02, 'structure': 'BCC', 'theta_D': 380, 'mass': 50.9, 'Z': 2},
    'Mo':  {'a': 3.15, 'structure': 'BCC', 'theta_D': 450, 'mass': 95.9, 'Z': 2},
    'W':   {'a': 3.16, 'structure': 'BCC', 'theta_D': 400, 'mass': 183.8, 'Z': 2},
    'Ta':  {'a': 3.30, 'structure': 'BCC', 'theta_D': 240, 'mass': 180.9, 'Z': 2},
    'Nb':  {'a': 3.30, 'structure': 'BCC', 'theta_D': 275, 'mass': 92.9, 'Z': 2},

    # HCP metals
    'Mg':  {'a': 3.21, 'structure': 'HCP', 'theta_D': 400, 'mass': 24.3, 'Z': 2, 'c': 5.21},
    'Zn':  {'a': 2.66, 'structure': 'HCP', 'theta_D': 327, 'mass': 65.4, 'Z': 2, 'c': 4.95},
    'Ti':  {'a': 2.95, 'structure': 'HCP', 'theta_D': 420, 'mass': 47.9, 'Z': 2, 'c': 4.69},
    'Co':  {'a': 2.51, 'structure': 'HCP', 'theta_D': 445, 'mass': 58.9, 'Z': 2, 'c': 4.07},
    'Be':  {'a': 2.29, 'structure': 'HCP', 'theta_D': 1440, 'mass': 9.01, 'Z': 2, 'c': 3.58},

    # Diamond structure
    'C':   {'a': 3.57, 'structure': 'Diamond', 'theta_D': 2230, 'mass': 12.0, 'Z': 8},
    'Si':  {'a': 5.43, 'structure': 'Diamond', 'theta_D': 645, 'mass': 28.1, 'Z': 8},
    'Ge':  {'a': 5.66, 'structure': 'Diamond', 'theta_D': 374, 'mass': 72.6, 'Z': 8},
    'Sn':  {'a': 6.49, 'structure': 'Diamond', 'theta_D': 200, 'mass': 118.7, 'Z': 8},  # gray Sn
}

# Calculate derived quantities
T = 300

for mat in materials:
    data = materials[mat]
    data['gamma_phonon'] = 2 * T / data['theta_D']

    # Calculate atomic volume from lattice parameter
    if data['structure'] == 'FCC':
        data['V_a'] = data['a']**3 / 4
    elif data['structure'] == 'BCC':
        data['V_a'] = data['a']**3 / 2
    elif data['structure'] == 'HCP':
        # V = √3/2 × a² × c, with 2 atoms per unit cell
        data['V_a'] = (np.sqrt(3)/2 * data['a']**2 * data['c']) / 2
    elif data['structure'] == 'Diamond':
        data['V_a'] = data['a']**3 / 8

    # Nearest neighbor distance
    if data['structure'] == 'FCC':
        data['d_nn'] = data['a'] / np.sqrt(2)
    elif data['structure'] == 'BCC':
        data['d_nn'] = data['a'] * np.sqrt(3) / 2
    elif data['structure'] == 'HCP':
        data['d_nn'] = data['a']  # Simplified
    elif data['structure'] == 'Diamond':
        data['d_nn'] = data['a'] * np.sqrt(3) / 4

# Extract arrays
names = list(materials.keys())
a = np.array([materials[m]['a'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in names])
V_a = np.array([materials[m]['V_a'] for m in names])
mass = np.array([materials[m]['mass'] for m in names])
d_nn = np.array([materials[m]['d_nn'] for m in names])

print("=" * 70)
print("CHEMISTRY SESSION #121: LATTICE PARAMETER AND COHERENCE")
print("=" * 70)

print(f"\nDataset: {len(names)} materials")
print(f"a₀ range: {np.min(a):.2f} - {np.max(a):.2f} Å")

# Test 1: a vs γ_phonon
r1, p1 = stats.pearsonr(a, gamma_phonon)
print(f"\n1. a₀ vs γ_phonon: r = {r1:.3f}, p = {p1:.2e}")

# Test 2: d_nn vs γ_phonon (nearest neighbor distance)
r2, p2 = stats.pearsonr(d_nn, gamma_phonon)
print(f"2. d_nn vs γ_phonon: r = {r2:.3f}")

# Test 3: V_a vs γ_phonon (confirms #114)
r3, p3 = stats.pearsonr(V_a, gamma_phonon)
print(f"3. V_a vs γ_phonon: r = {r3:.3f}")
if abs(r3) > 0.9:
    print("   Confirms Session #114 (V_a EXCELLENT)")

# Test 4: a vs θ_D
r4, p4 = stats.pearsonr(a, theta_D)
print(f"\n4. a₀ vs θ_D: r = {r4:.3f}")

# Test 5: d_nn/√M vs θ_D (Debye model expects correlation)
d_over_sqrtM = d_nn / np.sqrt(mass)
r5, p5 = stats.pearsonr(d_over_sqrtM, theta_D)
print(f"5. d_nn/√M vs θ_D: r = {r5:.3f}")

# Test 6: Power law
log_a = np.log(a)
log_theta = np.log(theta_D)
r6, _ = stats.pearsonr(log_a, log_theta)
slope6, intercept6, _, _, _ = stats.linregress(log_a, log_theta)
print(f"\n6. log(a) vs log(θ_D): r = {r6:.3f}")
print(f"   Power law: θ_D ∝ a^{slope6:.2f}")

# Structure analysis
print("\n" + "=" * 70)
print("STRUCTURE-DEPENDENT ANALYSIS")
print("=" * 70)

structures = ['FCC', 'BCC', 'HCP', 'Diamond']
print(f"\n{'Structure':<10} {'Mean a':<10} {'Mean θ_D':<10} {'Mean γ_ph':<10} {'Count'}")
print("-" * 50)

for struct in structures:
    members = [m for m in names if materials[m]['structure'] == struct]
    if len(members) >= 2:
        struct_a = [materials[m]['a'] for m in members]
        struct_theta = [materials[m]['theta_D'] for m in members]
        struct_gamma = [materials[m]['gamma_phonon'] for m in members]
        print(f"{struct:<10} {np.mean(struct_a):<10.2f} {np.mean(struct_theta):<10.0f} "
              f"{np.mean(struct_gamma):<10.2f} {len(members)}")

# Within-structure correlations
print("\n" + "-" * 70)
print("WITHIN-STRUCTURE CORRELATIONS (a vs γ_phonon)")
print("-" * 70)

for struct in structures:
    members = [m for m in names if materials[m]['structure'] == struct]
    if len(members) >= 4:
        struct_a = [materials[m]['a'] for m in members]
        struct_gamma = [materials[m]['gamma_phonon'] for m in members]
        r_ws, _ = stats.pearsonr(struct_a, struct_gamma)
        print(f"  {struct}: r = {r_ws:.3f} (n={len(members)})")

# Nearest neighbor vs lattice parameter
print("\n" + "=" * 70)
print("NEAREST NEIGHBOR ANALYSIS")
print("=" * 70)

print(f"\nNearest neighbor distance d_nn:")
print(f"{'Element':<5} {'a (Å)':<8} {'d_nn (Å)':<10} {'θ_D (K)':<8} {'γ_phonon':<10}")
print("-" * 50)

# Sort by d_nn
sorted_mats = sorted(names, key=lambda m: materials[m]['d_nn'])
for m in sorted_mats[:10]:
    print(f"{m:<5} {materials[m]['a']:<8.2f} {materials[m]['d_nn']:<10.2f} "
          f"{materials[m]['theta_D']:<8} {materials[m]['gamma_phonon']:<10.2f}")

print("\n...")
for m in sorted_mats[-5:]:
    print(f"{m:<5} {materials[m]['a']:<8.2f} {materials[m]['d_nn']:<10.2f} "
          f"{materials[m]['theta_D']:<8} {materials[m]['gamma_phonon']:<10.2f}")

# Summary
print("\n" + "=" * 70)
print("SESSION #121 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- a₀ vs γ_phonon: r = {r1:.3f}
- d_nn vs γ_phonon: r = {r2:.3f}
- V_a vs γ_phonon: r = {r3:.3f}
- θ_D ∝ a^{slope6:.2f}: r = {r6:.3f}

Comparison:
- V_a (#114): r = 0.956
- r_cov (#119): r = 0.796
- a₀ (#121): r = {r1:.3f}
- d_nn (#121): r = {r2:.3f}
""")

if abs(r2) > abs(r1):
    print("Nearest neighbor distance d_nn is BETTER than lattice parameter a₀")
    print("Bond length directly sets phonon scale")
else:
    print("Lattice parameter shows similar correlation to other length scales")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: a vs γ_phonon by structure
ax1 = axes[0, 0]
colors = {'FCC': 'blue', 'BCC': 'red', 'HCP': 'green', 'Diamond': 'purple'}

for struct in structures:
    members = [m for m in names if materials[m]['structure'] == struct]
    struct_a = [materials[m]['a'] for m in members]
    struct_gamma = [materials[m]['gamma_phonon'] for m in members]
    ax1.scatter(struct_a, struct_gamma, c=colors[struct], label=struct, s=80, alpha=0.7)
    for i, m in enumerate(members):
        ax1.annotate(m, (struct_a[i], struct_gamma[i]), fontsize=7)

ax1.set_xlabel('Lattice Parameter a₀ (Å)')
ax1.set_ylabel('γ_phonon = 2T/θ_D')
ax1.set_title(f'a₀ vs γ_phonon\nr = {r1:.3f}')
ax1.legend(loc='best')
ax1.grid(True, alpha=0.3)

# Plot 2: d_nn vs γ_phonon
ax2 = axes[0, 1]
for struct in structures:
    members = [m for m in names if materials[m]['structure'] == struct]
    struct_d = [materials[m]['d_nn'] for m in members]
    struct_gamma = [materials[m]['gamma_phonon'] for m in members]
    ax2.scatter(struct_d, struct_gamma, c=colors[struct], label=struct, s=80, alpha=0.7)

ax2.set_xlabel('Nearest Neighbor Distance d_nn (Å)')
ax2.set_ylabel('γ_phonon = 2T/θ_D')
ax2.set_title(f'd_nn vs γ_phonon\nr = {r2:.3f}')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

# Plot 3: V_a vs γ_phonon (confirms #114)
ax3 = axes[1, 0]
for struct in structures:
    members = [m for m in names if materials[m]['structure'] == struct]
    struct_V = [materials[m]['V_a'] for m in members]
    struct_gamma = [materials[m]['gamma_phonon'] for m in members]
    ax3.scatter(struct_V, struct_gamma, c=colors[struct], label=struct, s=80, alpha=0.7)

ax3.set_xlabel('Atomic Volume V_a (Å³)')
ax3.set_ylabel('γ_phonon = 2T/θ_D')
ax3.set_title(f'V_a vs γ_phonon\nr = {r3:.3f}')
ax3.legend(loc='best')
ax3.grid(True, alpha=0.3)

# Plot 4: log-log power law
ax4 = axes[1, 1]
for struct in structures:
    members = [m for m in names if materials[m]['structure'] == struct]
    struct_log_a = [np.log10(materials[m]['a']) for m in members]
    struct_log_theta = [np.log10(materials[m]['theta_D']) for m in members]
    ax4.scatter(struct_log_a, struct_log_theta, c=colors[struct], label=struct, s=80, alpha=0.7)
    for i, m in enumerate(members):
        ax4.annotate(m, (struct_log_a[i], struct_log_theta[i]), fontsize=6, alpha=0.8)

# Add power law fit
x_fit = np.linspace(0.3, 0.85, 50)
y_fit = slope6 * x_fit + intercept6/np.log(10)  # Convert intercept
ax4.plot(x_fit, y_fit, 'k--', alpha=0.5, label=f'θ_D ∝ a^{slope6:.1f}')

ax4.set_xlabel('log₁₀(a₀)')
ax4.set_ylabel('log₁₀(θ_D)')
ax4.set_title(f'Power Law: θ_D ∝ a^{slope6:.2f}\nr = {r6:.3f}')
ax4.legend(loc='best')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/lattice_parameter_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: lattice_parameter_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r3) > 0.9:
    print("\n✓ EXCELLENT VALIDATION (via V_a)")
    print(f"  V_a vs γ_phonon: r = {r3:.3f}")
    print("  Confirms Session #114")
    print("\n  Hierarchy of length scales:")
    print(f"    V_a: r = {r3:.3f} (best)")
    print(f"    d_nn: r = {r2:.3f}")
    print(f"    a₀: r = {r1:.3f}")
elif abs(r2) > 0.8:
    print("\n✓ GOOD VALIDATION")
    print(f"  d_nn vs γ_phonon: r = {r2:.3f}")
    print("  Nearest neighbor distance sets phonon scale")
else:
    print("\n○ MODERATE")
    print("  All length scales show similar correlations")
    print("  V_a (#114) remains best coherence determinant")
