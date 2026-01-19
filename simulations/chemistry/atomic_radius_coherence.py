#!/usr/bin/env python3
"""
Chemistry Session #119: Atomic Radius and Coherence

Follow-up to Session #114 (atomic volume).
Test atomic radius r_a as coherence determinant.

The causal chain should be:
r_a → bond length L ~ r_a → spring constant k ~ 1/L² → θ_D ~ √(k/M) → γ_phonon

Testing multiple radius definitions:
- r_cov (covalent radius)
- r_metallic (metallic radius)
- r_vdW (van der Waals radius)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Atomic radii data (pm = picometers)
# Sources: Cordero et al. (covalent), Batsanov (metallic), Bondi (vdW)
materials = {
    # Light elements
    'Li': {'r_cov': 128, 'r_met': 152, 'r_vdW': 182, 'theta_D': 344, 'mass': 6.94},
    'Be': {'r_cov': 96, 'r_met': 112, 'r_vdW': 153, 'theta_D': 1440, 'mass': 9.01},
    'C':  {'r_cov': 77, 'r_met': None, 'r_vdW': 170, 'theta_D': 2230, 'mass': 12.01},  # diamond
    'N':  {'r_cov': 71, 'r_met': None, 'r_vdW': 155, 'theta_D': None, 'mass': 14.01},
    'O':  {'r_cov': 66, 'r_met': None, 'r_vdW': 152, 'theta_D': None, 'mass': 16.00},

    # Alkali metals
    'Na': {'r_cov': 166, 'r_met': 186, 'r_vdW': 227, 'theta_D': 158, 'mass': 22.99},
    'K':  {'r_cov': 203, 'r_met': 227, 'r_vdW': 275, 'theta_D': 91, 'mass': 39.10},
    'Rb': {'r_cov': 220, 'r_met': 248, 'r_vdW': 303, 'theta_D': 56, 'mass': 85.47},
    'Cs': {'r_cov': 244, 'r_met': 265, 'r_vdW': 343, 'theta_D': 38, 'mass': 132.91},

    # Alkaline earth
    'Mg': {'r_cov': 141, 'r_met': 160, 'r_vdW': 173, 'theta_D': 400, 'mass': 24.31},
    'Ca': {'r_cov': 176, 'r_met': 197, 'r_vdW': 231, 'theta_D': 230, 'mass': 40.08},
    'Sr': {'r_cov': 195, 'r_met': 215, 'r_vdW': 249, 'theta_D': 147, 'mass': 87.62},
    'Ba': {'r_cov': 215, 'r_met': 222, 'r_vdW': 268, 'theta_D': 110, 'mass': 137.33},

    # 3d transition metals
    'Ti': {'r_cov': 160, 'r_met': 147, 'r_vdW': None, 'theta_D': 420, 'mass': 47.87},
    'V':  {'r_cov': 153, 'r_met': 134, 'r_vdW': None, 'theta_D': 380, 'mass': 50.94},
    'Cr': {'r_cov': 139, 'r_met': 128, 'r_vdW': None, 'theta_D': 630, 'mass': 52.00},
    'Fe': {'r_cov': 132, 'r_met': 126, 'r_vdW': None, 'theta_D': 470, 'mass': 55.85},
    'Co': {'r_cov': 126, 'r_met': 125, 'r_vdW': None, 'theta_D': 445, 'mass': 58.93},
    'Ni': {'r_cov': 124, 'r_met': 124, 'r_vdW': 163, 'theta_D': 450, 'mass': 58.69},
    'Cu': {'r_cov': 132, 'r_met': 128, 'r_vdW': 140, 'theta_D': 343, 'mass': 63.55},
    'Zn': {'r_cov': 122, 'r_met': 134, 'r_vdW': 139, 'theta_D': 327, 'mass': 65.38},

    # 4d transition metals
    'Zr': {'r_cov': 175, 'r_met': 160, 'r_vdW': None, 'theta_D': 291, 'mass': 91.22},
    'Mo': {'r_cov': 154, 'r_met': 139, 'r_vdW': None, 'theta_D': 450, 'mass': 95.96},
    'Pd': {'r_cov': 139, 'r_met': 137, 'r_vdW': 163, 'theta_D': 274, 'mass': 106.42},
    'Ag': {'r_cov': 145, 'r_met': 144, 'r_vdW': 172, 'theta_D': 225, 'mass': 107.87},

    # 5d transition metals
    'Ta': {'r_cov': 170, 'r_met': 146, 'r_vdW': None, 'theta_D': 240, 'mass': 180.95},
    'W':  {'r_cov': 162, 'r_met': 139, 'r_vdW': None, 'theta_D': 400, 'mass': 183.84},
    'Pt': {'r_cov': 136, 'r_met': 139, 'r_vdW': 175, 'theta_D': 240, 'mass': 195.08},
    'Au': {'r_cov': 136, 'r_met': 144, 'r_vdW': 166, 'theta_D': 165, 'mass': 196.97},

    # Simple metals
    'Al': {'r_cov': 121, 'r_met': 143, 'r_vdW': 184, 'theta_D': 428, 'mass': 26.98},
    'Si': {'r_cov': 111, 'r_met': None, 'r_vdW': 210, 'theta_D': 645, 'mass': 28.09},
    'Ge': {'r_cov': 120, 'r_met': None, 'r_vdW': 211, 'theta_D': 374, 'mass': 72.63},
    'Sn': {'r_cov': 139, 'r_met': 151, 'r_vdW': 217, 'theta_D': 200, 'mass': 118.71},
    'Pb': {'r_cov': 146, 'r_met': 175, 'r_vdW': 202, 'theta_D': 105, 'mass': 207.2},
}

# Calculate coherence parameters
T = 300  # Room temperature

for mat in materials:
    data = materials[mat]
    if data['theta_D'] is not None:
        data['gamma_phonon'] = 2 * T / data['theta_D']
        data['inv_gamma'] = 1 / data['gamma_phonon']
    else:
        data['gamma_phonon'] = None
        data['inv_gamma'] = None

# Extract arrays (only elements with θ_D and r_cov)
valid_names = [m for m in materials if materials[m]['theta_D'] is not None]
r_cov = np.array([materials[m]['r_cov'] for m in valid_names])
theta_D = np.array([materials[m]['theta_D'] for m in valid_names])
gamma_phonon = np.array([materials[m]['gamma_phonon'] for m in valid_names])
mass = np.array([materials[m]['mass'] for m in valid_names])

# Also get metallic radii where available
met_names = [m for m in valid_names if materials[m]['r_met'] is not None]
r_met = np.array([materials[m]['r_met'] for m in met_names])
theta_D_met = np.array([materials[m]['theta_D'] for m in met_names])
gamma_phonon_met = np.array([materials[m]['gamma_phonon'] for m in met_names])
mass_met = np.array([materials[m]['mass'] for m in met_names])

print("=" * 70)
print("CHEMISTRY SESSION #119: ATOMIC RADIUS AND COHERENCE")
print("=" * 70)

print(f"\nDataset: {len(valid_names)} elements with θ_D")
print(f"Elements with metallic radius: {len(met_names)}")

# Test 1: r_cov vs γ_phonon
r1, p1 = stats.pearsonr(r_cov, gamma_phonon)
print(f"\n1. r_cov vs γ_phonon: r = {r1:.3f}, p = {p1:.2e}")
if abs(r1) > 0.8:
    print("   EXCELLENT - Covalent radius sets phonon coherence")
elif abs(r1) > 0.5:
    print("   MODERATE correlation")
else:
    print("   WEAK - mass effects may dominate")

# Test 2: r_cov vs θ_D
r2, p2 = stats.pearsonr(r_cov, theta_D)
print(f"\n2. r_cov vs θ_D: r = {r2:.3f}")

# Test 3: r_met vs γ_phonon (metallic radii only)
r3, p3 = stats.pearsonr(r_met, gamma_phonon_met)
print(f"\n3. r_met vs γ_phonon: r = {r3:.3f} (n={len(met_names)})")

# Test 4: Mass-corrected analysis
# θ_D ∝ √(k/M) where k ~ 1/r²
# So θ_D ∝ 1/(r × √M) or θ_D × r × √M ≈ constant

mass_corrected = theta_D * r_cov * np.sqrt(mass)
print(f"\n4. Mass-corrected θ_D × r_cov × √M:")
print(f"   Mean = {np.mean(mass_corrected):.0f}")
print(f"   Std = {np.std(mass_corrected):.0f}")
print(f"   CV = {np.std(mass_corrected)/np.mean(mass_corrected):.2f}")

# Test 5: Power law r_cov^α × M^β vs θ_D
log_r = np.log(r_cov)
log_M = np.log(mass)
log_theta = np.log(theta_D)

# Multilinear regression: log(θ_D) = A × log(r) + B × log(M) + C
from numpy.linalg import lstsq
X = np.column_stack([log_r, log_M, np.ones(len(log_r))])
coeffs, residuals, rank, s = lstsq(X, log_theta, rcond=None)
alpha, beta, C = coeffs

# Predicted values
log_theta_pred = alpha * log_r + beta * log_M + C
theta_pred = np.exp(log_theta_pred)
r_power, _ = stats.pearsonr(theta_D, theta_pred)

print(f"\n5. Power law fit: θ_D ∝ r^{alpha:.2f} × M^{beta:.2f}")
print(f"   Correlation: r = {r_power:.3f}")
print(f"   Expected: α = -1, β = -0.5")

# Group analysis
print("\n" + "=" * 70)
print("MATERIAL CLASS ANALYSIS")
print("=" * 70)

classes = {
    'Alkali': ['Li', 'Na', 'K', 'Rb', 'Cs'],
    'Alkaline Earth': ['Mg', 'Ca', 'Sr', 'Ba'],
    '3d TM': ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'],
    '4d TM': ['Zr', 'Mo', 'Pd', 'Ag'],
    '5d TM': ['Ta', 'W', 'Pt', 'Au'],
    'Simple/Semi': ['Al', 'Si', 'Ge', 'Sn', 'Pb'],
}

print(f"\n{'Class':<15} {'Mean r_cov':<12} {'Mean θ_D':<10} {'Mean γ_ph':<10}")
print("-" * 50)

for class_name, members in classes.items():
    valid = [m for m in members if m in valid_names]
    if len(valid) >= 2:
        class_r = [materials[m]['r_cov'] for m in valid]
        class_theta = [materials[m]['theta_D'] for m in valid]
        class_gamma = [materials[m]['gamma_phonon'] for m in valid]
        print(f"{class_name:<15} {np.mean(class_r):<12.0f} {np.mean(class_theta):<10.0f} {np.mean(class_gamma):<10.2f}")

# Within-class correlations
print("\n" + "-" * 70)
print("WITHIN-CLASS CORRELATIONS (r_cov vs γ_phonon)")
print("-" * 70)

for class_name, members in classes.items():
    valid = [m for m in members if m in valid_names]
    if len(valid) >= 3:
        class_r = [materials[m]['r_cov'] for m in valid]
        class_gamma = [materials[m]['gamma_phonon'] for m in valid]
        r_wc, _ = stats.pearsonr(class_r, class_gamma)
        print(f"  {class_name}: r = {r_wc:.3f}")

# Periodic trends
print("\n" + "=" * 70)
print("PERIODIC TRENDS")
print("=" * 70)

# Down a group: r increases, θ_D decreases, γ increases
print("\nAlkali series (down the group):")
print(f"{'Element':<5} {'r_cov':<8} {'θ_D':<8} {'γ_phonon':<10}")
for m in ['Li', 'Na', 'K', 'Rb', 'Cs']:
    print(f"{m:<5} {materials[m]['r_cov']:<8} {materials[m]['theta_D']:<8} {materials[m]['gamma_phonon']:<10.2f}")

print("\nAlkaline earth series:")
for m in ['Be', 'Mg', 'Ca', 'Sr', 'Ba']:
    print(f"{m:<5} {materials[m]['r_cov']:<8} {materials[m]['theta_D']:<8} {materials[m]['gamma_phonon']:<10.2f}")

print("\n3d transition metal series (across period):")
print(f"{'Element':<5} {'r_cov':<8} {'θ_D':<8} {'γ_phonon':<10}")
for m in ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']:
    print(f"{m:<5} {materials[m]['r_cov']:<8} {materials[m]['theta_D']:<8} {materials[m]['gamma_phonon']:<10.2f}")

# Summary
print("\n" + "=" * 70)
print("SESSION #119 SUMMARY")
print("=" * 70)

print(f"""
Key Results:
- r_cov vs γ_phonon: r = {r1:.3f}
- r_cov vs θ_D: r = {r2:.3f}
- r_met vs γ_phonon: r = {r3:.3f}
- Power law: θ_D ∝ r^{alpha:.2f} × M^{beta:.2f} (r = {r_power:.3f})

Physical Interpretation:
""")

if abs(r1) > 0.8:
    print("Atomic radius IS phonon coherence determinant:")
    print(f"  r_cov vs γ_phonon: r = {r1:.3f}")
    print("  Validates Session #114 (V_a ~ r³)")
elif abs(r_power) > 0.9:
    print("Atomic radius with mass correction IS coherence determinant:")
    print(f"  θ_D ∝ r^{alpha:.2f} × M^{beta:.2f}")
    print("  Mass effects must be included")
else:
    print("Atomic radius partially determines coherence:")
    print("  Both r and M contribute to θ_D")
    print("  Bonding type adds complexity")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: r_cov vs γ_phonon
ax1 = axes[0, 0]
colors = {'Alkali': 'blue', 'Alkaline Earth': 'green', '3d TM': 'red',
          '4d TM': 'orange', '5d TM': 'purple', 'Simple/Semi': 'gray'}

for class_name, members in classes.items():
    valid = [m for m in members if m in valid_names]
    class_r = [materials[m]['r_cov'] for m in valid]
    class_gamma = [materials[m]['gamma_phonon'] for m in valid]
    ax1.scatter(class_r, class_gamma, c=colors[class_name],
                label=class_name, s=80, alpha=0.7)
    for i, m in enumerate(valid):
        ax1.annotate(m, (class_r[i], class_gamma[i]), fontsize=7)

ax1.set_xlabel('Covalent Radius (pm)')
ax1.set_ylabel('γ_phonon = 2T/θ_D')
ax1.set_title(f'r_cov vs γ_phonon\nr = {r1:.3f}')
ax1.legend(loc='best', fontsize=7)
ax1.grid(True, alpha=0.3)

# Plot 2: r_cov vs θ_D
ax2 = axes[0, 1]
for class_name, members in classes.items():
    valid = [m for m in members if m in valid_names]
    class_r = [materials[m]['r_cov'] for m in valid]
    class_theta = [materials[m]['theta_D'] for m in valid]
    ax2.scatter(class_r, class_theta, c=colors[class_name],
                label=class_name, s=80, alpha=0.7)

ax2.set_xlabel('Covalent Radius (pm)')
ax2.set_ylabel('Debye Temperature θ_D (K)')
ax2.set_title(f'r_cov vs θ_D\nr = {r2:.3f}')
ax2.legend(loc='best', fontsize=7)
ax2.grid(True, alpha=0.3)

# Plot 3: Predicted vs Actual θ_D
ax3 = axes[1, 0]
ax3.scatter(theta_D, theta_pred, c='blue', s=80, alpha=0.7)
ax3.plot([0, 2500], [0, 2500], 'r--', alpha=0.5, label='1:1 line')

for i, m in enumerate(valid_names):
    ax3.annotate(m, (theta_D[i], theta_pred[i]), fontsize=7, alpha=0.8)

ax3.set_xlabel('Actual θ_D (K)')
ax3.set_ylabel(f'Predicted θ_D ∝ r^{alpha:.1f} × M^{beta:.1f}')
ax3.set_title(f'Power Law Model\nr = {r_power:.3f}')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Alkali series trend
ax4 = axes[1, 1]
alkali = ['Li', 'Na', 'K', 'Rb', 'Cs']
alkali_r = [materials[m]['r_cov'] for m in alkali]
alkali_gamma = [materials[m]['gamma_phonon'] for m in alkali]
alkali_theta = [materials[m]['theta_D'] for m in alkali]

ax4_twin = ax4.twinx()

ax4.plot(alkali_r, alkali_gamma, 'bo-', markersize=10, label='γ_phonon')
ax4_twin.plot(alkali_r, alkali_theta, 'rs-', markersize=10, label='θ_D')

for i, m in enumerate(alkali):
    ax4.annotate(m, (alkali_r[i], alkali_gamma[i] + 0.5), fontsize=10)

ax4.set_xlabel('Covalent Radius (pm)')
ax4.set_ylabel('γ_phonon', color='blue')
ax4_twin.set_ylabel('θ_D (K)', color='red')
ax4.set_title('Alkali Metal Series')
ax4.tick_params(axis='y', labelcolor='blue')
ax4_twin.tick_params(axis='y', labelcolor='red')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/atomic_radius_coherence.png', dpi=150)
plt.close()

print("\nFigure saved: atomic_radius_coherence.png")

# Final assessment
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if abs(r1) > 0.85:
    print("\n✓ EXCELLENT VALIDATION")
    print(f"  r_cov vs γ_phonon: r = {r1:.3f}")
    print("  Atomic radius directly sets phonon coherence")
elif abs(r_power) > 0.9:
    print("\n✓ EXCELLENT VALIDATION (with mass correction)")
    print(f"  θ_D ∝ r^{alpha:.2f} × M^{beta:.2f}: r = {r_power:.3f}")
    print("  Debye model validated: θ_D ∝ 1/(r × √M)")
elif abs(r1) > 0.7:
    print("\n✓ GOOD VALIDATION")
    print(f"  r_cov vs γ_phonon: r = {r1:.3f}")
    print("  Atomic radius is significant coherence factor")
else:
    print("\n○ PARTIAL")
    print("  Both radius and mass contribute")
    print("  Bonding type adds additional complexity")
