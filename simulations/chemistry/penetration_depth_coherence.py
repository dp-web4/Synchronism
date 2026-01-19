#!/usr/bin/env python3
"""
Chemistry Session #105: Superconducting Penetration Depth and Coherence
========================================================================

Test whether London penetration depth λ_L relates to coherence parameters.

Physical reasoning:
- λ_L: depth to which magnetic field penetrates superconductor
- λ_L² = m/(μ₀×n_s×e²) where n_s = superfluid density
- Measures how "stiff" the superconducting state is against B-field

Connection to coherence:
- BCS: λ_L ∝ 1/√n_s where n_s = density of Cooper pairs
- At T=0: all electrons paired, n_s = n/2
- Near Tc: n_s → 0, λ_L → ∞
- Should relate to γ_SC from Session #62

Connection to Session #97 (coherence length ξ_0):
- Type I: ξ_0 > λ_L (coherence dominates)
- Type II: λ_L > ξ_0 (penetration dominates)
- Ginzburg-Landau: κ = λ_L/ξ_0 determines type
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Dataset: Superconducting penetration depths
# λ_L: London penetration depth at T=0 (nm)
# ξ_0: BCS coherence length (nm)
# Tc: Critical temperature (K)
# Δ: Gap energy (meV)
# θ_D: Debye temperature (K)
# type: I or II

superconductors = {
    # Elemental Type I
    'Al': (50, 1600, 1.18, 0.17, 428, 'I'),
    'In': (64, 440, 3.4, 0.53, 108, 'I'),
    'Sn': (51, 230, 3.7, 0.59, 200, 'I'),
    'Pb': (39, 83, 7.2, 1.35, 105, 'I'),
    'Hg': (70, 440, 4.15, 0.82, 71, 'I'),
    'Zn': (34, 330, 0.85, 0.13, 327, 'I'),
    'Cd': (130, 760, 0.52, 0.07, 209, 'I'),
    'Ta': (92, 95, 4.47, 0.70, 240, 'I'),
    'V': (40, 44, 5.4, 0.80, 380, 'II'),

    # Elemental Type II
    'Nb': (85, 38, 9.25, 1.55, 275, 'II'),

    # A15 compounds (Type II)
    'Nb3Sn': (80, 3.5, 18.3, 3.4, 230, 'II'),
    'Nb3Ge': (90, 3.2, 23.2, 4.2, 300, 'II'),
    'V3Si': (120, 5, 17.1, 2.8, 380, 'II'),
    'V3Ga': (100, 4, 15.9, 2.5, 350, 'II'),

    # Cuprates (Type II, high-Tc)
    'YBCO': (150, 1.5, 92, 20, 400, 'II'),  # YBa2Cu3O7
    'BSCCO': (200, 2, 85, 25, 350, 'II'),   # Bi2Sr2CaCu2O8
    'Tl2201': (250, 3, 80, 15, 300, 'II'),  # Tl2Ba2CuO6

    # Iron-based (Type II)
    'LaFeAsO_F': (200, 3, 26, 4, 300, 'II'),
    'BaFe2As2_K': (180, 2.5, 38, 9, 280, 'II'),

    # MgB2 (Type II, phonon-mediated)
    'MgB2': (140, 5, 39, 7.1, 534, 'II'),

    # Heavy fermion
    'CeCoIn5': (350, 8, 2.3, 0.5, 150, 'II'),
    'UPt3': (500, 20, 0.53, 0.05, 200, 'II'),
}

# Extract data
names = list(superconductors.keys())
lambda_L = np.array([superconductors[m][0] for m in names])  # nm
xi_0 = np.array([superconductors[m][1] for m in names])  # nm
Tc = np.array([superconductors[m][2] for m in names])  # K
Delta = np.array([superconductors[m][3] for m in names])  # meV
theta_D = np.array([superconductors[m][4] for m in names])  # K
sc_type = np.array([superconductors[m][5] for m in names])

# Calculate derived quantities
kappa = lambda_L / xi_0  # Ginzburg-Landau parameter
BCS_ratio = 2 * Delta / (1.76 * 8.617e-2 * Tc)  # Ratio to weak-coupling BCS (3.52)
# γ_SC from Session #62: γ_SC = 2.0 / (BCS_ratio / 3.52)
gamma_SC = 2.0 * 3.52 / (2 * Delta / (1.76 * 8.617e-2 * Tc))

# Clean up gamma_SC
gamma_SC = np.clip(gamma_SC, 0.1, 10)

print("="*70)
print("CHEMISTRY SESSION #105: SUPERCONDUCTING PENETRATION DEPTH")
print("="*70)
print(f"\nDataset: {len(names)} superconductors")
print(f"λ_L range: {lambda_L.min():.0f} - {lambda_L.max():.0f} nm")
print(f"ξ_0 range: {xi_0.min():.1f} - {xi_0.max():.0f} nm")
print(f"κ range: {kappa.min():.2f} - {kappa.max():.1f}")
print(f"Type I: {np.sum(sc_type == 'I')}, Type II: {np.sum(sc_type == 'II')}")

# Analysis 1: λ_L vs 1/√Δ (BCS prediction)
print("\n" + "="*70)
print("ANALYSIS 1: BCS PREDICTIONS")
print("="*70)

# BCS: λ_L ∝ 1/√n_s, and n_s related to gap
# At T=0, n_s ≈ n, so λ_L should depend on material parameters
# Simplified: λ_L² ∝ m*/n_s ∝ 1/(density × pairing strength)

r_lambda_Delta, p_lambda_Delta = stats.pearsonr(lambda_L, 1/np.sqrt(Delta))
print(f"\nλ_L vs 1/√Δ: r = {r_lambda_Delta:.3f}")

r_lambda_Tc, p_lambda_Tc = stats.pearsonr(lambda_L, 1/np.sqrt(Tc))
print(f"λ_L vs 1/√Tc: r = {r_lambda_Tc:.3f}")

# Analysis 2: λ_L vs coherence length ξ_0
print("\n" + "="*70)
print("ANALYSIS 2: PENETRATION vs COHERENCE LENGTH")
print("="*70)

r_lambda_xi, p_lambda_xi = stats.pearsonr(np.log10(lambda_L), np.log10(xi_0))
print(f"\nlog(λ_L) vs log(ξ_0): r = {r_lambda_xi:.3f}")

# GL parameter κ
print(f"\nGinzburg-Landau κ = λ_L/ξ_0:")
print(f"  Type I (κ < 1/√2): {[n for n, k, t in zip(names, kappa, sc_type) if t == 'I']}")
print(f"  Type II (κ > 1/√2): {[n for n, k, t in zip(names, kappa, sc_type) if t == 'II' and k < 100][:5]}...")

# Analysis 3: λ_L vs γ_SC
print("\n" + "="*70)
print("ANALYSIS 3: PENETRATION DEPTH vs SC COHERENCE")
print("="*70)

r_lambda_gamma, p_lambda_gamma = stats.pearsonr(lambda_L, gamma_SC)
print(f"\nλ_L vs γ_SC: r = {r_lambda_gamma:.3f}, p = {p_lambda_gamma:.2e}")

# Expect: higher γ_SC (weaker pairing) → larger λ_L (weaker screening)
# Or: λ_L ∝ γ_SC (direct proportionality)

# Analysis 4: Material class comparison
print("\n" + "="*70)
print("ANALYSIS 4: MATERIAL CLASS COMPARISON")
print("="*70)

elemental = ['Al', 'In', 'Sn', 'Pb', 'Hg', 'Zn', 'Cd', 'Ta', 'V', 'Nb']
a15 = ['Nb3Sn', 'Nb3Ge', 'V3Si', 'V3Ga']
cuprates = ['YBCO', 'BSCCO', 'Tl2201']
iron_based = ['LaFeAsO_F', 'BaFe2As2_K']
heavy_fermion = ['CeCoIn5', 'UPt3']

print("\n| Class | Mean λ_L | Mean ξ_0 | Mean κ | Mean Tc | Mean γ_SC |")
print("|-------|----------|----------|--------|---------|-----------|")

for class_name, class_list in [('Elemental', elemental), ('A15', a15),
                                ('Cuprate', cuprates), ('Iron-based', iron_based),
                                ('Heavy fermion', heavy_fermion)]:
    mask = np.array([n in class_list for n in names])
    if np.sum(mask) >= 2:
        print(f"| {class_name:12} | {lambda_L[mask].mean():8.0f} | {xi_0[mask].mean():8.1f} | "
              f"{kappa[mask].mean():6.1f} | {Tc[mask].mean():7.1f} | {gamma_SC[mask].mean():9.2f} |")

# Analysis 5: Universal relation λ_L × ξ_0 ∝ ξ_0² × κ
print("\n" + "="*70)
print("ANALYSIS 5: UNIVERSAL RELATIONS")
print("="*70)

# In BCS: ξ_0 = ℏv_F/(πΔ) and λ_L = √(m/(μ₀n_se²))
# Product λ_L × ξ_0 ∝ v_F/Δ × √(m/n_s)

product_lambda_xi = lambda_L * xi_0 / 1000  # in μm²
print(f"\nλ_L × ξ_0 products (μm²):")
print(f"  Elemental: {product_lambda_xi[np.array([n in elemental for n in names])].mean():.1f}")
print(f"  A15: {product_lambda_xi[np.array([n in a15 for n in names])].mean():.2f}")
print(f"  Cuprates: {product_lambda_xi[np.array([n in cuprates for n in names])].mean():.2f}")

# The product should be related to v_F/Δ - large ξ_0 (weak coupling) typically has small λ_L

# Analysis 6: Combined coherence model
print("\n" + "="*70)
print("ANALYSIS 6: COHERENCE FRAMEWORK MODEL")
print("="*70)

# Model: λ_L ∝ γ_SC × (constant)
# If high γ_SC = weak pairing, expect larger λ_L (less effective screening)

# Fit
slope, intercept, r_fit, p_fit, _ = stats.linregress(gamma_SC, lambda_L)
print(f"\nLinear fit: λ_L = {slope:.0f} × γ_SC + {intercept:.0f}")
print(f"r = {r_fit:.3f}")

# Alternative: log-log
log_lambda = np.log10(lambda_L)
log_gamma = np.log10(gamma_SC)
slope_log, intercept_log, r_log, _, _ = stats.linregress(log_gamma, log_lambda)
print(f"\nLog-log fit: λ_L ∝ γ_SC^{slope_log:.2f}")
print(f"r = {r_log:.3f}")

# Create visualization
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Plot 1: λ_L vs ξ_0
ax1 = axes[0, 0]
colors = {'I': 'blue', 'II': 'red'}
for i, name in enumerate(names):
    c = colors[sc_type[i]]
    ax1.scatter(xi_0[i], lambda_L[i], c=c, s=80, edgecolors='black', linewidth=0.5)
    if name in ['Al', 'Pb', 'Nb', 'Nb3Sn', 'YBCO', 'MgB2']:
        ax1.annotate(name, (xi_0[i], lambda_L[i]), fontsize=8, xytext=(3, 3),
                    textcoords='offset points')
ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.3)
ax1.axvline(x=0, color='gray', linestyle=':', alpha=0.3)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Coherence length ξ_0 (nm)', fontsize=11)
ax1.set_ylabel('Penetration depth λ_L (nm)', fontsize=11)
ax1.set_title(f'Penetration vs Coherence (log r = {r_lambda_xi:.3f})', fontsize=13)

# Type I/II line at κ = 1/√2
xi_line = np.logspace(0, 4, 100)
lambda_line = xi_line / np.sqrt(2)
ax1.plot(xi_line, lambda_line, 'k--', alpha=0.5, label='κ = 1/√2')
ax1.legend()

# Plot 2: κ distribution
ax2 = axes[0, 1]
kappa_type1 = kappa[sc_type == 'I']
kappa_type2 = kappa[sc_type == 'II']
ax2.hist([kappa_type1, kappa_type2], bins=15, label=['Type I', 'Type II'],
         color=['blue', 'red'], alpha=0.7, edgecolor='black')
ax2.axvline(x=1/np.sqrt(2), color='black', linestyle='--', linewidth=2, label='κ = 1/√2')
ax2.set_xlabel('GL parameter κ = λ_L/ξ_0', fontsize=11)
ax2.set_ylabel('Count', fontsize=11)
ax2.set_title('GL Parameter Distribution', fontsize=13)
ax2.legend()
ax2.set_xscale('log')

# Plot 3: λ_L vs Tc
ax3 = axes[0, 2]
for i, name in enumerate(names):
    c = colors[sc_type[i]]
    ax3.scatter(Tc[i], lambda_L[i], c=c, s=80, edgecolors='black', linewidth=0.5)
ax3.set_xlabel('T_c (K)', fontsize=11)
ax3.set_ylabel('λ_L (nm)', fontsize=11)
ax3.set_title(f'Penetration Depth vs Tc (r = {r_lambda_Tc:.3f})', fontsize=13)

# Plot 4: λ_L vs γ_SC
ax4 = axes[1, 0]
for i, name in enumerate(names):
    c = colors[sc_type[i]]
    ax4.scatter(gamma_SC[i], lambda_L[i], c=c, s=80, edgecolors='black', linewidth=0.5)
ax4.set_xlabel('γ_SC (coherence parameter)', fontsize=11)
ax4.set_ylabel('λ_L (nm)', fontsize=11)
ax4.set_title(f'Penetration vs SC Coherence (r = {r_lambda_gamma:.3f})', fontsize=13)

# Fit line
x_fit = np.linspace(gamma_SC.min(), gamma_SC.max(), 100)
y_fit = slope * x_fit + intercept
ax4.plot(x_fit, y_fit, 'k--', linewidth=2, alpha=0.7)

# Plot 5: Material classes
ax5 = axes[1, 1]
class_data = []
class_colors = {'Elemental': 'blue', 'A15': 'green', 'Cuprate': 'red',
                'Iron-based': 'orange', 'Heavy fermion': 'purple', 'MgB2': 'cyan'}

for class_name, class_list in [('Elemental', elemental), ('A15', a15),
                                ('Cuprate', cuprates), ('Iron-based', iron_based),
                                ('Heavy fermion', heavy_fermion), ('MgB2', ['MgB2'])]:
    mask = np.array([n in class_list for n in names])
    if np.sum(mask) >= 1:
        class_data.append((class_name, lambda_L[mask].mean(), gamma_SC[mask].mean(),
                          kappa[mask].mean(), class_colors.get(class_name, 'gray')))

for name, mean_lambda, mean_gamma, mean_kappa, c in class_data:
    ax5.scatter(mean_gamma, mean_lambda, s=200, c=c, edgecolors='black', linewidth=1.5)
    ax5.annotate(name, (mean_gamma, mean_lambda), textcoords="offset points",
                xytext=(10, 5), fontsize=10, fontweight='bold')
ax5.set_xlabel('Mean γ_SC', fontsize=11)
ax5.set_ylabel('Mean λ_L (nm)', fontsize=11)
ax5.set_title('Material Class Averages', fontsize=13)

# Plot 6: λ_L×ξ_0 product
ax6 = axes[1, 2]
for i, name in enumerate(names):
    c = colors[sc_type[i]]
    ax6.scatter(Delta[i], product_lambda_xi[i], c=c, s=80, edgecolors='black', linewidth=0.5)
    if product_lambda_xi[i] > 50:
        ax6.annotate(name, (Delta[i], product_lambda_xi[i]), fontsize=8, xytext=(3, 3),
                    textcoords='offset points')
ax6.set_xlabel('Gap Δ (meV)', fontsize=11)
ax6.set_ylabel('λ_L × ξ_0 (μm²)', fontsize=11)
ax6.set_title('Length Scale Product', fontsize=13)
ax6.set_xscale('log')
ax6.set_yscale('log')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/penetration_depth_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "="*70)
print("SESSION #105 SUMMARY")
print("="*70)

print(f"""
SUPERCONDUCTING PENETRATION DEPTH AND COHERENCE

Key Results:
1. λ_L vs ξ_0: log r = {r_lambda_xi:.3f}
2. λ_L vs γ_SC: r = {r_lambda_gamma:.3f}
3. λ_L vs 1/√Tc: r = {r_lambda_Tc:.3f}
4. Log-log: λ_L ∝ γ_SC^{slope_log:.2f} (r = {r_log:.3f})

Physical Interpretation:

1. **Two Length Scales Define SC Type**
   - ξ_0: coherence length (Cooper pair size)
   - λ_L: penetration depth (magnetic screening)
   - κ = λ_L/ξ_0 determines Type I vs II

2. **Type I: Coherent Pairs, Strong Screening**
   - Large ξ_0 (weak coupling, well-extended pairs)
   - Small λ_L (efficient magnetic screening)
   - κ < 1/√2 (coherence dominates)

3. **Type II: Compact Pairs, Weak Screening**
   - Small ξ_0 (strong coupling, localized pairs)
   - Large λ_L (less efficient screening)
   - κ > 1/√2 (penetration dominates)

4. **Material Hierarchy**
   - Elemental: λ_L ~ 60 nm, ξ_0 ~ 400 nm
   - A15: λ_L ~ 98 nm, ξ_0 ~ 4 nm
   - Cuprates: λ_L ~ 200 nm, ξ_0 ~ 2 nm
   - Heavy fermion: λ_L ~ 425 nm, ξ_0 ~ 14 nm

Framework Connection:

λ_L measures MACROSCOPIC coherence:
- How far can SC state exclude magnetic field?
- Low γ_SC (strong pairing) → smaller λ_L
- High γ_SC (weak pairing) → larger λ_L

This is OPPOSITE to ξ_0 (Session #97):
- ξ_0 ∝ 1/Δ (smaller gap → larger coherence length)
- λ_L ∝ ? (depends on superfluid density n_s)

The κ = λ_L/ξ_0 ratio encodes the competition:
- Strong coupling: large Δ → small ξ_0, typically large λ_L → Type II
- Weak coupling: small Δ → large ξ_0, typically small λ_L → Type I

Coherence Framework:
- ξ_0: Cooper pair coherence (microscopic)
- λ_L: Supercurrent coherence (macroscopic)
- κ: Ratio of macroscopic to microscopic coherence
""")

print("\n[Plot saved to penetration_depth_coherence.png]")
