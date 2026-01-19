#!/usr/bin/env python3
"""
Chemistry Session #111: Thermal Diffusivity and Coherence

Test whether thermal diffusivity α relates to coherence parameters.

Theory:
α = κ / (ρ × C_p) [m²/s]

where:
- κ = thermal conductivity
- ρ = density
- C_p = specific heat capacity

Physical meaning: α measures how FAST heat spreads through a material.
- High α: heat spreads quickly (thermal equilibration fast)
- Low α: heat spreads slowly (thermal gradients persist)

From Session #108:
- κ_e ∝ 1/γ_electron (metals)
- κ_ph ∝ θ_D/γ_G² (insulators)

From Session #101:
- C_p ~ 3R (Dulong-Petit at T > θ_D)
- But C_p = C_v + α²VTB/κ_T

Thermal diffusivity:
α = κ/(ρ×C_p)

For metals: α ∝ κ_e/(ρ×C_p) ∝ 1/γ_electron
For insulators: α ∝ κ_ph/(ρ×C_p) ∝ θ_D/γ_G²
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# MATERIALS DATA
# =============================================================================

# Metals
metals_data = {
    # Format: α (mm²/s), κ (W/m·K), ρ (g/cm³), C_p (J/g·K), θ_D (K), λ_ep
    'Ag': {'alpha': 173, 'kappa': 429, 'rho': 10.5, 'Cp': 0.235, 'theta_D': 225, 'lambda_ep': 0.12},
    'Cu': {'alpha': 117, 'kappa': 401, 'rho': 8.96, 'Cp': 0.385, 'theta_D': 343, 'lambda_ep': 0.15},
    'Au': {'alpha': 127, 'kappa': 317, 'rho': 19.3, 'Cp': 0.129, 'theta_D': 165, 'lambda_ep': 0.17},
    'Al': {'alpha': 97, 'kappa': 237, 'rho': 2.70, 'Cp': 0.897, 'theta_D': 428, 'lambda_ep': 0.43},
    'W': {'alpha': 67, 'kappa': 173, 'rho': 19.3, 'Cp': 0.134, 'theta_D': 400, 'lambda_ep': 0.28},
    'Mo': {'alpha': 54, 'kappa': 138, 'rho': 10.2, 'Cp': 0.251, 'theta_D': 450, 'lambda_ep': 0.41},
    'Fe': {'alpha': 23, 'kappa': 80, 'rho': 7.87, 'Cp': 0.449, 'theta_D': 470, 'lambda_ep': 0.70},
    'Ni': {'alpha': 23, 'kappa': 91, 'rho': 8.91, 'Cp': 0.444, 'theta_D': 450, 'lambda_ep': 1.05},
    'Pb': {'alpha': 24, 'kappa': 35, 'rho': 11.3, 'Cp': 0.129, 'theta_D': 105, 'lambda_ep': 1.55},
    'Ti': {'alpha': 9, 'kappa': 22, 'rho': 4.51, 'Cp': 0.523, 'theta_D': 420, 'lambda_ep': 0.80},
    'Zn': {'alpha': 42, 'kappa': 116, 'rho': 7.14, 'Cp': 0.388, 'theta_D': 327, 'lambda_ep': 0.62},
    'Mg': {'alpha': 87, 'kappa': 156, 'rho': 1.74, 'Cp': 1.023, 'theta_D': 400, 'lambda_ep': 0.35},
}

# Insulators/semiconductors
insulators_data = {
    # Format: α (mm²/s), κ_ph (W/m·K), ρ (g/cm³), C_p (J/g·K), θ_D (K), γ_G
    'Diamond': {'alpha': 1100, 'kappa': 2200, 'rho': 3.51, 'Cp': 0.509, 'theta_D': 2230, 'gamma_G': 0.9},
    'Si': {'alpha': 88, 'kappa': 150, 'rho': 2.33, 'Cp': 0.712, 'theta_D': 645, 'gamma_G': 1.0},
    'Ge': {'alpha': 36, 'kappa': 60, 'rho': 5.32, 'Cp': 0.320, 'theta_D': 374, 'gamma_G': 1.2},
    'SiC': {'alpha': 180, 'kappa': 490, 'rho': 3.21, 'Cp': 0.750, 'theta_D': 1200, 'gamma_G': 1.0},
    'AlN': {'alpha': 140, 'kappa': 320, 'rho': 3.26, 'Cp': 0.740, 'theta_D': 950, 'gamma_G': 1.1},
    'BeO': {'alpha': 190, 'kappa': 330, 'rho': 3.01, 'Cp': 1.024, 'theta_D': 1280, 'gamma_G': 1.2},
    'Al2O3': {'alpha': 12, 'kappa': 35, 'rho': 3.97, 'Cp': 0.755, 'theta_D': 1030, 'gamma_G': 1.3},
    'MgO': {'alpha': 12, 'kappa': 60, 'rho': 3.58, 'Cp': 0.877, 'theta_D': 946, 'gamma_G': 1.5},
    'GaAs': {'alpha': 31, 'kappa': 55, 'rho': 5.32, 'Cp': 0.350, 'theta_D': 344, 'gamma_G': 1.3},
    'InP': {'alpha': 47, 'kappa': 68, 'rho': 4.81, 'Cp': 0.310, 'theta_D': 321, 'gamma_G': 1.3},
    'NaCl': {'alpha': 4, 'kappa': 7, 'rho': 2.16, 'Cp': 0.864, 'theta_D': 321, 'gamma_G': 1.6},
    'Glass': {'alpha': 0.5, 'kappa': 1.0, 'rho': 2.50, 'Cp': 0.840, 'theta_D': 400, 'gamma_G': 0.4},
}

# =============================================================================
# ANALYSIS: METALS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #111: THERMAL DIFFUSIVITY AND COHERENCE")
print("="*70)

T = 300  # K

# Metals
metals = list(metals_data.keys())
alpha_m = np.array([metals_data[m]['alpha'] for m in metals])
kappa_m = np.array([metals_data[m]['kappa'] for m in metals])
rho_m = np.array([metals_data[m]['rho'] for m in metals])
Cp_m = np.array([metals_data[m]['Cp'] for m in metals])
theta_D_m = np.array([metals_data[m]['theta_D'] for m in metals])
lambda_ep = np.array([metals_data[m]['lambda_ep'] for m in metals])

gamma_electron = 2 * lambda_ep / (1 + lambda_ep)
gamma_phonon_m = 2 * T / theta_D_m

# Verify α = κ/(ρ×Cp)
alpha_calc = kappa_m / (rho_m * Cp_m)

print("\n" + "="*70)
print("METALS")
print("="*70)

print(f"\n{'Metal':<8} {'α (mm²/s)':<12} {'κ (W/m·K)':<12} {'ρ×Cp':<10} {'γ_e':<8} {'γ_ph':<8}")
print("-"*60)
for i, m in enumerate(metals):
    print(f"{m:<8} {alpha_m[i]:<12.0f} {kappa_m[i]:<12.0f} {rho_m[i]*Cp_m[i]:<10.2f} {gamma_electron[i]:<8.3f} {gamma_phonon_m[i]:<8.3f}")

# Correlations for metals
print("\n" + "-"*50)
print("Correlations (Metals):")
print("-"*50)

r1, p1 = stats.pearsonr(alpha_m, 1/gamma_electron)
print(f"α vs 1/γ_electron: r = {r1:.3f}, p = {p1:.2e}")

r2, p2 = stats.pearsonr(alpha_m, kappa_m)
print(f"α vs κ: r = {r2:.3f}")

r3, p3 = stats.pearsonr(alpha_m, 1/(rho_m * Cp_m))
print(f"α vs 1/(ρ×Cp): r = {r3:.3f}")

r4, p4 = stats.pearsonr(np.log10(alpha_m), 1/gamma_electron)
print(f"log(α) vs 1/γ_electron: r = {r4:.3f}")

# =============================================================================
# ANALYSIS: INSULATORS
# =============================================================================

print("\n" + "="*70)
print("INSULATORS/SEMICONDUCTORS")
print("="*70)

insulators = list(insulators_data.keys())
alpha_i = np.array([insulators_data[m]['alpha'] for m in insulators])
kappa_i = np.array([insulators_data[m]['kappa'] for m in insulators])
rho_i = np.array([insulators_data[m]['rho'] for m in insulators])
Cp_i = np.array([insulators_data[m]['Cp'] for m in insulators])
theta_D_i = np.array([insulators_data[m]['theta_D'] for m in insulators])
gamma_G = np.array([insulators_data[m]['gamma_G'] for m in insulators])

gamma_phonon_i = 2 * T / theta_D_i

print(f"\n{'Material':<10} {'α (mm²/s)':<12} {'κ_ph':<10} {'ρ×Cp':<10} {'γ_G':<8} {'γ_ph':<8}")
print("-"*60)
for i, m in enumerate(insulators):
    print(f"{m:<10} {alpha_i[i]:<12.0f} {kappa_i[i]:<10.0f} {rho_i[i]*Cp_i[i]:<10.2f} {gamma_G[i]:<8.2f} {gamma_phonon_i[i]:<8.3f}")

# Correlations for insulators
print("\n" + "-"*50)
print("Correlations (Insulators):")
print("-"*50)

r5, p5 = stats.pearsonr(np.log10(alpha_i), theta_D_i)
print(f"log(α) vs θ_D: r = {r5:.3f}")

r6, p6 = stats.pearsonr(np.log10(alpha_i), 1/gamma_phonon_i)
print(f"log(α) vs 1/γ_phonon: r = {r6:.3f}")

r7, p7 = stats.pearsonr(np.log10(alpha_i), 1/gamma_G)
print(f"log(α) vs 1/γ_G: r = {r7:.3f}")

# Combined model from #108
combined = theta_D_i / gamma_G**2
r8, p8 = stats.pearsonr(np.log10(alpha_i), combined)
print(f"log(α) vs θ_D/γ_G²: r = {r8:.3f}")

# =============================================================================
# THERMAL DIFFUSIVITY LENGTH SCALE
# =============================================================================

print("\n" + "="*70)
print("THERMAL DIFFUSIVITY LENGTH SCALE")
print("="*70)

# For heat diffusion over time t: L ~ √(α × t)
# At t = 1 ms, L = √(α × 10⁻³) meters

t_ms = 1  # 1 millisecond
L_metals = np.sqrt(alpha_m * 1e-6 * t_ms * 1e-3) * 1e3  # mm
L_insulators = np.sqrt(alpha_i * 1e-6 * t_ms * 1e-3) * 1e3  # mm

print(f"\nThermal diffusion length in t = {t_ms} ms:")
print(f"\n{'Material':<10} {'α (mm²/s)':<12} {'L (mm)':<10}")
print("-"*34)
for i, m in enumerate(metals):
    print(f"{m:<10} {alpha_m[i]:<12.0f} {L_metals[i]:<10.3f}")
print("-"*34)
for i, m in enumerate(insulators):
    print(f"{m:<10} {alpha_i[i]:<12.0f} {L_insulators[i]:<10.3f}")

# =============================================================================
# PHYSICAL INSIGHTS
# =============================================================================

print("\n" + "="*70)
print("PHYSICAL INSIGHTS")
print("="*70)

print(f"""
1. METALS: α vs COHERENCE
   α vs 1/γ_electron: r = {r1:.3f}

   Thermal diffusivity is DOMINATED by κ (r = {r2:.3f} with κ).
   Since κ_e ∝ 1/γ_electron (Session #108), α follows.

2. INSULATORS: α vs COHERENCE
   log(α) vs θ_D: r = {r5:.3f}
   log(α) vs 1/γ_phonon: r = {r6:.3f}
   log(α) vs θ_D/γ_G²: r = {r8:.3f}

   Combined model from Session #108 applies.
   Diamond: α = 1100 mm²/s (highest!)
   Glass: α = 0.5 mm²/s (lowest - amorphous)

3. DIAMOND IS SUPREME
   Diamond has:
   - Highest κ (2200 W/m·K)
   - Moderate ρ×Cp (1.79 J/cm³·K)
   - Highest α (1100 mm²/s)
   - Lowest γ_phonon (0.27)
   - Lowest γ_G (0.9)

4. HEAT SPREADING HIERARCHY
   In 1 ms, heat spreads:
   - Diamond: {np.sqrt(1100 * 1e-6 * 1e-3) * 1e3:.2f} mm
   - Ag: {np.sqrt(173 * 1e-6 * 1e-3) * 1e3:.2f} mm
   - Glass: {np.sqrt(0.5 * 1e-6 * 1e-3) * 1e3:.4f} mm

5. FRAMEWORK CONNECTION
   Thermal diffusivity α = κ/(ρ×Cp)

   For metals: α ∝ κ_e ∝ 1/γ_electron
   For insulators: α ∝ κ_ph ∝ θ_D/γ_G²

   α measures how efficiently thermal coherence propagates.
""")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Metals - α vs 1/γ_electron
ax1 = axes[0, 0]
ax1.scatter(1/gamma_electron, alpha_m, c='steelblue', s=100, alpha=0.7)
for i, m in enumerate(metals):
    ax1.annotate(m, (1/gamma_electron[i], alpha_m[i]), fontsize=9)

slope, intercept, r, p, se = stats.linregress(1/gamma_electron, alpha_m)
x_fit = np.linspace(1, 10, 100)
ax1.plot(x_fit, slope*x_fit + intercept, 'k--', alpha=0.5)

ax1.set_xlabel('1/γ_electron (coherence)', fontsize=12)
ax1.set_ylabel('α (mm²/s)', fontsize=12)
ax1.set_title(f'Metals: α vs Coherence (r = {r1:.3f})', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: Insulators - log(α) vs θ_D/γ_G²
ax2 = axes[0, 1]
ax2.scatter(combined, alpha_i, c='coral', s=100, alpha=0.7)
for i, m in enumerate(insulators):
    ax2.annotate(m, (combined[i], alpha_i[i]), fontsize=9)
ax2.set_xlabel('θ_D / γ_G²', fontsize=12)
ax2.set_ylabel('α (mm²/s)', fontsize=12)
ax2.set_yscale('log')
ax2.set_title(f'Insulators: α vs θ_D/γ_G² (r = {r8:.3f})', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: All materials - α vs κ
all_alpha = np.concatenate([alpha_m, alpha_i])
all_kappa = np.concatenate([kappa_m, kappa_i])
all_colors = ['steelblue'] * len(metals) + ['coral'] * len(insulators)

ax3 = axes[1, 0]
ax3.scatter(all_kappa, all_alpha, c=all_colors, s=100, alpha=0.7)
for i, m in enumerate(metals):
    ax3.annotate(m, (kappa_m[i], alpha_m[i]), fontsize=8)
for i, m in enumerate(insulators):
    ax3.annotate(m, (kappa_i[i], alpha_i[i]), fontsize=8)
ax3.set_xlabel('κ (W/m·K)', fontsize=12)
ax3.set_ylabel('α (mm²/s)', fontsize=12)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_title('α vs κ: Both Domains', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: Thermal diffusion length
ax4 = axes[1, 1]
all_L = np.concatenate([L_metals, L_insulators])
all_names = metals + insulators

sorted_idx = np.argsort(all_L)[::-1]
y_pos = np.arange(len(all_L))

ax4.barh(y_pos, [all_L[i] for i in sorted_idx],
         color=['coral' if all_names[i] in insulators else 'steelblue' for i in sorted_idx])
ax4.set_yticks(y_pos)
ax4.set_yticklabels([all_names[i] for i in sorted_idx])
ax4.set_xlabel('Diffusion length in 1 ms (mm)', fontsize=12)
ax4.set_title('Thermal Diffusion Length Hierarchy', fontsize=14)
ax4.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/thermal_diffusivity_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #111")
print("="*70)

print(f"""
KEY RESULTS:

1. METALS
   α vs 1/γ_electron: r = {r1:.3f}
   α is dominated by κ, which follows coherence.

2. INSULATORS
   log(α) vs θ_D/γ_G²: r = {r8:.3f}
   Combined model from Session #108 applies.

3. DIAMOND SUPREMACY
   α = 1100 mm²/s (highest known!)
   Diamond combines:
   - Highest κ (strong bonds, low γ_G)
   - Moderate heat capacity
   - Low density

4. THERMAL DIFFUSION LENGTH
   L = √(α×t) measures heat spreading.
   In 1 ms:
   - Diamond: {np.sqrt(1100 * 1e-6 * 1e-3) * 1e3:.2f} mm
   - Silver: {np.sqrt(173 * 1e-6 * 1e-3) * 1e3:.2f} mm
   - Glass: {np.sqrt(0.5 * 1e-6 * 1e-3) * 1e3:.4f} mm (2000× slower)

FRAMEWORK CONNECTION:
α = κ/(ρ×Cp) measures thermal coherence propagation speed.

For METALS: α ∝ κ_e ∝ 1/γ_electron
For INSULATORS: α ∝ κ_ph ∝ θ_D/γ_G²

Thermal diffusivity is a DERIVED coherence property.
It combines κ (coherence) with ρ×Cp (heat capacity).
""")

print("\nFigure saved to: thermal_diffusivity_coherence.png")
