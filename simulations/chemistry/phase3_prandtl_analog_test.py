#!/usr/bin/env python3
"""
Phase 3 Session #1: Prandtl Analog Test

If coherence channels behave like viscosity modes in a multi-component N-S
fluid, then the RATIOS between channel viscosities (Prandtl analogs) should
be approximately constant within material classes, even though absolute
values vary widely.

In fluid dynamics:
  Pr = ν/α ≈ 0.71 for all ideal gases (remarkable!)
  Pm = ν/η varies by ~10^20 across materials (not a universal constant)

The test: Does γ_phonon/γ_electron behave more like Pr (class-invariant)
or more like raw γ (unconstrained)?

Failure criterion: If the ratio varies as much WITHIN a class as BETWEEN
classes, the Prandtl analog hypothesis fails.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("PHASE 3 SESSION #1: PRANDTL ANALOG TEST")
print("Channel Ratio Stability Within Material Classes")
print("=" * 70)

# ==============================================================================
# Extended dataset: 3d transition metals, noble metals, post-transition metals
# All at T = 300K
# ==============================================================================

T = 300  # K

# Data: material, class, θ_D(K), σ(10^7 S/m), λ_ep
# θ_D from standard tables, σ from CRC Handbook, λ_ep from McMillan/Allen-Dynes
materials = {
    # 3d transition metals
    'Ti': ('3d', 420, 0.238, 0.38),
    'V':  ('3d', 380, 0.500, 0.82),
    'Cr': ('3d', 630, 0.774, 0.40),
    'Mn': ('3d', 410, 0.069, 0.60),  # α-Mn, poor conductor
    'Fe': ('3d', 470, 1.000, 0.40),
    'Co': ('3d', 445, 1.720, 0.38),
    'Ni': ('3d', 450, 1.430, 0.60),
    'Cu': ('3d', 343, 5.960, 0.13),
    'Zn': ('3d', 327, 1.690, 0.38),

    # Noble metals
    'Cu_n': ('noble', 343, 5.960, 0.13),
    'Ag':   ('noble', 225, 6.300, 0.12),
    'Au':   ('noble', 165, 4.520, 0.17),

    # 4d transition metals
    'Zr': ('4d', 291, 0.238, 0.41),
    'Nb': ('4d', 275, 0.693, 1.26),
    'Mo': ('4d', 450, 1.870, 0.41),
    'Ru': ('4d', 600, 1.350, 0.38),
    'Rh': ('4d', 480, 2.080, 0.30),
    'Pd': ('4d', 274, 0.953, 0.40),

    # 5d transition metals
    'Hf': ('5d', 252, 0.312, 0.34),
    'Ta': ('5d', 258, 0.761, 0.65),
    'W':  ('5d', 400, 1.790, 0.28),
    'Re': ('5d', 430, 0.543, 0.46),
    'Os': ('5d', 500, 1.100, 0.27),
    'Ir': ('5d', 420, 1.960, 0.34),
    'Pt': ('5d', 240, 0.952, 0.66),

    # Post-transition (soft, heavy)
    'In': ('post', 108, 1.190, 0.81),
    'Sn': ('post', 200, 0.917, 0.72),
    'Pb': ('post', 105, 0.481, 1.55),
    'Bi': ('post', 119, 0.086, 0.30),  # semimetal

    # Alkali metals
    'Na': ('alkali', 158, 2.100, 0.17),
    'K':  ('alkali', 91,  1.390, 0.13),
    'Rb': ('alkali', 56,  0.800, 0.15),
    'Cs': ('alkali', 38,  0.500, 0.16),
}

# ==============================================================================
# Compute channel-specific γ and ratios
# ==============================================================================

print(f"\n{'Mat':<5} {'Class':<8} {'θ_D':<6} {'σ(10⁷)':<9} {'γ_ph':<7} {'γ_el':<7} {'Ratio':<8} {'λ_ep':<6}")
print("-" * 65)

classes = {}
for name, (cls, theta_D, sigma_7, lam) in materials.items():
    gamma_phonon = 2 * T / theta_D
    # γ_electron: normalize conductivity to define a "disorder" parameter
    # Use resistivity ρ = 1/σ as analog; normalize by typical metallic ρ_0
    # Higher resistivity = more electron scattering = higher γ_electron
    rho = 1 / (sigma_7 * 1e7)  # Ω·m
    rho_0 = 1e-7  # typical metallic resistivity scale
    gamma_electron = rho / rho_0  # dimensionless electron disorder

    ratio = gamma_phonon / gamma_electron if gamma_electron > 0 else np.nan

    display_name = name.replace('_n', '*')
    print(f"{display_name:<5} {cls:<8} {theta_D:<6} {sigma_7:<9.3f} {gamma_phonon:<7.2f} {gamma_electron:<7.2f} {ratio:<8.2f} {lam:<6.2f}")

    if cls not in classes:
        classes[cls] = {'names': [], 'gamma_ph': [], 'gamma_el': [], 'ratio': [], 'lambda': []}
    classes[cls]['names'].append(name)
    classes[cls]['gamma_ph'].append(gamma_phonon)
    classes[cls]['gamma_el'].append(gamma_electron)
    classes[cls]['ratio'].append(ratio)
    classes[cls]['lambda'].append(lam)

# ==============================================================================
# Within-class vs between-class variance of RATIO
# ==============================================================================
print("\n" + "=" * 70)
print("PRANDTL ANALOG TEST: RATIO STABILITY WITHIN CLASSES")
print("=" * 70)

print(f"\n{'Class':<10} {'N':<4} {'Mean γ_ph':<10} {'Std γ_ph':<10} {'CV γ_ph':<9} {'Mean ratio':<12} {'Std ratio':<10} {'CV ratio':<9}")
print("-" * 80)

class_means_ratio = []
class_cvs_ratio = []
class_cvs_gamma = []
all_ratios = []

for cls in ['3d', '4d', '5d', 'noble', 'post', 'alkali']:
    if cls not in classes:
        continue
    data = classes[cls]
    gp = np.array(data['gamma_ph'])
    ratios = np.array(data['ratio'])
    all_ratios.extend(ratios)

    mean_gp = np.mean(gp)
    std_gp = np.std(gp)
    cv_gp = std_gp / mean_gp if mean_gp > 0 else np.nan

    mean_r = np.mean(ratios)
    std_r = np.std(ratios)
    cv_r = std_r / mean_r if mean_r > 0 else np.nan

    class_means_ratio.append(mean_r)
    class_cvs_ratio.append(cv_r)
    class_cvs_gamma.append(cv_gp)

    print(f"{cls:<10} {len(gp):<4} {mean_gp:<10.2f} {std_gp:<10.2f} {cv_gp:<9.2f} {mean_r:<12.2f} {std_r:<10.2f} {cv_r:<9.2f}")

# The key test: Is within-class CV of ratio < within-class CV of raw gamma?
mean_cv_ratio = np.mean(class_cvs_ratio)
mean_cv_gamma = np.mean(class_cvs_gamma)

print(f"\nMean within-class CV of γ_phonon:            {mean_cv_gamma:.3f}")
print(f"Mean within-class CV of γ_ph/γ_el ratio:     {mean_cv_ratio:.3f}")
print(f"Ratio stability improvement:                  {mean_cv_gamma/mean_cv_ratio:.2f}×")

# Between-class variance of ratio
between_class_cv = np.std(class_means_ratio) / np.mean(class_means_ratio)
print(f"\nBetween-class CV of ratio:                   {between_class_cv:.3f}")

# ==============================================================================
# F-test: within-class vs between-class variance
# ==============================================================================
print("\n" + "=" * 70)
print("ANOVA-STYLE TEST: DO CLASSES CLUSTER IN RATIO SPACE?")
print("=" * 70)

# Collect all ratios by class for ANOVA
class_groups = []
class_labels = []
for cls in ['3d', '4d', '5d', 'noble', 'post', 'alkali']:
    if cls in classes:
        class_groups.append(np.array(classes[cls]['ratio']))
        class_labels.append(cls)

# One-way ANOVA
f_stat, p_val = stats.f_oneway(*class_groups)
print(f"\nOne-way ANOVA: F = {f_stat:.2f}, p = {p_val:.4f}")
print(f"{'Classes have SIGNIFICANTLY different ratios (p < 0.05)' if p_val < 0.05 else 'Classes do NOT have significantly different ratios (p >= 0.05)'}")

if p_val < 0.05:
    print("\n→ SUPPORTS Prandtl analog: ratios cluster by class")
    print("  (Different material classes have characteristic channel-viscosity ratios)")
else:
    print("\n→ FAILS Prandtl analog: ratios do not cluster by class")
    print("  (The ratio is as variable within classes as between them)")

# ==============================================================================
# Comparison: How variable are real Prandtl numbers?
# ==============================================================================
print("\n" + "=" * 70)
print("COMPARISON: REAL FLUID PRANDTL NUMBERS")
print("=" * 70)

print("""
For reference, real Prandtl numbers Pr = ν/α (kinematic viscosity / thermal diffusivity):

  Ideal gases:    Pr ≈ 0.71 ± 0.02  (CV ≈ 0.03)  ← remarkably constant!
  Liquid metals:  Pr ≈ 0.01-0.03    (CV ≈ 0.5)    ← class-constant
  Water at 20°C:  Pr ≈ 7.0
  Glycerol:       Pr ≈ 12,500
  Molten salt:    Pr ≈ 3-10

Key insight: Pr is constant WITHIN a class (all gases, all liquid metals)
but varies enormously BETWEEN classes (gases vs liquids vs molten salts).

This is exactly the multi-component N-S prediction: transport coefficient
ratios are determined by the dominant scattering mechanism, which is
class-dependent but material-independent within a class.
""")

# ==============================================================================
# Visualization
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Phase 3 Session #1: Prandtl Analog Test\nγ_phonon/γ_electron Ratio Stability Within Material Classes',
             fontsize=13, fontweight='bold')

# Plot 1: Raw γ_phonon by class
ax = axes[0, 0]
colors = {'3d': 'blue', '4d': 'green', '5d': 'purple', 'noble': 'gold', 'post': 'red', 'alkali': 'cyan'}
for cls in ['3d', '4d', '5d', 'noble', 'post', 'alkali']:
    if cls in classes:
        data = classes[cls]
        y = data['gamma_ph']
        x = [cls] * len(y)
        ax.scatter(x, y, c=colors.get(cls, 'gray'), s=80, alpha=0.7, edgecolors='k', linewidths=0.5)
ax.set_ylabel('γ_phonon = 2T/θ_D')
ax.set_title('Raw γ_phonon by Class')
ax.grid(True, alpha=0.3)

# Plot 2: Ratio by class (the Prandtl analog test)
ax = axes[0, 1]
for cls in ['3d', '4d', '5d', 'noble', 'post', 'alkali']:
    if cls in classes:
        data = classes[cls]
        y = data['ratio']
        x = [cls] * len(y)
        ax.scatter(x, y, c=colors.get(cls, 'gray'), s=80, alpha=0.7, edgecolors='k', linewidths=0.5)
ax.set_ylabel('γ_phonon / γ_electron')
ax.set_title(f'Prandtl Analog Ratio by Class\nANOVA: F={f_stat:.1f}, p={p_val:.4f}')
ax.grid(True, alpha=0.3)

# Plot 3: CV comparison (ratio vs raw)
ax = axes[1, 0]
x_pos = np.arange(len(class_labels))
width = 0.35
bars1 = ax.bar(x_pos - width/2, class_cvs_gamma, width, label='CV(γ_phonon)', color='steelblue', alpha=0.7)
bars2 = ax.bar(x_pos + width/2, class_cvs_ratio, width, label='CV(γ_ph/γ_el)', color='coral', alpha=0.7)
ax.set_xticks(x_pos)
ax.set_xticklabels(class_labels)
ax.set_ylabel('Coefficient of Variation')
ax.set_title('Within-Class Variability: Raw γ vs Ratio')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

# Plot 4: γ_phonon vs γ_electron scatter, colored by class
ax = axes[1, 1]
for cls in ['3d', '4d', '5d', 'noble', 'post', 'alkali']:
    if cls in classes:
        data = classes[cls]
        ax.scatter(data['gamma_el'], data['gamma_ph'], c=colors.get(cls, 'gray'),
                  s=80, alpha=0.7, label=cls, edgecolors='k', linewidths=0.5)
        # Draw constant-ratio lines for each class mean
        mean_ratio = np.mean(data['ratio'])
        x_range = np.linspace(0.01, max(max(d['gamma_el']) for d in classes.values()) * 1.1, 50)
        ax.plot(x_range, mean_ratio * x_range, '--', c=colors.get(cls, 'gray'), alpha=0.3)

ax.set_xlabel('γ_electron (resistivity-based)')
ax.set_ylabel('γ_phonon = 2T/θ_D')
ax.set_title('Channel Phase Space (dashed = class mean ratio)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phase3_prandtl_analog_test.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nFigure saved: phase3_prandtl_analog_test.png")

# ==============================================================================
# Final Assessment
# ==============================================================================
print("\n" + "=" * 70)
print("FINAL ASSESSMENT")
print("=" * 70)

if p_val < 0.05 and mean_cv_ratio < mean_cv_gamma:
    verdict = "PARTIALLY SUPPORTED"
    detail = ("Ratios cluster by class AND are more stable than raw γ values.\n"
              "This is consistent with the Prandtl analog hypothesis: material classes\n"
              "have characteristic channel-viscosity ratios, as predicted by multi-component N-S.")
elif p_val < 0.05:
    verdict = "MIXED"
    detail = ("Ratios cluster by class but are not more stable than raw γ.\n"
              "The clustering may reflect trivial correlations rather than N-S structure.")
else:
    verdict = "FAILED"
    detail = ("Ratios do NOT cluster by class. The γ_ph/γ_el ratio varies as much\n"
              "within classes as between them. The Prandtl analog hypothesis is not supported.")

print(f"\nVerdict: {verdict}")
print(f"\n{detail}")

print(f"""
KEY NUMBERS:
  Within-class CV of raw γ_phonon:  {mean_cv_gamma:.3f}
  Within-class CV of ratio:         {mean_cv_ratio:.3f}
  ANOVA F-statistic:                {f_stat:.2f}
  ANOVA p-value:                    {p_val:.4f}
  Between-class CV of mean ratio:   {between_class_cv:.3f}

INTERPRETATION:
  If CV(ratio) << CV(γ_phonon): ratio is more fundamental than raw γ
  If p < 0.05: classes have characteristic ratios (Prandtl-like)
  If between-class CV >> within-class CV: strong class signature

STATUS: {'Validated' if verdict == 'PARTIALLY SUPPORTED' else 'Untested' if verdict == 'MIXED' else 'Failed'} — {verdict}
""")
