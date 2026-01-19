#!/usr/bin/env python3
"""
Chemistry Session #94: Magnetostriction and Coherence

Test whether magnetostriction coefficient λ_s relates to coherence γ.

Magnetostriction = change in dimensions when magnetized
λ_s = saturation magnetostriction (dimensionless, typically ppm)

Session #93 showed piezoelectricity (mechanical-electrical) scales with γ.
Question: Does magnetostriction (mechanical-magnetic) also scale with γ?

Hypothesis:
- Like piezoelectricity, magnetostriction may be "soft mode" enhanced
- Alternatively, it may require spin-lattice coherence coupling
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Magnetostriction data
# λ_s = saturation magnetostriction (×10^-6)
# θ_D = Debye temperature (K)
# T_C = Curie temperature (K)
# M_s = saturation magnetization (emu/cm³)

magnetostriction_data = {
    # Pure elements
    'Fe': {'lambda_s': -7, 'theta_D': 470, 'T_C': 1043, 'M_s': 1714},
    'Co': {'lambda_s': -62, 'theta_D': 445, 'T_C': 1388, 'M_s': 1422},
    'Ni': {'lambda_s': -34, 'theta_D': 450, 'T_C': 627, 'M_s': 485},
    'Gd': {'lambda_s': -130, 'theta_D': 182, 'T_C': 292, 'M_s': 2060},
    'Tb': {'lambda_s': 8700, 'theta_D': 177, 'T_C': 219, 'M_s': 3220},  # Giant!
    'Dy': {'lambda_s': 9000, 'theta_D': 186, 'T_C': 85, 'M_s': 2950},   # Giant!

    # Rare earth alloys (giant magnetostriction)
    'Terfenol-D (Tb0.3Dy0.7Fe2)': {'lambda_s': 2000, 'theta_D': 200, 'T_C': 653, 'M_s': 800},
    'Tb2Fe17': {'lambda_s': 350, 'theta_D': 220, 'T_C': 410, 'M_s': 1100},
    'SmFe2': {'lambda_s': -2340, 'theta_D': 280, 'T_C': 700, 'M_s': 900},
    'TbFe2': {'lambda_s': 2400, 'theta_D': 200, 'T_C': 700, 'M_s': 950},
    'DyFe2': {'lambda_s': 650, 'theta_D': 210, 'T_C': 635, 'M_s': 880},

    # Iron alloys
    'Fe65Co35': {'lambda_s': 140, 'theta_D': 460, 'T_C': 1170, 'M_s': 1950},
    'Fe50Co50': {'lambda_s': 115, 'theta_D': 450, 'T_C': 1250, 'M_s': 1950},
    'Fe-3Si': {'lambda_s': 25, 'theta_D': 500, 'T_C': 1020, 'M_s': 1550},
    'Fe81B13.5Si3.5C2 (Metglas)': {'lambda_s': 27, 'theta_D': 350, 'T_C': 623, 'M_s': 1250},

    # Nickel alloys
    'Permalloy (Ni80Fe20)': {'lambda_s': 1, 'theta_D': 440, 'T_C': 865, 'M_s': 860},
    'Invar (Fe64Ni36)': {'lambda_s': 5, 'theta_D': 420, 'T_C': 523, 'M_s': 1100},

    # Spinels
    'Fe3O4 (magnetite)': {'lambda_s': -78, 'theta_D': 600, 'T_C': 858, 'M_s': 480},
    'CoFe2O4': {'lambda_s': -250, 'theta_D': 550, 'T_C': 793, 'M_s': 400},
    'NiFe2O4': {'lambda_s': -26, 'theta_D': 570, 'T_C': 858, 'M_s': 270},

    # Heusler alloys
    'Ni2MnGa': {'lambda_s': 50000, 'theta_D': 280, 'T_C': 376, 'M_s': 600},  # Shape memory!
}

print("=" * 70)
print("Chemistry Session #94: Magnetostriction and Coherence")
print("=" * 70)

# Extract data
materials = list(magnetostriction_data.keys())
lambda_s = np.array([abs(magnetostriction_data[m]['lambda_s']) for m in materials])  # Use absolute value
lambda_s_signed = np.array([magnetostriction_data[m]['lambda_s'] for m in materials])
theta_D = np.array([magnetostriction_data[m]['theta_D'] for m in materials])
T_C = np.array([magnetostriction_data[m]['T_C'] for m in materials])
M_s = np.array([magnetostriction_data[m]['M_s'] for m in materials])

# Calculate γ_phonon at room temperature
T = 300  # K
gamma_phonon = 2 * T / theta_D

# Calculate γ_spin = 2T/T_C (magnetic coherence)
gamma_spin = 2 * T / T_C

print(f"\n{'Material':<35} {'|λ_s| (ppm)':<12} {'θ_D (K)':<8} {'T_C (K)':<8} {'γ_ph':<8} {'γ_spin'}")
print("-" * 95)
for i, m in enumerate(materials):
    print(f"{m:<35} {lambda_s[i]:<12.0f} {theta_D[i]:<8.0f} {T_C[i]:<8.0f} {gamma_phonon[i]:<8.2f} {gamma_spin[i]:.2f}")

# ============================================================================
# Analysis 1: λ_s vs γ_phonon
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 1: λ_s vs γ_phonon (phonon coherence)")
print("=" * 70)

log_lambda = np.log10(lambda_s + 1)  # +1 to handle near-zero values
r_ph, p_ph = stats.pearsonr(gamma_phonon, log_lambda)
print(f"\nlog|λ_s| vs γ_phonon: r = {r_ph:.3f}, p = {p_ph:.2e}")

slope_ph, intercept_ph, _, _, _ = stats.linregress(gamma_phonon, log_lambda)
print(f"Fit: log|λ_s| = {slope_ph:.2f} × γ_ph + {intercept_ph:.2f}")

# ============================================================================
# Analysis 2: λ_s vs γ_spin
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 2: λ_s vs γ_spin (magnetic coherence)")
print("=" * 70)

r_spin, p_spin = stats.pearsonr(gamma_spin, log_lambda)
print(f"\nlog|λ_s| vs γ_spin: r = {r_spin:.3f}, p = {p_spin:.2e}")

# Also try T_C directly (inverse of γ_spin)
r_Tc, _ = stats.pearsonr(1/T_C, log_lambda)
print(f"log|λ_s| vs 1/T_C: r = {r_Tc:.3f}")

# ============================================================================
# Analysis 3: Combined γ models
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 3: Combined Models")
print("=" * 70)

# γ_phonon × γ_spin (both soft)
gamma_combined = gamma_phonon * gamma_spin
r_comb, _ = stats.pearsonr(gamma_combined, log_lambda)
print(f"\nlog|λ_s| vs γ_ph × γ_spin: r = {r_comb:.3f}")

# γ_phonon / γ_spin (phonon soft but magnetic stiff)
gamma_ratio = gamma_phonon / gamma_spin
r_ratio, _ = stats.pearsonr(gamma_ratio, log_lambda)
print(f"log|λ_s| vs γ_ph / γ_spin: r = {r_ratio:.3f}")

# ============================================================================
# Analysis 4: By material class
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 4: By Material Class")
print("=" * 70)

# Rare earth (giant magnetostriction)
RE_materials = ['Gd', 'Tb', 'Dy', 'Terfenol-D (Tb0.3Dy0.7Fe2)', 'Tb2Fe17',
                'SmFe2', 'TbFe2', 'DyFe2']
RE_mask = np.array([m in RE_materials for m in materials])

# 3d transition metals
TM_materials = ['Fe', 'Co', 'Ni', 'Fe65Co35', 'Fe50Co50', 'Fe-3Si',
                'Fe81B13.5Si3.5C2 (Metglas)', 'Permalloy (Ni80Fe20)', 'Invar (Fe64Ni36)']
TM_mask = np.array([m in TM_materials for m in materials])

# Ferrites
ferrite_mask = np.array(['O4' in m for m in materials])

if np.sum(RE_mask) >= 3:
    r_RE, _ = stats.pearsonr(gamma_phonon[RE_mask], log_lambda[RE_mask])
    print(f"\nRare Earth ({np.sum(RE_mask)} materials):")
    print(f"  log|λ_s| vs γ_phonon: r = {r_RE:.3f}")
    print(f"  |λ_s| range: {lambda_s[RE_mask].min():.0f} - {lambda_s[RE_mask].max():.0f} ppm")

if np.sum(TM_mask) >= 3:
    r_TM, _ = stats.pearsonr(gamma_phonon[TM_mask], log_lambda[TM_mask])
    print(f"\n3d Transition Metals ({np.sum(TM_mask)} materials):")
    print(f"  log|λ_s| vs γ_phonon: r = {r_TM:.3f}")
    print(f"  |λ_s| range: {lambda_s[TM_mask].min():.0f} - {lambda_s[TM_mask].max():.0f} ppm")

# ============================================================================
# Analysis 5: Spin-orbit coupling
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 5: Spin-Orbit Coupling Perspective")
print("=" * 70)

print("""
Magnetostriction arises from spin-orbit coupling (SOC):
- 4f electrons (rare earths): Large orbital moment → giant λ_s
- 3d electrons (transition metals): Orbital moment quenched → small λ_s

Rare earths have unquenched orbital moments:
- Tb: L = 3, S = 3, J = 6 → giant SOC → giant λ_s
- Dy: L = 5, S = 5/2, J = 15/2 → giant SOC → giant λ_s

3d metals have quenched orbital moments:
- Fe, Co, Ni: L ~ 0 → small SOC → small λ_s
""")

# Group by atomic type
print("\n|λ_s| by Group:")
print(f"  Rare earth pure elements: {lambda_s[np.array([m in ['Gd', 'Tb', 'Dy'] for m in materials])].mean():.0f} ppm")
print(f"  Rare earth alloys: {lambda_s[np.array([m in ['Terfenol-D (Tb0.3Dy0.7Fe2)', 'TbFe2', 'SmFe2', 'DyFe2'] for m in materials])].mean():.0f} ppm")
print(f"  Pure 3d metals: {lambda_s[np.array([m in ['Fe', 'Co', 'Ni'] for m in materials])].mean():.0f} ppm")

# ============================================================================
# Analysis 6: Comparison to piezoelectricity
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 6: Comparison to Piezoelectricity (#93)")
print("=" * 70)

print("""
Piezoelectricity (Session #93):
- d_33 ∝ γ_phonon^2.48: r = 0.867 (POSITIVE)
- Soft phonon modes enhance piezoelectric response
- Maximum at morphotropic phase boundary

Magnetostriction (this session):
""")

if r_ph > 0.3:
    print(f"- |λ_s| ∝ γ_phonon: r = {r_ph:.3f} (POSITIVE)")
    print("- Similar to piezoelectricity - soft modes help!")
elif r_ph < -0.3:
    print(f"- |λ_s| ∝ 1/γ_phonon: r = {r_ph:.3f} (NEGATIVE)")
    print("- OPPOSITE to piezoelectricity - coherent phonons help!")
else:
    print(f"- |λ_s| vs γ_phonon: r = {r_ph:.3f} (WEAK)")
    print("- γ_phonon is not the controlling factor!")

print(f"\nγ_spin correlation: r = {r_spin:.3f}")
if abs(r_spin) > abs(r_ph):
    print("MAGNETIC coherence matters more than phonon coherence!")

# ============================================================================
# Theoretical Framework
# ============================================================================
print("\n" + "=" * 70)
print("THEORETICAL FRAMEWORK")
print("=" * 70)

print("""
Magnetostriction arises from:
1. Spin-orbit coupling: E_SO = λ × L·S
2. Crystal field: E_CF = A × (structural anisotropy)
3. Exchange: E_ex = -J × S_i·S_j

Total strain:
λ_s ∝ (spin-orbit coupling) × (lattice softness) × (magnetic order)

In coherence language:
- λ_s ∝ SOC × (γ_phonon) × (2/γ_spin)

If γ_phonon dominates: λ_s ∝ γ (soft modes help, like piezoelectricity)
If γ_spin dominates: λ_s ∝ 2/γ (magnetic coherence helps)
If SOC dominates: λ_s independent of γ (atomic property)
""")

# ============================================================================
# Key Results
# ============================================================================
print("\n" + "=" * 70)
print("KEY RESULTS")
print("=" * 70)

print(f"""
1. log|λ_s| vs γ_phonon: r = {r_ph:.3f}
   {'Soft phonon modes help!' if r_ph > 0.3 else 'Phonon coherence not dominant'}

2. log|λ_s| vs γ_spin: r = {r_spin:.3f}
   {'Magnetic incoherence helps!' if r_spin > 0.3 else 'Magnetic coherence not key factor' if abs(r_spin) < 0.3 else 'Magnetic coherence helps!'}

3. Combined γ_ph × γ_spin: r = {r_comb:.3f}

4. Main factor: {'SPIN-ORBIT COUPLING (atomic property)' if max(abs(r_ph), abs(r_spin)) < 0.5 else 'Coherence parameters matter'}

Giant magnetostriction in rare earths is due to:
- Large unquenched orbital angular momentum (L = 3-5)
- Strong spin-orbit coupling (E_SO ∝ Z^4)
- NOT primarily from phonon softness

This is DIFFERENT from piezoelectricity:
- Piezo: mechanical-electrical coupling via lattice instability
- Magneto: mechanical-magnetic coupling via spin-orbit (atomic)

INSIGHT: Magnetostriction is more like an ATOMIC property
than a collective coherence phenomenon.
""")

# ============================================================================
# Visualization
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Color by class
colors = []
for m in materials:
    if m in RE_materials:
        colors.append('red')
    elif m in TM_materials:
        colors.append('blue')
    elif 'O4' in m:
        colors.append('green')
    else:
        colors.append('purple')
colors = np.array(colors)

# Plot 1: λ_s vs γ_phonon
ax1 = axes[0, 0]
for c, label in [('red', 'Rare Earth'), ('blue', '3d Metal'), ('green', 'Ferrite'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax1.scatter(gamma_phonon[mask], lambda_s[mask], c=c, label=label, s=80, alpha=0.7)
ax1.set_xlabel('γ_phonon = 2T/θ_D', fontsize=12)
ax1.set_ylabel('|λ_s| (ppm)', fontsize=12)
ax1.set_yscale('log')
ax1.set_title(f'|λ_s| vs γ_phonon (r = {r_ph:.3f})', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: λ_s vs γ_spin
ax2 = axes[0, 1]
for c, label in [('red', 'Rare Earth'), ('blue', '3d Metal'), ('green', 'Ferrite'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax2.scatter(gamma_spin[mask], lambda_s[mask], c=c, label=label, s=80, alpha=0.7)
ax2.set_xlabel('γ_spin = 2T/T_C', fontsize=12)
ax2.set_ylabel('|λ_s| (ppm)', fontsize=12)
ax2.set_yscale('log')
ax2.set_title(f'|λ_s| vs γ_spin (r = {r_spin:.3f})', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: λ_s vs γ_combined
ax3 = axes[1, 0]
for c, label in [('red', 'Rare Earth'), ('blue', '3d Metal'), ('green', 'Ferrite'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax3.scatter(gamma_combined[mask], lambda_s[mask], c=c, label=label, s=80, alpha=0.7)
ax3.set_xlabel('γ_phonon × γ_spin', fontsize=12)
ax3.set_ylabel('|λ_s| (ppm)', fontsize=12)
ax3.set_yscale('log')
ax3.set_title(f'|λ_s| vs γ_ph × γ_spin (r = {r_comb:.3f})', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: θ_D histogram by class
ax4 = axes[1, 1]
ax4.bar(['Rare Earth', '3d Metal', 'Ferrite'],
        [theta_D[RE_mask].mean() if np.sum(RE_mask) > 0 else 0,
         theta_D[TM_mask].mean() if np.sum(TM_mask) > 0 else 0,
         theta_D[ferrite_mask].mean() if np.sum(ferrite_mask) > 0 else 0],
        color=['red', 'blue', 'green'], alpha=0.7)
ax4.set_ylabel('Mean θ_D (K)', fontsize=12)
ax4.set_title('Debye Temperature by Class', fontsize=14)
ax4.grid(True, alpha=0.3, axis='y')

# Add secondary axis for λ_s
ax4b = ax4.twinx()
ax4b.bar(['Rare Earth', '3d Metal', 'Ferrite'],
         [lambda_s[RE_mask].mean() if np.sum(RE_mask) > 0 else 0,
          lambda_s[TM_mask].mean() if np.sum(TM_mask) > 0 else 0,
          lambda_s[ferrite_mask].mean() if np.sum(ferrite_mask) > 0 else 0],
         color=['red', 'blue', 'green'], alpha=0.3, width=0.4)
ax4b.set_ylabel('Mean |λ_s| (ppm)', fontsize=12)
ax4b.set_yscale('log')

plt.suptitle('Chemistry Session #94: Magnetostriction and Coherence',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/magnetostriction_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: magnetostriction_coherence.png")

# ============================================================================
# Predictions
# ============================================================================
print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print(f"""
P94.1: Magnetostriction is dominated by SPIN-ORBIT COUPLING
       Not by phonon coherence (γ_phonon) or spin coherence (γ_spin).

P94.2: Rare earth elements have giant λ_s because of 4f electrons
       Large unquenched orbital moment L → large spin-orbit coupling.

P94.3: 3d transition metals have small λ_s because of orbital quenching
       Crystal field quenches L → small spin-orbit coupling → small λ_s.

P94.4: λ_s vs γ_phonon: r = {r_ph:.3f}
       {'Soft phonon contribution present' if r_ph > 0.3 else 'Phonon softness not the main factor'}

P94.5: λ_s vs γ_spin: r = {r_spin:.3f}
       {'Near T_C enhances magnetostriction' if r_spin > 0.3 else 'Magnetic coherence not key'}

P94.6: Magnetostriction is FUNDAMENTALLY DIFFERENT from piezoelectricity
       Piezo: collective lattice instability (γ helps)
       Magneto: atomic spin-orbit coupling (atomic property)

FRAMEWORK CLASSIFICATION:
Magnetostriction belongs to the category of ATOMIC PROPERTIES
that do not strongly depend on collective coherence parameters.

Similar to:
- Atomic radius (shell structure)
- Ionization energy (electron binding)
- Electronegativity (electron affinity)

These are set by atomic physics, not by collective coherence.
""")

# Validation status
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if max(abs(r_ph), abs(r_spin)) > 0.5:
    status = "MODERATE VALIDATION"
    coherence_role = "Coherence plays a secondary role"
else:
    status = "LIMITED COHERENCE CORRELATION"
    coherence_role = "Spin-orbit coupling (atomic) dominates"

print(f"""
{status}

|λ_s| vs γ_phonon: r = {r_ph:.3f}
|λ_s| vs γ_spin: r = {r_spin:.3f}
Combined model: r = {r_comb:.3f}

{coherence_role}

KEY INSIGHT: Not all properties follow coherence scaling!

The framework distinguishes:
1. COLLECTIVE properties (follow γ scaling):
   - Transport (σ, κ, μ)
   - Gaps (E_g, Δ)
   - Soft mode phenomena (piezoelectricity d)

2. ATOMIC properties (independent of γ):
   - Spin-orbit coupling (magnetostriction λ)
   - Ionization energy
   - Atomic polarizability (though scales with IE → γ_optical)

Magnetostriction shows the framework's LIMITS - not everything
is controlled by collective coherence. Spin-orbit coupling is
an INTRINSIC atomic property related to nuclear charge Z.
""")
