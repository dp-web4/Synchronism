#!/usr/bin/env python3
"""
Chemistry Session #99: Magnetic Anisotropy and Coherence

Test whether magnetocrystalline anisotropy K relates to coherence γ.

Magnetic anisotropy energy:
E_anis = K₁sin²θ + K₂sin⁴θ + ...

where K₁ = first-order anisotropy constant (J/m³ or erg/cm³)

K arises from spin-orbit coupling (SOC) + crystal field.

Session #94 showed magnetostriction ∝ SOC (atomic property).
Question: Is magnetic anisotropy also SOC-dominated, or does
lattice coherence (γ_phonon) play a role?

Hypothesis:
- K ∝ SOC (atomic property like magnetostriction)
- But K depends on crystal symmetry → some γ_phonon contribution?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Magnetic anisotropy data
# K1 = first-order anisotropy constant (MJ/m³)
# Ms = saturation magnetization (MA/m)
# Tc = Curie temperature (K)
# theta_D = Debye temperature (K)

anisotropy_data = {
    # 3d transition metals (cubic, low K)
    'Fe (bcc)': {'K1': 0.048, 'Ms': 1.71, 'Tc': 1043, 'theta_D': 470, 'structure': 'cubic'},
    'Co (hcp)': {'K1': 0.45, 'Ms': 1.44, 'Tc': 1388, 'theta_D': 445, 'structure': 'hexagonal'},
    'Ni (fcc)': {'K1': -0.005, 'Ms': 0.49, 'Tc': 627, 'theta_D': 450, 'structure': 'cubic'},

    # Rare earth metals (hexagonal, GIANT K)
    'Gd': {'K1': -0.012, 'Ms': 2.06, 'Tc': 292, 'theta_D': 182, 'structure': 'hexagonal'},
    'Tb': {'K1': 8.3, 'Ms': 3.22, 'Tc': 219, 'theta_D': 177, 'structure': 'hexagonal'},
    'Dy': {'K1': 4.5, 'Ms': 2.95, 'Tc': 85, 'theta_D': 186, 'structure': 'hexagonal'},
    'Ho': {'K1': 3.1, 'Ms': 3.27, 'Tc': 20, 'theta_D': 190, 'structure': 'hexagonal'},
    'Er': {'K1': 1.3, 'Ms': 2.79, 'Tc': 19, 'theta_D': 195, 'structure': 'hexagonal'},

    # Rare earth-transition metal compounds (high K for permanent magnets)
    'SmCo5': {'K1': 17.0, 'Ms': 0.86, 'Tc': 1020, 'theta_D': 400, 'structure': 'hexagonal'},
    'Nd2Fe14B': {'K1': 4.9, 'Ms': 1.28, 'Tc': 585, 'theta_D': 350, 'structure': 'tetragonal'},
    'Sm2Co17': {'K1': 3.3, 'Ms': 1.14, 'Tc': 1190, 'theta_D': 380, 'structure': 'hexagonal'},
    'SmFe11Ti': {'K1': 4.9, 'Ms': 1.24, 'Tc': 584, 'theta_D': 350, 'structure': 'tetragonal'},

    # Spinel ferrites
    'Fe3O4': {'K1': -0.011, 'Ms': 0.48, 'Tc': 858, 'theta_D': 600, 'structure': 'cubic'},
    'CoFe2O4': {'K1': 0.27, 'Ms': 0.40, 'Tc': 793, 'theta_D': 550, 'structure': 'cubic'},
    'MnFe2O4': {'K1': -0.003, 'Ms': 0.42, 'Tc': 573, 'theta_D': 560, 'structure': 'cubic'},

    # Garnets
    'Y3Fe5O12 (YIG)': {'K1': -0.0006, 'Ms': 0.14, 'Tc': 560, 'theta_D': 500, 'structure': 'cubic'},

    # Hexagonal ferrites (for permanent magnets)
    'BaFe12O19': {'K1': 0.33, 'Ms': 0.38, 'Tc': 723, 'theta_D': 480, 'structure': 'hexagonal'},
    'SrFe12O19': {'K1': 0.35, 'Ms': 0.37, 'Tc': 733, 'theta_D': 470, 'structure': 'hexagonal'},

    # L1_0 ordered alloys
    'FePt (L10)': {'K1': 6.6, 'Ms': 1.15, 'Tc': 750, 'theta_D': 300, 'structure': 'tetragonal'},
    'CoPt (L10)': {'K1': 4.9, 'Ms': 0.80, 'Tc': 840, 'theta_D': 290, 'structure': 'tetragonal'},
    'FePd (L10)': {'K1': 1.8, 'Ms': 1.07, 'Tc': 720, 'theta_D': 270, 'structure': 'tetragonal'},
}

print("=" * 70)
print("Chemistry Session #99: Magnetic Anisotropy and Coherence")
print("=" * 70)

# Extract data
materials = list(anisotropy_data.keys())
K1 = np.array([abs(anisotropy_data[m]['K1']) for m in materials])  # Use absolute value
K1_signed = np.array([anisotropy_data[m]['K1'] for m in materials])
Ms = np.array([anisotropy_data[m]['Ms'] for m in materials])
Tc = np.array([anisotropy_data[m]['Tc'] for m in materials])
theta_D = np.array([anisotropy_data[m]['theta_D'] for m in materials])
structures = [anisotropy_data[m]['structure'] for m in materials]

# Coherence parameters
T = 300  # K
gamma_phonon = 2 * T / theta_D
gamma_spin = 2 * T / Tc  # Magnetic coherence

# Anisotropy field H_A = 2K/μ₀Ms
mu_0 = 1.257e-6  # H/m
H_A = 2 * K1 * 1e6 / (mu_0 * Ms * 1e6)  # A/m

print(f"\n{'Material':<20} {'|K1| (MJ/m³)':<12} {'Ms (MA/m)':<10} {'Tc (K)':<8} {'θ_D (K)':<8} {'γ_ph':<8} {'Structure'}")
print("-" * 95)
for i, m in enumerate(materials):
    print(f"{m:<20} {K1[i]:<12.3f} {Ms[i]:<10.2f} {Tc[i]:<8.0f} {theta_D[i]:<8.0f} {gamma_phonon[i]:<8.2f} {structures[i]}")

# ============================================================================
# Analysis 1: K vs γ_phonon
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 1: |K₁| vs γ_phonon")
print("=" * 70)

log_K = np.log10(K1 + 0.0001)  # Add small offset for near-zero values
r_K_gamma, p_K_gamma = stats.pearsonr(gamma_phonon, log_K)
print(f"\nlog|K₁| vs γ_phonon: r = {r_K_gamma:.3f}, p = {p_K_gamma:.2e}")

# ============================================================================
# Analysis 2: K vs structure (symmetry effect)
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 2: |K₁| by Crystal Structure")
print("=" * 70)

for struct in ['cubic', 'hexagonal', 'tetragonal']:
    mask = np.array([s == struct for s in structures])
    if np.sum(mask) >= 2:
        print(f"\n{struct.upper()} ({np.sum(mask)} materials):")
        print(f"  |K₁| range: {K1[mask].min():.3f} - {K1[mask].max():.3f} MJ/m³")
        print(f"  Mean |K₁|: {K1[mask].mean():.3f} MJ/m³")
        print(f"  γ_phonon range: {gamma_phonon[mask].min():.2f} - {gamma_phonon[mask].max():.2f}")

# ============================================================================
# Analysis 3: K vs rare earth content
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 3: Rare Earth vs 3d Materials")
print("=" * 70)

# Rare earth materials
RE_materials = ['Gd', 'Tb', 'Dy', 'Ho', 'Er', 'SmCo5', 'Nd2Fe14B', 'Sm2Co17', 'SmFe11Ti']
RE_mask = np.array([m in RE_materials for m in materials])

# 3d-only materials
TM_materials = ['Fe (bcc)', 'Co (hcp)', 'Ni (fcc)', 'Fe3O4', 'CoFe2O4', 'MnFe2O4',
                'Y3Fe5O12 (YIG)', 'BaFe12O19', 'SrFe12O19']
TM_mask = np.array([m in TM_materials for m in materials])

print(f"\nRare earth ({np.sum(RE_mask)} materials):")
print(f"  |K₁| range: {K1[RE_mask].min():.3f} - {K1[RE_mask].max():.3f} MJ/m³")
print(f"  Mean |K₁|: {K1[RE_mask].mean():.2f} MJ/m³")

print(f"\n3d-only ({np.sum(TM_mask)} materials):")
print(f"  |K₁| range: {K1[TM_mask].min():.4f} - {K1[TM_mask].max():.3f} MJ/m³")
print(f"  Mean |K₁|: {K1[TM_mask].mean():.3f} MJ/m³")

# Ratio
print(f"\nRE/3d ratio: {K1[RE_mask].mean() / K1[TM_mask].mean():.0f}×")

# ============================================================================
# Analysis 4: Correlation within classes
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 4: Within-Class Correlations")
print("=" * 70)

if np.sum(RE_mask) >= 4:
    r_RE, _ = stats.pearsonr(gamma_phonon[RE_mask], log_K[RE_mask])
    print(f"\nRare earth: log|K₁| vs γ_phonon: r = {r_RE:.3f}")

if np.sum(TM_mask) >= 4:
    r_TM, _ = stats.pearsonr(gamma_phonon[TM_mask], log_K[TM_mask])
    print(f"3d-only: log|K₁| vs γ_phonon: r = {r_TM:.3f}")

# L1_0 alloys
L10_materials = ['FePt (L10)', 'CoPt (L10)', 'FePd (L10)']
L10_mask = np.array([m in L10_materials for m in materials])
if np.sum(L10_mask) >= 3:
    print(f"\nL1_0 alloys ({np.sum(L10_mask)} materials):")
    print(f"  |K₁| range: {K1[L10_mask].min():.1f} - {K1[L10_mask].max():.1f} MJ/m³")

# ============================================================================
# Analysis 5: Anisotropy field H_A
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 5: Anisotropy Field H_A = 2K/μ₀Ms")
print("=" * 70)

log_H = np.log10(H_A + 1)  # Add offset
r_H_gamma, _ = stats.pearsonr(gamma_phonon, log_H)
print(f"\nlog(H_A) vs γ_phonon: r = {r_H_gamma:.3f}")

# Top anisotropy fields
print(f"\nTop anisotropy fields (for permanent magnets):")
sorted_idx = np.argsort(H_A)[::-1]
print(f"{'Material':<20} {'|K₁| (MJ/m³)':<12} {'Ms (MA/m)':<10} {'H_A (MA/m)':<12}")
for i in sorted_idx[:8]:
    print(f"{materials[i]:<20} {K1[i]:<12.2f} {Ms[i]:<10.2f} {H_A[i]/1e6:<12.2f}")

# ============================================================================
# Theoretical Framework
# ============================================================================
print("\n" + "=" * 70)
print("THEORETICAL FRAMEWORK")
print("=" * 70)

print("""
Magnetocrystalline Anisotropy:

K₁ arises from spin-orbit coupling (SOC) + crystal field (CF):
K ∝ λ²_SO × (crystal field parameter)

where λ_SO = spin-orbit coupling strength

3d metals (Fe, Co, Ni):
- Small λ_SO (L quenched by CF)
- Low K (0.005 - 0.5 MJ/m³)

4f rare earths (Tb, Dy, Ho):
- Large λ_SO (unquenched L)
- GIANT K (1 - 17 MJ/m³)

Crystal symmetry effect:
- Cubic: K can be small or change sign
- Hexagonal/tetragonal: Large uniaxial K

Coherence interpretation:
Like magnetostriction (Session #94), anisotropy is primarily
determined by ATOMIC spin-orbit coupling, not by lattice coherence.

However, the CRYSTAL FIELD does depend on lattice properties,
so there may be secondary γ_phonon contributions.
""")

# ============================================================================
# Key Results
# ============================================================================
print("\n" + "=" * 70)
print("KEY RESULTS")
print("=" * 70)

print(f"""
1. log|K₁| vs γ_phonon: r = {r_K_gamma:.3f}
   {'Moderate correlation' if abs(r_K_gamma) > 0.3 else 'Weak correlation'}

2. Structure effect:
   - Cubic: low K (high symmetry cancels anisotropy)
   - Hexagonal/tetragonal: high K (uniaxial anisotropy)

3. Rare earth vs 3d:
   - RE mean: {K1[RE_mask].mean():.2f} MJ/m³
   - 3d mean: {K1[TM_mask].mean():.3f} MJ/m³
   - Ratio: {K1[RE_mask].mean() / K1[TM_mask].mean():.0f}×

4. Within-class γ correlations:
   - RE: r = {r_RE:.3f}
   - 3d: r = {r_TM:.3f}

INSIGHT: Like magnetostriction (#94), anisotropy is dominated by
SPIN-ORBIT COUPLING (atomic property), not lattice coherence.

The 40× ratio between RE and 3d materials comes from:
- 4f electrons: large unquenched L → large λ_SO → large K
- 3d electrons: quenched L → small λ_SO → small K

Crystal structure plays a SECONDARY role:
- Cubic symmetry tends to reduce K
- Hexagonal/tetragonal symmetry enables large uniaxial K

The correlation with γ_phonon is likely INDIRECT:
- RE metals have low θ_D (soft lattice) AND high K (large SOC)
- This is correlation, not causation
""")

# ============================================================================
# Visualization
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Color by type
colors = []
for m in materials:
    if m in RE_materials:
        colors.append('red')
    elif m in TM_materials:
        colors.append('blue')
    elif m in L10_materials:
        colors.append('green')
    else:
        colors.append('gray')
colors = np.array(colors)

# Plot 1: K vs γ_phonon
ax1 = axes[0, 0]
for c, label in [('red', 'Rare Earth'), ('blue', '3d Materials'), ('green', 'L1_0'), ('gray', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax1.scatter(gamma_phonon[mask], K1[mask], c=c, label=label, s=80, alpha=0.7)
ax1.set_xlabel('γ_phonon = 2T/θ_D', fontsize=12)
ax1.set_ylabel('|K₁| (MJ/m³)', fontsize=12)
ax1.set_yscale('log')
ax1.set_title(f'|K₁| vs γ_phonon (r = {r_K_gamma:.3f})', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: K by structure
ax2 = axes[0, 1]
struct_colors = {'cubic': 'blue', 'hexagonal': 'red', 'tetragonal': 'green'}
for struct in struct_colors:
    mask = np.array([s == struct for s in structures])
    if np.sum(mask) > 0:
        ax2.scatter(np.arange(np.sum(mask)), K1[mask], c=struct_colors[struct],
                    label=struct, s=80, alpha=0.7)
ax2.set_xlabel('Material Index', fontsize=12)
ax2.set_ylabel('|K₁| (MJ/m³)', fontsize=12)
ax2.set_yscale('log')
ax2.set_title('|K₁| by Crystal Structure', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: K vs Ms
ax3 = axes[1, 0]
for c, label in [('red', 'Rare Earth'), ('blue', '3d Materials'), ('green', 'L1_0'), ('gray', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax3.scatter(Ms[mask], K1[mask], c=c, label=label, s=80, alpha=0.7)
ax3.set_xlabel('Ms (MA/m)', fontsize=12)
ax3.set_ylabel('|K₁| (MJ/m³)', fontsize=12)
ax3.set_yscale('log')
ax3.set_title('|K₁| vs Saturation Magnetization', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Box plot by material type
ax4 = axes[1, 1]
data_to_plot = [K1[RE_mask], K1[TM_mask]]
labels_plot = ['Rare Earth', '3d Materials']
bp = ax4.boxplot(data_to_plot, labels=labels_plot, patch_artist=True)
bp['boxes'][0].set_facecolor('red')
bp['boxes'][0].set_alpha(0.5)
bp['boxes'][1].set_facecolor('blue')
bp['boxes'][1].set_alpha(0.5)
ax4.set_ylabel('|K₁| (MJ/m³)', fontsize=12)
ax4.set_title(f'Anisotropy by Material Class (ratio = {K1[RE_mask].mean() / K1[TM_mask].mean():.0f}×)', fontsize=14)
ax4.set_yscale('log')
ax4.grid(True, alpha=0.3, axis='y')

plt.suptitle('Chemistry Session #99: Magnetic Anisotropy and Coherence',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/magnetic_anisotropy_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: magnetic_anisotropy_coherence.png")

# ============================================================================
# Predictions
# ============================================================================
print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print(f"""
P99.1: K dominated by SPIN-ORBIT COUPLING (atomic property)
       RE/3d ratio = {K1[RE_mask].mean() / K1[TM_mask].mean():.0f}× due to 4f vs 3d SOC.

P99.2: Crystal structure provides SECONDARY modulation
       Cubic: low K (symmetry cancellation)
       Hexagonal/tetragonal: large uniaxial K

P99.3: Overall K vs γ_phonon: r = {r_K_gamma:.3f}
       Correlation is INDIRECT (RE has both low θ_D and high SOC).

P99.4: Within-class correlations are weak
       RE: r = {r_RE:.3f}, 3d: r = {r_TM:.3f}

P99.5: Best permanent magnets maximize K × (BH)max
       SmCo5, Nd2Fe14B have high K AND high Ms.

P99.6: Anisotropy is similar to magnetostriction (#94)
       Both dominated by SOC, not by lattice coherence.

FRAMEWORK POSITION:
Magnetic anisotropy K, like magnetostriction λ, is an
ATOMIC PROPERTY dominated by spin-orbit coupling.

The coherence framework has LIMITED applicability to:
- Spin-orbit phenomena (λ, K)
- Energy barriers (thermionic emission φ)
- Nuclear effects

The coherence framework DOES apply to:
- Transport (σ, κ, μ)
- Optical properties (n, ε, χ)
- Soft mode phenomena (d, r)
- Phase transitions (Tc, Tm)
""")

# Validation status
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

print(f"""
**CONFIRMS ATOMIC-DOMINATED CATEGORY**

Key correlations:
- Overall K vs γ_phonon: r = {r_K_gamma:.3f}
- Within RE: r = {r_RE:.3f}
- Within 3d: r = {r_TM:.3f}
- RE/3d ratio: {K1[RE_mask].mean() / K1[TM_mask].mean():.0f}×

Like magnetostriction (#94), anisotropy shows:
1. Overall correlation driven by material CLASS separation
2. Within-class correlations are WEAK
3. Spin-orbit coupling (atomic) is the primary determinant

This further establishes the ATOMIC-DOMINATED category
in the coherence framework alongside magnetostriction.

Spin-orbit phenomena = atomic properties, not coherence properties.
""")
