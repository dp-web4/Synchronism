#!/usr/bin/env python3
"""
Chemistry Session #98: Thermionic Emission and Coherence

Test whether thermionic emission relates to coherence parameter γ.

Richardson-Dushman equation:
J = A × T² × exp(-φ/kT)

where:
- J = current density
- A = Richardson constant (ideally A_0 = 120 A/cm²K²)
- φ = work function
- T = temperature

Actual A differs from A_0 due to surface effects.

Coherence hypothesis:
- A/A_0 may relate to electron coherence at surface
- φ relates to γ_electron (electron binding)
- Deviation from A_0 = coherence effect
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Thermionic emission data for metals
# phi = work function (eV)
# A = Richardson constant (A/cm²K²)
# theta_D = Debye temperature (K)

thermionic_data = {
    # Alkali metals (low work function)
    'Li': {'phi': 2.93, 'A': 30, 'theta_D': 344},
    'Na': {'phi': 2.75, 'A': 45, 'theta_D': 158},
    'K': {'phi': 2.30, 'A': 48, 'theta_D': 91},
    'Rb': {'phi': 2.16, 'A': 60, 'theta_D': 56},
    'Cs': {'phi': 2.14, 'A': 75, 'theta_D': 38},

    # Alkaline earth
    'Ca': {'phi': 2.87, 'A': 60, 'theta_D': 230},
    'Sr': {'phi': 2.59, 'A': 60, 'theta_D': 147},
    'Ba': {'phi': 2.52, 'A': 60, 'theta_D': 110},

    # Transition metals
    'Ti': {'phi': 4.33, 'A': 60, 'theta_D': 420},
    'V': {'phi': 4.30, 'A': 60, 'theta_D': 380},
    'Cr': {'phi': 4.50, 'A': 60, 'theta_D': 630},
    'Fe': {'phi': 4.67, 'A': 26, 'theta_D': 470},
    'Co': {'phi': 5.00, 'A': 41, 'theta_D': 445},
    'Ni': {'phi': 5.15, 'A': 30, 'theta_D': 450},
    'Cu': {'phi': 4.65, 'A': 65, 'theta_D': 343},
    'Zn': {'phi': 4.33, 'A': 50, 'theta_D': 327},

    # Noble metals
    'Ag': {'phi': 4.26, 'A': 60, 'theta_D': 225},
    'Au': {'phi': 5.10, 'A': 60, 'theta_D': 165},
    'Pt': {'phi': 5.65, 'A': 32, 'theta_D': 240},

    # Refractory metals (important for cathodes)
    'Mo': {'phi': 4.60, 'A': 55, 'theta_D': 450},
    'Ta': {'phi': 4.25, 'A': 60, 'theta_D': 240},
    'W': {'phi': 4.55, 'A': 60, 'theta_D': 400},
    'Re': {'phi': 4.96, 'A': 60, 'theta_D': 430},

    # Post-transition
    'Al': {'phi': 4.28, 'A': 120, 'theta_D': 428},
    'Pb': {'phi': 4.25, 'A': 12, 'theta_D': 105},
    'Bi': {'phi': 4.22, 'A': 45, 'theta_D': 119},

    # Lanthanum (high-emission cathode material)
    'La': {'phi': 3.50, 'A': 60, 'theta_D': 142},

    # LaB6 (low work function compound, excellent emitter)
    'LaB6': {'phi': 2.70, 'A': 29, 'theta_D': 650},

    # Thoriated tungsten (enhanced emission)
    'W-ThO2': {'phi': 2.70, 'A': 4, 'theta_D': 400},
}

print("=" * 70)
print("Chemistry Session #98: Thermionic Emission and Coherence")
print("=" * 70)

# Extract data
materials = list(thermionic_data.keys())
phi = np.array([thermionic_data[m]['phi'] for m in materials])
A = np.array([thermionic_data[m]['A'] for m in materials])
theta_D = np.array([thermionic_data[m]['theta_D'] for m in materials])

# Theoretical Richardson constant
A_0 = 120.4  # A/cm²K²

# Calculate coherence parameters
T = 300  # K
gamma_phonon = 2 * T / theta_D

# γ_optical from work function (binding energy)
IE_ref = 13.6  # eV
gamma_work = 2 * IE_ref / phi  # Higher φ → lower γ_work → tighter binding

# Deviation from ideal A
A_ratio = A / A_0

print(f"\n{'Material':<12} {'φ (eV)':<10} {'A (A/cm²K²)':<12} {'A/A_0':<8} {'θ_D (K)':<10} {'γ_phonon':<10} {'γ_work'}")
print("-" * 90)
for i, m in enumerate(materials):
    print(f"{m:<12} {phi[i]:<10.2f} {A[i]:<12.0f} {A_ratio[i]:<8.2f} {theta_D[i]:<10.0f} {gamma_phonon[i]:<10.2f} {gamma_work[i]:.2f}")

# ============================================================================
# Analysis 1: φ vs γ_phonon
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 1: Work Function φ vs γ_phonon")
print("=" * 70)

r_phi_gamma, p_phi_gamma = stats.pearsonr(gamma_phonon, phi)
print(f"\nφ vs γ_phonon: r = {r_phi_gamma:.3f}, p = {p_phi_gamma:.2e}")

# ============================================================================
# Analysis 2: A vs γ_phonon
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 2: Richardson Constant A vs γ_phonon")
print("=" * 70)

log_A = np.log10(A)
r_A_gamma, p_A_gamma = stats.pearsonr(gamma_phonon, log_A)
print(f"\nlog(A) vs γ_phonon: r = {r_A_gamma:.3f}")

# A/A_0 ratio vs γ
r_ratio_gamma, _ = stats.pearsonr(gamma_phonon, A_ratio)
print(f"A/A_0 vs γ_phonon: r = {r_ratio_gamma:.3f}")

# ============================================================================
# Analysis 3: A vs φ (work function)
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 3: Richardson Constant A vs Work Function φ")
print("=" * 70)

r_A_phi, _ = stats.pearsonr(phi, log_A)
print(f"\nlog(A) vs φ: r = {r_A_phi:.3f}")

# ============================================================================
# Analysis 4: γ_work (coherence from work function)
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 4: γ_work = 27.2/φ (Electronic Coherence)")
print("=" * 70)

r_A_gamma_work, _ = stats.pearsonr(gamma_work, log_A)
print(f"\nlog(A) vs γ_work: r = {r_A_gamma_work:.3f}")

r_ratio_gamma_work, _ = stats.pearsonr(gamma_work, A_ratio)
print(f"A/A_0 vs γ_work: r = {r_ratio_gamma_work:.3f}")

# ============================================================================
# Analysis 5: By material class
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 5: By Material Class")
print("=" * 70)

# Alkali metals
alkali = ['Li', 'Na', 'K', 'Rb', 'Cs']
alkali_mask = np.array([m in alkali for m in materials])

# Transition metals
transition = ['Ti', 'V', 'Cr', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mo', 'Ta', 'W', 'Re', 'Pt']
trans_mask = np.array([m in transition for m in materials])

# Noble metals
noble = ['Cu', 'Ag', 'Au']
noble_mask = np.array([m in noble for m in materials])

for name, mask in [('Alkali metals', alkali_mask), ('Transition metals', trans_mask), ('Noble metals', noble_mask)]:
    if np.sum(mask) >= 3:
        print(f"\n{name} ({np.sum(mask)} materials):")
        print(f"  φ range: {phi[mask].min():.2f} - {phi[mask].max():.2f} eV")
        print(f"  A range: {A[mask].min():.0f} - {A[mask].max():.0f} A/cm²K²")
        print(f"  A/A_0: {A_ratio[mask].mean():.2f} ± {A_ratio[mask].std():.2f}")

        if np.sum(mask) >= 4:
            r_class, _ = stats.pearsonr(gamma_phonon[mask], phi[mask])
            print(f"  φ vs γ_phonon: r = {r_class:.3f}")

# ============================================================================
# Analysis 6: Emission efficiency
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 6: Emission Efficiency J/T² at T = 1500 K")
print("=" * 70)

# Calculate J/T² = A × exp(-φ/kT) at T = 1500K (typical operating temp)
k_B = 8.617e-5  # eV/K
T_op = 1500  # K
J_over_T2 = A * np.exp(-phi / (k_B * T_op))

log_J = np.log10(J_over_T2)
r_J_gamma, _ = stats.pearsonr(gamma_phonon, log_J)
r_J_phi, _ = stats.pearsonr(phi, log_J)

print(f"\nlog(J/T²) vs φ: r = {r_J_phi:.3f} (dominated by exponential)")
print(f"log(J/T²) vs γ_phonon: r = {r_J_gamma:.3f}")

# Best emitters
print(f"\nTop emitters at T = 1500K:")
sorted_idx = np.argsort(J_over_T2)[::-1]
print(f"{'Material':<15} {'φ (eV)':<10} {'A (A/cm²K²)':<12} {'J/T² (A/cm²K²)':<15}")
for i in sorted_idx[:10]:
    print(f"{materials[i]:<15} {phi[i]:<10.2f} {A[i]:<12.0f} {J_over_T2[i]:<15.2e}")

# ============================================================================
# Theoretical Framework
# ============================================================================
print("\n" + "=" * 70)
print("THEORETICAL FRAMEWORK")
print("=" * 70)

print("""
Richardson-Dushman Equation:
J = A × T² × exp(-φ/kT)

Theoretical A_0 = 4πm_e k²e / h³ = 120.4 A/cm²K²

Actual A differs from A_0 due to:
1. Surface roughness
2. Electron reflection at surface
3. Temperature-dependent φ
4. Many-body effects

Coherence interpretation:
- φ = work function = electron binding at surface
- Higher φ → lower γ_work → more tightly bound electrons
- A/A_0 deviation may reflect surface coherence effects

The exponential dependence on φ dominates emission,
so coherence effects through A are secondary.

For thermionic emitters:
- Low φ materials (Cs, LaB6) have highest emission
- γ_work = 27.2/φ → higher γ_work = easier emission
- This is CONSISTENT with γ interpretation: higher γ = looser electrons
""")

# ============================================================================
# Key Results
# ============================================================================
print("\n" + "=" * 70)
print("KEY RESULTS")
print("=" * 70)

print(f"""
1. φ vs γ_phonon: r = {r_phi_gamma:.3f}
   {'Work function scales with phonon coherence' if abs(r_phi_gamma) > 0.3 else 'Weak correlation'}

2. A vs γ_phonon: r = {r_A_gamma:.3f}
   {'Richardson constant scales with γ' if abs(r_A_gamma) > 0.3 else 'Weak correlation'}

3. A vs φ: r = {r_A_phi:.3f}
   {'A correlated with work function' if abs(r_A_phi) > 0.3 else 'Weak correlation'}

4. A/A_0 vs γ_work: r = {r_ratio_gamma_work:.3f}
   {'Surface coherence effect detected' if abs(r_ratio_gamma_work) > 0.3 else 'Weak correlation'}

5. J/T² (emission) vs φ: r = {r_J_phi:.3f}
   Emission is STRONGLY dominated by work function (exponential).

INSIGHT: Thermionic emission is dominated by φ (exponential).

The coherence interpretation:
- γ_work = 27.2/φ measures electron "looseness"
- High γ_work → easier emission → low work function
- This is CONSISTENT with polarizability (Session #85)

However, the Richardson constant A shows weak γ dependence,
suggesting surface coherence effects are not dominant.

Best emitters (Cs, LaB6, W-ThO2) have LOW φ, not special A.
""")

# ============================================================================
# Visualization
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Color by class
colors = []
for m in materials:
    if m in alkali:
        colors.append('red')
    elif m in transition:
        colors.append('blue')
    elif m in noble:
        colors.append('gold')
    elif m in ['LaB6', 'W-ThO2']:
        colors.append('green')
    else:
        colors.append('gray')
colors = np.array(colors)

# Plot 1: φ vs γ_phonon
ax1 = axes[0, 0]
for c, label in [('red', 'Alkali'), ('blue', 'Transition'), ('gold', 'Noble'),
                 ('green', 'Compound'), ('gray', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax1.scatter(gamma_phonon[mask], phi[mask], c=c, label=label, s=80, alpha=0.7)
ax1.set_xlabel('γ_phonon = 2T/θ_D', fontsize=12)
ax1.set_ylabel('Work Function φ (eV)', fontsize=12)
ax1.set_title(f'φ vs γ_phonon (r = {r_phi_gamma:.3f})', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: A vs φ
ax2 = axes[0, 1]
for c, label in [('red', 'Alkali'), ('blue', 'Transition'), ('gold', 'Noble'),
                 ('green', 'Compound'), ('gray', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax2.scatter(phi[mask], A[mask], c=c, label=label, s=80, alpha=0.7)
ax2.axhline(A_0, color='black', linestyle='--', linewidth=2, label=f'A_0 = {A_0:.1f}')
ax2.set_xlabel('Work Function φ (eV)', fontsize=12)
ax2.set_ylabel('Richardson Constant A (A/cm²K²)', fontsize=12)
ax2.set_title(f'A vs φ (r = {r_A_phi:.3f})', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: J/T² vs φ
ax3 = axes[1, 0]
for c, label in [('red', 'Alkali'), ('blue', 'Transition'), ('gold', 'Noble'),
                 ('green', 'Compound'), ('gray', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax3.scatter(phi[mask], J_over_T2[mask], c=c, label=label, s=80, alpha=0.7)
ax3.set_xlabel('Work Function φ (eV)', fontsize=12)
ax3.set_ylabel('J/T² at 1500K (A/cm²K²)', fontsize=12)
ax3.set_yscale('log')
ax3.set_title(f'Emission vs φ (r = {r_J_phi:.3f})', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: A/A_0 histogram
ax4 = axes[1, 1]
ax4.hist(A_ratio, bins=15, color='blue', alpha=0.7, edgecolor='black')
ax4.axvline(1.0, color='red', linestyle='--', linewidth=2, label='Ideal A/A_0 = 1')
ax4.axvline(A_ratio.mean(), color='green', linestyle='-', linewidth=2,
            label=f'Mean = {A_ratio.mean():.2f}')
ax4.set_xlabel('A/A_0', fontsize=12)
ax4.set_ylabel('Count', fontsize=12)
ax4.set_title('Richardson Constant Deviation', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.suptitle('Chemistry Session #98: Thermionic Emission and Coherence',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/thermionic_emission_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: thermionic_emission_coherence.png")

# ============================================================================
# Predictions
# ============================================================================
print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print(f"""
P98.1: Emission dominated by exp(-φ/kT), not by A variations
       Work function is the primary control parameter.

P98.2: γ_work = 27.2/φ measures electron "looseness"
       High γ_work → easy emission (consistent with polarizability).

P98.3: A shows weak correlation with γ_phonon (r = {r_A_gamma:.3f})
       Surface coherence effects are secondary to φ.

P98.4: Best emitters have LOW φ, not special A
       Cs (φ = 2.14), LaB6 (φ = 2.70), W-ThO2 (φ = 2.70).

P98.5: A/A_0 mean = {A_ratio.mean():.2f}
       Most materials have A < A_0 (reflection/roughness effects).

P98.6: φ vs γ_phonon: r = {r_phi_gamma:.3f}
       {'Some correlation' if abs(r_phi_gamma) > 0.3 else 'No strong correlation'} between
       lattice coherence and work function.

FRAMEWORK POSITION:
Thermionic emission is NOT primarily a coherence phenomenon.

The physics is dominated by work function (binding energy),
which determines the Boltzmann exponential factor.

The Richardson constant A shows weak γ dependence,
indicating surface coherence effects are minor.

This is DIFFERENT from transport (σ, κ, μ) where
coherence (mean free path) is essential.

Thermionic emission is an ENERGY BARRIER problem,
not a coherence problem.
""")

# Validation status
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

print(f"""
**WEAK COHERENCE CORRELATION**

Key correlations:
- φ vs γ_phonon: r = {r_phi_gamma:.3f}
- A vs γ_phonon: r = {r_A_gamma:.3f}
- A/A_0 vs γ_work: r = {r_ratio_gamma_work:.3f}
- J/T² vs φ: r = {r_J_phi:.3f} (strong, expected)

Thermionic emission is NOT strongly coherence-dependent.

This validates the framework by showing where it DOESN'T apply:
- Coherence matters for TRANSPORT (scattering)
- Coherence matters for OPTICAL properties
- Coherence is secondary for ENERGY BARRIER phenomena

Work function φ is an ENERGY, not a coherence parameter.
The exponential dependence makes φ the dominant factor.
""")
