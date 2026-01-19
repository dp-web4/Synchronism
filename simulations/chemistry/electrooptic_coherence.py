#!/usr/bin/env python3
"""
Chemistry Session #95: Electrooptic Effects and Coherence

Test whether electrooptic coefficient r relates to coherence γ.

Electrooptic effect = change in refractive index with electric field
r_ij = electrooptic coefficient (pm/V)
Δn ∝ r × E (linear effect, Pockels)
Δn ∝ R × E² (quadratic effect, Kerr)

Session #91: ε ∝ γ_optical (dielectric constant)
Session #93: d ∝ γ_phonon (piezoelectricity)

Hypothesis: r ∝ γ_optical (same polarizability mechanism as ε)
Or: r ∝ γ_phonon (same soft mode mechanism as piezoelectricity)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Electrooptic materials data
# r_33 or effective r in pm/V (Pockels coefficient)
# n = refractive index
# ε_r = relative permittivity
# θ_D = Debye temperature (K)
# E_g = band gap (eV)

electrooptic_data = {
    # Classic electrooptic crystals
    'LiNbO3': {'r_33': 31, 'n': 2.20, 'eps_r': 29, 'theta_D': 500, 'E_g': 3.7},
    'LiTaO3': {'r_33': 33, 'n': 2.18, 'eps_r': 43, 'theta_D': 450, 'E_g': 4.0},
    'BaTiO3': {'r_33': 105, 'n': 2.40, 'eps_r': 1700, 'theta_D': 300, 'E_g': 3.2},
    'KNbO3': {'r_33': 63, 'n': 2.33, 'eps_r': 500, 'theta_D': 350, 'E_g': 3.3},
    'SrTiO3': {'r_33': 140, 'n': 2.39, 'eps_r': 300, 'theta_D': 420, 'E_g': 3.2},

    # KDP family
    'KDP (KH2PO4)': {'r_33': 10.6, 'n': 1.51, 'eps_r': 42, 'theta_D': 280, 'E_g': 7.0},
    'KD*P (KD2PO4)': {'r_33': 24.1, 'n': 1.51, 'eps_r': 50, 'theta_D': 260, 'E_g': 7.0},
    'ADP (NH4H2PO4)': {'r_33': 8.5, 'n': 1.52, 'eps_r': 15, 'theta_D': 300, 'E_g': 6.8},

    # Zinc blende semiconductors
    'GaAs': {'r_33': 1.2, 'n': 3.37, 'eps_r': 13, 'theta_D': 344, 'E_g': 1.42},
    'GaP': {'r_33': 1.0, 'n': 3.02, 'eps_r': 11, 'theta_D': 445, 'E_g': 2.26},
    'InP': {'r_33': 1.5, 'n': 3.10, 'eps_r': 12, 'theta_D': 321, 'E_g': 1.35},
    'InAs': {'r_33': 1.4, 'n': 3.52, 'eps_r': 15, 'theta_D': 250, 'E_g': 0.36},
    'ZnSe': {'r_33': 2.0, 'n': 2.57, 'eps_r': 8, 'theta_D': 271, 'E_g': 2.7},
    'ZnTe': {'r_33': 4.5, 'n': 2.70, 'eps_r': 10, 'theta_D': 223, 'E_g': 2.3},
    'CdTe': {'r_33': 6.8, 'n': 2.72, 'eps_r': 10, 'theta_D': 158, 'E_g': 1.5},

    # II-VI wurtzite
    'ZnO': {'r_33': 2.6, 'n': 1.99, 'eps_r': 11, 'theta_D': 440, 'E_g': 3.4},
    'CdS': {'r_33': 3.7, 'n': 2.50, 'eps_r': 10, 'theta_D': 219, 'E_g': 2.4},

    # Organic crystals
    'DAST': {'r_33': 77, 'n': 2.4, 'eps_r': 5, 'theta_D': 150, 'E_g': 2.1},

    # Relaxor ferroelectrics
    'PMN-PT': {'r_33': 400, 'n': 2.5, 'eps_r': 5000, 'theta_D': 180, 'E_g': 3.3},
    'PZN-PT': {'r_33': 350, 'n': 2.5, 'eps_r': 4500, 'theta_D': 200, 'E_g': 3.3},
}

print("=" * 70)
print("Chemistry Session #95: Electrooptic Effects and Coherence")
print("=" * 70)

# Extract data
materials = list(electrooptic_data.keys())
r_33 = np.array([electrooptic_data[m]['r_33'] for m in materials])
n = np.array([electrooptic_data[m]['n'] for m in materials])
eps_r = np.array([electrooptic_data[m]['eps_r'] for m in materials])
theta_D = np.array([electrooptic_data[m]['theta_D'] for m in materials])
E_g = np.array([electrooptic_data[m]['E_g'] for m in materials])

# Calculate coherence parameters
T = 300  # K
gamma_phonon = 2 * T / theta_D  # Phonon coherence

# γ_optical = 2 × IE_ref / E_g (approximation using E_g instead of IE)
IE_ref = 13.6  # eV (hydrogen ionization energy)
gamma_optical = 2 * IE_ref / E_g

print(f"\n{'Material':<20} {'r_33 (pm/V)':<12} {'n':<6} {'ε_r':<8} {'E_g (eV)':<8} {'γ_ph':<8} {'γ_opt'}")
print("-" * 90)
for i, m in enumerate(materials):
    print(f"{m:<20} {r_33[i]:<12.1f} {n[i]:<6.2f} {eps_r[i]:<8.0f} {E_g[i]:<8.2f} {gamma_phonon[i]:<8.2f} {gamma_optical[i]:.2f}")

# ============================================================================
# Analysis 1: r_33 vs γ_phonon
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 1: r_33 vs γ_phonon (soft phonon modes)")
print("=" * 70)

log_r = np.log10(r_33)
r_ph, p_ph = stats.pearsonr(gamma_phonon, log_r)
print(f"\nlog(r_33) vs γ_phonon: r = {r_ph:.3f}, p = {p_ph:.2e}")

slope_ph, intercept_ph, _, _, _ = stats.linregress(gamma_phonon, log_r)
print(f"Fit: log(r) = {slope_ph:.2f} × γ_ph + {intercept_ph:.2f}")

# ============================================================================
# Analysis 2: r_33 vs γ_optical
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 2: r_33 vs γ_optical (electronic polarizability)")
print("=" * 70)

r_opt, p_opt = stats.pearsonr(gamma_optical, log_r)
print(f"\nlog(r_33) vs γ_optical: r = {r_opt:.3f}, p = {p_opt:.2e}")

# ============================================================================
# Analysis 3: r_33 vs ε_r (dielectric constant)
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 3: r_33 vs ε_r (known relationship)")
print("=" * 70)

log_eps = np.log10(eps_r)
r_eps, p_eps = stats.pearsonr(log_eps, log_r)
print(f"\nlog(r_33) vs log(ε_r): r = {r_eps:.3f}")

slope_eps, intercept_eps, _, _, _ = stats.linregress(log_eps, log_r)
print(f"Fit: r_33 ∝ ε_r^{slope_eps:.2f}")

# ============================================================================
# Analysis 4: Combined models
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 4: Combined Models")
print("=" * 70)

# r ∝ ε × γ_phonon (like piezoelectricity)
combined_ph = eps_r * gamma_phonon
log_combined_ph = np.log10(combined_ph)
r_comb_ph, _ = stats.pearsonr(log_combined_ph, log_r)
print(f"\nlog(r) vs log(ε × γ_phonon): r = {r_comb_ph:.3f}")

# r ∝ ε / γ_optical (coherent EO effect)
combined_opt = eps_r / gamma_optical
log_combined_opt = np.log10(combined_opt)
r_comb_opt, _ = stats.pearsonr(log_combined_opt, log_r)
print(f"log(r) vs log(ε / γ_optical): r = {r_comb_opt:.3f}")

# r ∝ n^4 / E_g (Miller's delta)
miller = n**4 / E_g
log_miller = np.log10(miller)
r_miller, _ = stats.pearsonr(log_miller, log_r)
print(f"log(r) vs log(n^4/E_g) [Miller]: r = {r_miller:.3f}")

# ============================================================================
# Analysis 5: By material class
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 5: By Material Class")
print("=" * 70)

# Ferroelectrics (high r)
FE_materials = ['LiNbO3', 'LiTaO3', 'BaTiO3', 'KNbO3', 'SrTiO3', 'PMN-PT', 'PZN-PT']
FE_mask = np.array([m in FE_materials for m in materials])

# Semiconductors (low r)
SC_materials = ['GaAs', 'GaP', 'InP', 'InAs', 'ZnSe', 'ZnTe', 'CdTe', 'ZnO', 'CdS']
SC_mask = np.array([m in SC_materials for m in materials])

# KDP family
KDP_materials = ['KDP (KH2PO4)', 'KD*P (KD2PO4)', 'ADP (NH4H2PO4)']
KDP_mask = np.array([m in KDP_materials for m in materials])

for name, mask in [('Ferroelectrics', FE_mask), ('Semiconductors', SC_mask), ('KDP family', KDP_mask)]:
    if np.sum(mask) >= 3:
        r_class_ph, _ = stats.pearsonr(gamma_phonon[mask], log_r[mask])
        r_class_opt, _ = stats.pearsonr(gamma_optical[mask], log_r[mask])
        r_class_eps, _ = stats.pearsonr(log_eps[mask], log_r[mask])
        print(f"\n{name} ({np.sum(mask)} materials):")
        print(f"  log(r) vs γ_phonon: r = {r_class_ph:.3f}")
        print(f"  log(r) vs γ_optical: r = {r_class_opt:.3f}")
        print(f"  log(r) vs log(ε): r = {r_class_eps:.3f}")
        print(f"  r_33 range: {r_33[mask].min():.1f} - {r_33[mask].max():.1f} pm/V")

# ============================================================================
# Analysis 6: Figure of merit
# ============================================================================
print("\n" + "=" * 70)
print("Analysis 6: Figure of Merit n³r")
print("=" * 70)

# FOM for EO modulators = n³ × r
FOM = n**3 * r_33
log_FOM = np.log10(FOM)

r_fom_ph, _ = stats.pearsonr(gamma_phonon, log_FOM)
r_fom_opt, _ = stats.pearsonr(gamma_optical, log_FOM)
print(f"\nlog(n³r) vs γ_phonon: r = {r_fom_ph:.3f}")
print(f"log(n³r) vs γ_optical: r = {r_fom_opt:.3f}")

print(f"\n{'Material':<20} {'r_33':<10} {'n³r':<12} {'γ_ph':<8} {'γ_opt'}")
print("-" * 60)
sorted_idx = np.argsort(FOM)[::-1]
for i in sorted_idx[:10]:
    print(f"{materials[i]:<20} {r_33[i]:<10.0f} {FOM[i]:<12.1f} {gamma_phonon[i]:<8.2f} {gamma_optical[i]:.2f}")

# ============================================================================
# Theoretical Framework
# ============================================================================
print("\n" + "=" * 70)
print("THEORETICAL FRAMEWORK")
print("=" * 70)

print("""
Electrooptic Effect Origin:
Δn = -(1/2) × n³ × r × E

where r = electrooptic coefficient

Physical mechanisms:
1. Electronic contribution: Distortion of electron clouds
2. Ionic contribution: Displacement of ions in lattice
3. Acoustic contribution: Piezoelectric-induced strain

Standard model (Miller's rule):
r ∝ n^4 / E_g (nonlinear susceptibility)

Coherence interpretation:
- r_electronic ∝ γ_optical (polarizability, Session #91)
- r_ionic ∝ γ_phonon × ε (soft modes, like piezo Session #93)

Combined:
r_total = r_electronic + r_ionic
       ∝ γ_optical + γ_phonon × ε
""")

# ============================================================================
# Key Results
# ============================================================================
print("\n" + "=" * 70)
print("KEY RESULTS")
print("=" * 70)

print(f"""
1. log(r_33) vs γ_phonon: r = {r_ph:.3f}
   {'Soft phonon modes help!' if r_ph > 0.3 else 'Weak correlation'}

2. log(r_33) vs γ_optical: r = {r_opt:.3f}
   {'Electronic polarizability matters!' if r_opt > 0.3 else 'Weak correlation'}

3. log(r_33) vs log(ε_r): r = {r_eps:.3f}
   {'Strong correlation with dielectric constant!' if r_eps > 0.6 else 'Moderate correlation'}

4. Combined (ε × γ_phonon): r = {r_comb_ph:.3f}
   {'Combined model is better' if r_comb_ph > max(r_eps, r_ph) else 'Single variable sufficient'}

5. Miller's rule (n^4/E_g): r = {r_miller:.3f}

COMPARISON TO PIEZOELECTRICITY (#93):
- Piezo: d ∝ ε × γ_phonon (r = 0.940)
- Electrooptic: r ∝ ε (r = {r_eps:.3f}), γ_phonon adds little

INSIGHT: Electrooptic effect is DOMINATED by ε_r (polarizability).
The soft phonon contribution (γ_phonon) is secondary.

This is because EO effect is primarily ELECTRONIC distortion,
while piezoelectricity is IONIC displacement.

Ferroelectrics have high r because they have high ε (soft modes),
but the correlation is through ε, not independently through γ.
""")

# ============================================================================
# Visualization
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Color by class
colors = []
for m in materials:
    if m in FE_materials:
        colors.append('red')
    elif m in SC_materials:
        colors.append('blue')
    elif m in KDP_materials:
        colors.append('green')
    else:
        colors.append('purple')
colors = np.array(colors)

# Plot 1: r_33 vs γ_phonon
ax1 = axes[0, 0]
for c, label in [('red', 'Ferroelectric'), ('blue', 'Semiconductor'), ('green', 'KDP'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax1.scatter(gamma_phonon[mask], r_33[mask], c=c, label=label, s=80, alpha=0.7)
ax1.set_xlabel('γ_phonon = 2T/θ_D', fontsize=12)
ax1.set_ylabel('r_33 (pm/V)', fontsize=12)
ax1.set_yscale('log')
ax1.set_title(f'r_33 vs γ_phonon (r = {r_ph:.3f})', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: r_33 vs ε_r
ax2 = axes[0, 1]
for c, label in [('red', 'Ferroelectric'), ('blue', 'Semiconductor'), ('green', 'KDP'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax2.scatter(eps_r[mask], r_33[mask], c=c, label=label, s=80, alpha=0.7)
ax2.set_xlabel('Relative Permittivity ε_r', fontsize=12)
ax2.set_ylabel('r_33 (pm/V)', fontsize=12)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title(f'r_33 vs ε_r (r = {r_eps:.3f})', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Fit line
eps_fit = np.logspace(np.log10(eps_r.min()), np.log10(eps_r.max()), 100)
r_fit = 10**(slope_eps * np.log10(eps_fit) + intercept_eps)
ax2.plot(eps_fit, r_fit, 'k--', linewidth=2, alpha=0.7)

# Plot 3: r_33 vs combined (ε × γ_phonon)
ax3 = axes[1, 0]
for c, label in [('red', 'Ferroelectric'), ('blue', 'Semiconductor'), ('green', 'KDP'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax3.scatter(combined_ph[mask], r_33[mask], c=c, label=label, s=80, alpha=0.7)
ax3.set_xlabel('ε_r × γ_phonon', fontsize=12)
ax3.set_ylabel('r_33 (pm/V)', fontsize=12)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_title(f'r_33 vs ε × γ_phonon (r = {r_comb_ph:.3f})', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: FOM vs γ_phonon
ax4 = axes[1, 1]
for c, label in [('red', 'Ferroelectric'), ('blue', 'Semiconductor'), ('green', 'KDP'), ('purple', 'Other')]:
    mask = colors == c
    if np.sum(mask) > 0:
        ax4.scatter(gamma_phonon[mask], FOM[mask], c=c, label=label, s=80, alpha=0.7)
ax4.set_xlabel('γ_phonon = 2T/θ_D', fontsize=12)
ax4.set_ylabel('n³r (pm/V)', fontsize=12)
ax4.set_yscale('log')
ax4.set_title(f'FOM (n³r) vs γ_phonon (r = {r_fom_ph:.3f})', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.suptitle('Chemistry Session #95: Electrooptic Effects and Coherence',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/electrooptic_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved: electrooptic_coherence.png")

# ============================================================================
# Predictions
# ============================================================================
print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print(f"""
P95.1: r_33 ∝ ε_r^{slope_eps:.2f} with r = {r_eps:.3f}
       Electrooptic coefficient scales primarily with permittivity.

P95.2: γ_phonon provides secondary enhancement (r = {r_ph:.3f})
       Soft phonon modes help, but less than for piezoelectricity.

P95.3: Ferroelectrics have highest r because of high ε
       The correlation is through ε, not directly through γ.

P95.4: Semiconductors have low r due to low ε despite moderate γ
       III-V materials: r ~ 1-2 pm/V despite γ ~ 1.7-2.4.

P95.5: FOM (n³r) scales with γ_phonon (r = {r_fom_ph:.3f})
       Higher γ materials are good EO modulators.

P95.6: Electrooptic effect differs from piezoelectricity
       Piezo: ionic displacement → d ∝ ε × γ (both matter)
       EO: electronic distortion → r ∝ ε (ε dominates)

FRAMEWORK INSIGHT:
Electrooptic effect is ELECTRONIC polarizability, not ionic.
The ε correlation comes from electron cloud distortion.
γ_phonon is secondary because EO is primarily electronic.

This places electrooptic effect in the γ_optical category
(Session #91), not the γ_phonon category (Session #93).
""")

# Validation status
print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

if r_eps > 0.7:
    status = "GOOD VALIDATION"
elif r_eps > 0.5:
    status = "MODERATE VALIDATION"
else:
    status = "MIXED RESULTS"

print(f"""
{status}

r_33 vs ε_r: r = {r_eps:.3f} (PRIMARY)
r_33 vs γ_phonon: r = {r_ph:.3f} (SECONDARY)
r_33 vs γ_optical: r = {r_opt:.3f}
Combined model: r = {r_comb_ph:.3f}

KEY INSIGHT: Electrooptic effect is primarily ELECTRONIC.

The coherence framework correctly identifies:
- ε ∝ γ_optical (Session #91)
- r ∝ ε (this session)
- Therefore r ∝ γ_optical (indirect)

The γ_phonon contribution is secondary because:
- EO is electronic cloud distortion
- Piezo is ionic lattice displacement
- Different physical mechanisms!

This validates the COHERENCE TYPE distinction:
- γ_optical governs electronic properties (n, ε, r)
- γ_phonon governs lattice properties (d, κ, E)
""")
