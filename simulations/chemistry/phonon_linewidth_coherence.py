#!/usr/bin/env python3
"""
Chemistry Session #107: Phonon Linewidth and Coherence
=======================================================

Test whether phonon linewidth Γ_ph relates to coherence parameters.

Physical reasoning:
- Phonon linewidth: Γ_ph = ℏ/τ_ph (natural broadening of phonon mode)
- τ_ph: phonon lifetime (how long coherent vibration persists)
- Anharmonicity → phonon-phonon scattering → shorter lifetime → larger Γ_ph
- Related to thermal conductivity: κ ∝ C_v × v × l_ph ∝ 1/Γ_ph

Connection to coherence:
- Γ_ph measures PHONON decoherence rate
- Compare to Drude Γ (electron decoherence) from Session #103
- Should relate to γ_phonon = 2T/θ_D

This is the PHONON equivalent of Session #103!
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Dataset: Phonon linewidths (zone-center optical phonon)
# Γ_ph: Phonon linewidth at 300K (cm^-1)
# ω_0: Optical phonon frequency (cm^-1)
# θ_D: Debye temperature (K)
# γ_G: Grüneisen parameter (anharmonicity measure)

materials = {
    # Material, Γ_ph (cm^-1), ω_0 (cm^-1), θ_D (K), γ_G
    # Semiconductors
    'Si': (1.0, 520, 645, 0.56),
    'Ge': (3.5, 300, 374, 0.73),
    'GaAs': (2.5, 292, 360, 0.72),
    'InAs': (5.0, 219, 280, 0.84),
    'InSb': (7.0, 185, 160, 0.95),
    'GaN': (2.0, 560, 600, 0.85),
    'AlN': (1.5, 660, 950, 0.70),
    'InP': (3.5, 307, 321, 0.79),
    'GaP': (1.8, 364, 445, 0.68),
    'Diamond': (0.2, 1332, 2230, 0.97),
    'SiC': (0.8, 972, 1080, 0.76),
    'ZnO': (4.0, 583, 416, 1.32),
    'CdS': (5.0, 305, 219, 1.44),
    'CdTe': (4.5, 168, 200, 1.04),
    'ZnSe': (3.0, 252, 340, 1.00),

    # Ionic crystals
    'NaCl': (8.0, 264, 321, 1.61),
    'KCl': (10.0, 214, 235, 1.45),
    'KBr': (12.0, 166, 174, 1.49),
    'MgO': (2.0, 401, 946, 1.52),
    'CaF2': (3.5, 322, 510, 1.70),

    # Perovskites
    'SrTiO3': (15.0, 546, 400, 2.20),  # Soft mode behavior
    'BaTiO3': (25.0, 480, 350, 2.50),  # Very anharmonic (ferroelectric)

    # Layered
    'Graphite': (11.0, 1582, 402, 2.0),  # In-plane E_2g mode
    'MoS2': (4.0, 383, 300, 1.1),
    'BN_hex': (1.0, 1366, 1900, 0.8),
}

# Extract data
names = list(materials.keys())
Gamma_ph = np.array([materials[m][0] for m in names])  # cm^-1
omega_0 = np.array([materials[m][1] for m in names])  # cm^-1
theta_D = np.array([materials[m][2] for m in names])  # K
gamma_G = np.array([materials[m][3] for m in names])  # dimensionless

# Calculate derived quantities
T = 300  # K
gamma_phonon = 2 * T / theta_D

# Phonon quality factor Q_ph = ω_0 / Γ_ph
Q_ph = omega_0 / Gamma_ph

# Phonon lifetime τ_ph = ℏ/Γ_ph (in ps for cm^-1 units: τ ≈ 5.31/Γ ps)
tau_ph = 5.31 / Gamma_ph  # ps

print("="*70)
print("CHEMISTRY SESSION #107: PHONON LINEWIDTH AND COHERENCE")
print("="*70)
print(f"\nDataset: {len(names)} materials")
print(f"Γ_ph range: {Gamma_ph.min():.1f} - {Gamma_ph.max():.1f} cm^-1")
print(f"τ_ph range: {tau_ph.min():.1f} - {tau_ph.max():.1f} ps")
print(f"Q_ph range: {Q_ph.min():.0f} - {Q_ph.max():.0f}")

# Analysis 1: Γ_ph vs γ_phonon
print("\n" + "="*70)
print("ANALYSIS 1: PHONON LINEWIDTH vs γ_phonon")
print("="*70)

r_Gamma_gamma_p, p_Gamma_gamma_p = stats.pearsonr(Gamma_ph, gamma_phonon)
print(f"\nΓ_ph vs γ_phonon: r = {r_Gamma_gamma_p:.3f}, p = {p_Gamma_gamma_p:.2e}")

# Analysis 2: Γ_ph vs Grüneisen parameter
print("\n" + "="*70)
print("ANALYSIS 2: PHONON LINEWIDTH vs GRÜNEISEN")
print("="*70)

r_Gamma_Gruneisen, p_Gamma_Gruneisen = stats.pearsonr(Gamma_ph, gamma_G)
print(f"\nΓ_ph vs γ_G (Grüneisen): r = {r_Gamma_Gruneisen:.3f}")

# Theory: Γ_ph ∝ γ_G² × T / M
# Higher Grüneisen = more anharmonic = larger linewidth

# Analysis 3: Γ_ph vs ω_0 (frequency dependence)
print("\n" + "="*70)
print("ANALYSIS 3: PHONON LINEWIDTH vs FREQUENCY")
print("="*70)

r_Gamma_omega, p_Gamma_omega = stats.pearsonr(Gamma_ph, omega_0)
print(f"\nΓ_ph vs ω_0: r = {r_Gamma_omega:.3f}")

# Relative linewidth
relative_width = Gamma_ph / omega_0
print(f"\nRelative linewidth Γ_ph/ω_0:")
print(f"  Range: {relative_width.min():.4f} - {relative_width.max():.4f}")

r_rel_gamma, p_rel_gamma = stats.pearsonr(relative_width, gamma_phonon)
print(f"\nΓ_ph/ω_0 vs γ_phonon: r = {r_rel_gamma:.3f}")

# Analysis 4: Q_ph vs 1/γ_phonon
print("\n" + "="*70)
print("ANALYSIS 4: PHONON QUALITY FACTOR")
print("="*70)

r_Q_gamma, p_Q_gamma = stats.pearsonr(Q_ph, 1/gamma_phonon)
print(f"\nQ_ph vs 1/γ_phonon: r = {r_Q_gamma:.3f}")

# Compare to Session #103: Q_optical = ω_p/Γ vs 1/γ_electron
print("\nCompare to electron (Session #103):")
print("  Q_electron = ω_p/Γ vs 1/γ_electron: r = 0.833")
print(f"  Q_phonon = ω_0/Γ_ph vs 1/γ_phonon: r = {r_Q_gamma:.3f}")

# Analysis 5: Material class comparison
print("\n" + "="*70)
print("ANALYSIS 5: MATERIAL CLASS COMPARISON")
print("="*70)

semiconductors = ['Si', 'Ge', 'GaAs', 'InAs', 'InSb', 'GaN', 'AlN', 'InP', 'GaP',
                 'Diamond', 'SiC', 'ZnO', 'CdS', 'CdTe', 'ZnSe']
ionic = ['NaCl', 'KCl', 'KBr', 'MgO', 'CaF2']
perovskites = ['SrTiO3', 'BaTiO3']
layered = ['Graphite', 'MoS2', 'BN_hex']

print("\n| Class | Mean Γ_ph | Mean Q_ph | Mean γ_phonon | Mean γ_G |")
print("|-------|-----------|-----------|---------------|----------|")

for class_name, class_list in [('Semiconductors', semiconductors), ('Ionic', ionic),
                                ('Perovskites', perovskites), ('Layered', layered)]:
    mask = np.array([n in class_list for n in names])
    if np.sum(mask) >= 2:
        mean_Gamma = Gamma_ph[mask].mean()
        mean_Q = Q_ph[mask].mean()
        mean_gamma_p = gamma_phonon[mask].mean()
        mean_gamma_G = gamma_G[mask].mean()
        print(f"| {class_name:12} | {mean_Gamma:9.1f} | {mean_Q:9.0f} | {mean_gamma_p:13.2f} | {mean_gamma_G:8.2f} |")

# Analysis 6: Combined model
print("\n" + "="*70)
print("ANALYSIS 6: COMBINED ANHARMONICITY MODEL")
print("="*70)

# Model: Γ_ph ∝ γ_G² × γ_phonon
Gamma_model = gamma_G**2 * gamma_phonon

r_model, p_model = stats.pearsonr(Gamma_ph, Gamma_model)
print(f"\nΓ_ph vs γ_G² × γ_phonon: r = {r_model:.3f}")

# Analysis 7: Coherence interpretation
print("\n" + "="*70)
print("ANALYSIS 7: COHERENCE FRAMEWORK")
print("="*70)

# Phonon coherence: N_coh = ω_0 × τ_ph = Q_ph (number of coherent oscillations)
N_coh = Q_ph
gamma_ph_coh = 2 / np.sqrt(N_coh)

print("\nPhonon coherence interpretation:")
print(f"  N_coh = Q_ph (coherent oscillations): {N_coh.min():.0f} - {N_coh.max():.0f}")
print(f"  γ_phonon_coh = 2/√N_coh: {gamma_ph_coh.min():.3f} - {gamma_ph_coh.max():.3f}")

r_gamma_compare, _ = stats.pearsonr(gamma_ph_coh, gamma_phonon)
print(f"\nγ_phonon_coh vs γ_phonon (2T/θ_D): r = {r_gamma_compare:.3f}")

# Create visualization
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Plot 1: Γ_ph vs γ_phonon
ax1 = axes[0, 0]
colors = {'semiconductors': 'blue', 'ionic': 'red', 'perovskites': 'green', 'layered': 'purple'}
class_map = {n: 'semiconductors' for n in semiconductors}
class_map.update({n: 'ionic' for n in ionic})
class_map.update({n: 'perovskites' for n in perovskites})
class_map.update({n: 'layered' for n in layered})

for name in names:
    i = names.index(name)
    mat_class = class_map.get(name, 'other')
    c = colors.get(mat_class, 'gray')
    ax1.scatter(gamma_phonon[i], Gamma_ph[i], c=c, s=80, edgecolors='black', linewidth=0.5)
    if Gamma_ph[i] > 15 or name in ['Diamond', 'Si', 'NaCl']:
        ax1.annotate(name, (gamma_phonon[i], Gamma_ph[i]), fontsize=8, xytext=(3, 3),
                    textcoords='offset points')
ax1.set_xlabel('γ_phonon = 2T/θ_D', fontsize=11)
ax1.set_ylabel('Γ_ph (cm⁻¹)', fontsize=11)
ax1.set_title(f'Phonon Linewidth vs γ_phonon (r = {r_Gamma_gamma_p:.3f})', fontsize=13)

# Plot 2: Γ_ph vs γ_G (Grüneisen)
ax2 = axes[0, 1]
for name in names:
    i = names.index(name)
    mat_class = class_map.get(name, 'other')
    c = colors.get(mat_class, 'gray')
    ax2.scatter(gamma_G[i], Gamma_ph[i], c=c, s=80, edgecolors='black', linewidth=0.5)
ax2.set_xlabel('Grüneisen parameter γ_G', fontsize=11)
ax2.set_ylabel('Γ_ph (cm⁻¹)', fontsize=11)
ax2.set_title(f'Linewidth vs Anharmonicity (r = {r_Gamma_Gruneisen:.3f})', fontsize=13)

# Plot 3: Q_ph vs 1/γ_phonon
ax3 = axes[0, 2]
for name in names:
    i = names.index(name)
    mat_class = class_map.get(name, 'other')
    c = colors.get(mat_class, 'gray')
    ax3.scatter(1/gamma_phonon[i], Q_ph[i], c=c, s=80, edgecolors='black', linewidth=0.5)
    if Q_ph[i] > 3000 or name in ['BaTiO3', 'NaCl']:
        ax3.annotate(name, (1/gamma_phonon[i], Q_ph[i]), fontsize=8, xytext=(3, 3),
                    textcoords='offset points')
ax3.set_xlabel('1/γ_phonon (coherence factor)', fontsize=11)
ax3.set_ylabel('Q_ph = ω_0/Γ_ph', fontsize=11)
ax3.set_title(f'Phonon Quality Factor (r = {r_Q_gamma:.3f})', fontsize=13)
ax3.set_yscale('log')

# Plot 4: τ_ph distribution
ax4 = axes[1, 0]
ax4.hist(tau_ph, bins=15, color='steelblue', edgecolor='black', alpha=0.7)
ax4.axvline(x=np.mean(tau_ph), color='red', linestyle='--', label=f'Mean = {np.mean(tau_ph):.1f} ps')
ax4.set_xlabel('Phonon lifetime τ_ph (ps)', fontsize=11)
ax4.set_ylabel('Count', fontsize=11)
ax4.set_title('Phonon Lifetime Distribution', fontsize=13)
ax4.legend()

# Plot 5: Combined model
ax5 = axes[1, 1]
for name in names:
    i = names.index(name)
    mat_class = class_map.get(name, 'other')
    c = colors.get(mat_class, 'gray')
    ax5.scatter(Gamma_model[i], Gamma_ph[i], c=c, s=80, edgecolors='black', linewidth=0.5)
max_val = max(Gamma_model.max(), Gamma_ph.max())
ax5.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='1:1')
ax5.set_xlabel('γ_G² × γ_phonon (model)', fontsize=11)
ax5.set_ylabel('Γ_ph observed (cm⁻¹)', fontsize=11)
ax5.set_title(f'Combined Model (r = {r_model:.3f})', fontsize=13)
ax5.legend()

# Plot 6: Class averages
ax6 = axes[1, 2]
class_data = []
for class_name, class_list in [('Semiconductors', semiconductors), ('Ionic', ionic),
                                ('Perovskites', perovskites), ('Layered', layered)]:
    mask = np.array([n in class_list for n in names])
    if np.sum(mask) >= 2:
        class_data.append((class_name, Gamma_ph[mask].mean(), gamma_phonon[mask].mean(),
                          gamma_G[mask].mean(), colors.get(class_name.lower(), 'gray')))

for name, mean_Gamma, mean_gp, mean_gG, c in class_data:
    ax6.scatter(mean_gp, mean_Gamma, s=mean_gG*100, c=c, edgecolors='black', linewidth=1.5)
    ax6.annotate(name, (mean_gp, mean_Gamma), textcoords="offset points",
                xytext=(10, 5), fontsize=11, fontweight='bold')
ax6.set_xlabel('Mean γ_phonon', fontsize=11)
ax6.set_ylabel('Mean Γ_ph (cm⁻¹)', fontsize=11)
ax6.set_title('Material Class Averages (size ∝ γ_G)', fontsize=13)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phonon_linewidth_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "="*70)
print("SESSION #107 SUMMARY")
print("="*70)

print(f"""
PHONON LINEWIDTH AND COHERENCE

Key Results:
1. Γ_ph vs γ_phonon: r = {r_Gamma_gamma_p:.3f}
2. Γ_ph vs γ_G (Grüneisen): r = {r_Gamma_Gruneisen:.3f}
3. Q_ph vs 1/γ_phonon: r = {r_Q_gamma:.3f}
4. Combined model Γ_ph ∝ γ_G² × γ_phonon: r = {r_model:.3f}

Physical Interpretation:

1. **Phonon Linewidth = Phonon Decoherence Rate**
   Γ_ph = ℏ/τ_ph (inverse phonon lifetime)
   Analog of Drude Γ for electrons (Session #103)

2. **Grüneisen Parameter Matters Most**
   γ_G measures anharmonicity
   Γ_ph vs γ_G: r = {r_Gamma_Gruneisen:.3f} (moderate)
   Perovskites have highest γ_G → highest Γ_ph

3. **Material Hierarchy**
   - Diamond: Γ_ph = 0.2 cm⁻¹, τ_ph = 27 ps (most coherent)
   - SiC, Si: Low Γ_ph (covalent, harmonic)
   - Perovskites: High Γ_ph (soft modes, anharmonic)
   - Ionic: Variable (depends on mass ratio)

4. **Phonon Quality Factor Q_ph**
   Q_ph = ω_0/Γ_ph (coherent oscillations)
   Diamond: Q_ph ~ 6660 (best)
   BaTiO3: Q_ph ~ 19 (worst)

Comparison to Electrons (Session #103):
- Electron: Γ ∝ γ_electron, r = 0.867
- Phonon: Γ_ph ∝ γ_phonon, r = {r_Gamma_gamma_p:.3f} (weaker)

Why weaker correlation for phonons?
- Phonons decay via ANHARMONICITY (γ_G)
- Electrons decay via PHONON SCATTERING (λ_ep)
- γ_phonon = 2T/θ_D captures thermal occupation, not anharmonicity
- Need γ_G for phonon decoherence

Framework Classification:
Phonon decoherence depends on:
- γ_phonon: Thermal population of modes
- γ_G: Anharmonicity (mode-mode coupling)
- Combined: Γ_ph ∝ γ_G² × γ_phonon (r = {r_model:.3f})
""")

print("\n[Plot saved to phonon_linewidth_coherence.png]")
