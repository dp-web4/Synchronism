"""
Session #91: Dielectric Constant and Optical Coherence

Hypothesis: Static dielectric constant ε relates to coherence through
the Clausius-Mossotti relation and γ_optical.

Key relationships:
- Clausius-Mossotti: (ε-1)/(ε+2) = (4π/3) × n × α
- Where n = number density, α = polarizability
- From Session #85: α ∝ γ_optical^3.4

For optical frequency:
- Refractive index n² = ε (at optical freq)
- Session #76: n ∝ γ^(1/4) via Moss's rule

Coherence prediction:
- ε_static should correlate with γ_optical
- Higher polarizability (higher γ_optical) → higher ε
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# DIELECTRIC CONSTANT DATA
# Sources: CRC Handbook, various literature
# =============================================================================

# ε_r = relative permittivity (static dielectric constant)
materials = {
    # Elemental semiconductors
    'Diamond': {'eps_r': 5.7, 'n': 2.42, 'E_g': 5.47, 'theta_D': 1860},
    'Si': {'eps_r': 11.9, 'n': 3.44, 'E_g': 1.12, 'theta_D': 645},
    'Ge': {'eps_r': 16.0, 'n': 4.00, 'E_g': 0.67, 'theta_D': 374},

    # III-V compounds
    'GaAs': {'eps_r': 12.9, 'n': 3.30, 'E_g': 1.42, 'theta_D': 344},
    'GaP': {'eps_r': 11.1, 'n': 2.90, 'E_g': 2.26, 'theta_D': 445},
    'GaSb': {'eps_r': 15.7, 'n': 3.82, 'E_g': 0.73, 'theta_D': 265},
    'InAs': {'eps_r': 14.6, 'n': 3.52, 'E_g': 0.36, 'theta_D': 249},
    'InP': {'eps_r': 12.5, 'n': 3.10, 'E_g': 1.35, 'theta_D': 321},
    'InSb': {'eps_r': 17.9, 'n': 3.96, 'E_g': 0.17, 'theta_D': 202},
    'AlAs': {'eps_r': 10.1, 'n': 2.87, 'E_g': 2.16, 'theta_D': 417},
    'AlSb': {'eps_r': 14.4, 'n': 3.30, 'E_g': 1.58, 'theta_D': 292},

    # II-VI compounds
    'ZnS': {'eps_r': 8.9, 'n': 2.36, 'E_g': 3.68, 'theta_D': 340},
    'ZnSe': {'eps_r': 9.2, 'n': 2.57, 'E_g': 2.70, 'theta_D': 271},
    'ZnTe': {'eps_r': 10.1, 'n': 2.72, 'E_g': 2.26, 'theta_D': 223},
    'CdS': {'eps_r': 8.7, 'n': 2.50, 'E_g': 2.42, 'theta_D': 215},
    'CdSe': {'eps_r': 10.2, 'n': 2.55, 'E_g': 1.74, 'theta_D': 181},
    'CdTe': {'eps_r': 10.4, 'n': 2.67, 'E_g': 1.49, 'theta_D': 158},

    # Lead chalcogenides
    'PbS': {'eps_r': 17.0, 'n': 4.10, 'E_g': 0.41, 'theta_D': 227},
    'PbSe': {'eps_r': 23.6, 'n': 4.80, 'E_g': 0.27, 'theta_D': 130},
    'PbTe': {'eps_r': 32.8, 'n': 5.60, 'E_g': 0.31, 'theta_D': 130},

    # Ionic crystals
    'NaCl': {'eps_r': 5.9, 'n': 1.54, 'E_g': 8.5, 'theta_D': 321},
    'KCl': {'eps_r': 4.8, 'n': 1.49, 'E_g': 8.4, 'theta_D': 235},
    'KBr': {'eps_r': 4.9, 'n': 1.56, 'E_g': 7.6, 'theta_D': 174},
    'MgO': {'eps_r': 9.6, 'n': 1.74, 'E_g': 7.8, 'theta_D': 946},

    # Oxides
    'SiO2': {'eps_r': 3.9, 'n': 1.46, 'E_g': 9.0, 'theta_D': 470},
    'Al2O3': {'eps_r': 9.3, 'n': 1.76, 'E_g': 8.8, 'theta_D': 1047},
    'TiO2': {'eps_r': 86, 'n': 2.61, 'E_g': 3.2, 'theta_D': 760},
    'SrTiO3': {'eps_r': 300, 'n': 2.39, 'E_g': 3.2, 'theta_D': 513},
    'BaTiO3': {'eps_r': 1200, 'n': 2.40, 'E_g': 3.2, 'theta_D': 490},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*60)
print("Session #91: Dielectric Constant and Optical Coherence")
print("="*60)

# Extract data
names = list(materials.keys())
eps_r = np.array([materials[m]['eps_r'] for m in names])
n_opt = np.array([materials[m]['n'] for m in names])
E_g = np.array([materials[m]['E_g'] for m in names])
theta_D = np.array([materials[m]['theta_D'] for m in names])

print(f"\nDataset: {len(names)} materials")
print(f"ε_r range: {eps_r.min():.1f} - {eps_r.max():.0f}")
print(f"n range: {n_opt.min():.2f} - {n_opt.max():.2f}")
print(f"E_g range: {E_g.min():.2f} - {E_g.max():.1f} eV")

# =============================================================================
# TEST 1: ε vs n² (optical-static relation)
# =============================================================================
print("\n" + "="*60)
print("TEST 1: ε_static vs n² (optical dielectric)")
print("="*60)

# At optical frequency: ε_optical = n²
eps_optical = n_opt ** 2

r_eps_opt, p = stats.pearsonr(eps_optical, eps_r)
print(f"ε_static vs ε_optical (n²): r = {r_eps_opt:.3f}, p = {p:.4f}")

# The difference is due to ionic contribution at low frequency
eps_ionic = eps_r - eps_optical
print(f"\nIonic contribution (ε_static - n²):")
print(f"  Mean: {eps_ionic.mean():.1f}")
print(f"  Range: {eps_ionic.min():.1f} to {eps_ionic.max():.0f}")

# =============================================================================
# TEST 2: ε vs band gap (Penn model)
# =============================================================================
print("\n" + "="*60)
print("TEST 2: Penn model ε ∝ 1 + (ℏω_p/E_g)²")
print("="*60)

# Penn model: ε ≈ 1 + (E_p/E_g)² where E_p ~ 15-20 eV for most semiconductors
# Simplified: ε ∝ 1/E_g² (for E_p constant)

r_eps_Eg, _ = stats.pearsonr(1/E_g**2, eps_r)
print(f"ε vs 1/E_g²: r = {r_eps_Eg:.3f}")

# Log-log
r_log, _ = stats.pearsonr(np.log(E_g), np.log(eps_r))
slope, intercept, r_pow, _, _ = stats.linregress(np.log(E_g), np.log(eps_r))
print(f"log(ε) vs log(E_g): r = {r_log:.3f}, slope = {slope:.2f}")
print(f"This means ε ∝ E_g^{slope:.2f}")

# =============================================================================
# TEST 3: γ_optical from E_g
# =============================================================================
print("\n" + "="*60)
print("TEST 3: Coherence interpretation via γ_optical")
print("="*60)

# From Session #85: γ_optical = 2 × E_ref / E
# Here E can be E_g (band gap as optical coherence scale)

E_ref = 13.6  # eV (hydrogen)
gamma_optical = 2 * E_ref / E_g
gamma_optical = np.clip(gamma_optical, 0.5, 50)  # Reasonable range

print(f"γ_optical = 2 × 13.6 / E_g")
print(f"γ_optical range: {gamma_optical.min():.2f} - {gamma_optical.max():.1f}")

r_eps_gamma, _ = stats.pearsonr(gamma_optical, eps_r)
print(f"ε vs γ_optical: r = {r_eps_gamma:.3f}")

# ε ∝ γ^n ?
slope_g, _, r_pow_g, _, _ = stats.linregress(np.log(gamma_optical), np.log(eps_r))
print(f"ε ∝ γ_optical^{slope_g:.2f} (r = {r_pow_g:.3f})")

# =============================================================================
# TEST 4: Clausius-Mossotti analysis
# =============================================================================
print("\n" + "="*60)
print("TEST 4: Clausius-Mossotti function")
print("="*60)

# CM = (ε-1)/(ε+2) ∝ N×α (polarizability per unit volume)
CM = (eps_r - 1) / (eps_r + 2)

# From Session #85: α ∝ γ^3.4
# So CM should correlate with γ^3.4

r_CM_gamma, _ = stats.pearsonr(gamma_optical**2, CM)
print(f"CM function vs γ²: r = {r_CM_gamma:.3f}")

# =============================================================================
# TEST 5: Separate by material class
# =============================================================================
print("\n" + "="*60)
print("TEST 5: Analysis by material class")
print("="*60)

classes = {
    'Semiconductors': ['Diamond', 'Si', 'Ge', 'GaAs', 'GaP', 'GaSb', 'InAs', 'InP', 'InSb'],
    'II-VI': ['ZnS', 'ZnSe', 'ZnTe', 'CdS', 'CdSe', 'CdTe'],
    'Lead salts': ['PbS', 'PbSe', 'PbTe'],
    'Ionic': ['NaCl', 'KCl', 'KBr', 'MgO', 'SiO2', 'Al2O3'],
    'Perovskites': ['TiO2', 'SrTiO3', 'BaTiO3'],
}

print(f"\n{'Class':<15} {'n':>4} {'r(ε,1/E_g²)':>12} {'r(ε,γ)':>10}")
print("-" * 45)

for cls_name, mats in classes.items():
    idx = [names.index(m) for m in mats if m in names]
    if len(idx) < 3:
        print(f"{cls_name:<15} {len(idx):>4} (too few)")
        continue

    eps_cls = eps_r[idx]
    Eg_cls = E_g[idx]
    gamma_cls = gamma_optical[idx]

    r1, _ = stats.pearsonr(1/Eg_cls**2, eps_cls)
    r2, _ = stats.pearsonr(gamma_cls, eps_cls)
    print(f"{cls_name:<15} {len(idx):>4} {r1:>12.3f} {r2:>10.3f}")

# =============================================================================
# TEST 6: Ferroelectric anomaly (BaTiO3, SrTiO3)
# =============================================================================
print("\n" + "="*60)
print("TEST 6: Ferroelectric anomaly")
print("="*60)

print("\nPerovskite dielectrics:")
print(f"{'Material':<10} {'ε_r':>8} {'E_g':>6} {'γ_optical':>10} {'Type':<15}")
print("-" * 55)
perovskites = ['TiO2', 'SrTiO3', 'BaTiO3']
for mat in perovskites:
    i = names.index(mat)
    mat_type = 'Ferroelectric' if mat == 'BaTiO3' else 'Paraelectric'
    print(f"{mat:<10} {eps_r[i]:>8.0f} {E_g[i]:>6.1f} {gamma_optical[i]:>10.2f} {mat_type:<15}")

print("""
Ferroelectric anomaly:
- BaTiO3: ε ~ 1200-3000 (temperature dependent)
- SrTiO3: ε ~ 300 (quantum paraelectric)
- TiO2: ε ~ 86 (normal dielectric)

These cannot be explained by simple γ_optical model!
Ferroelectricity involves COLLECTIVE coherence:
- Soft phonon modes (structural instability)
- Domain structure (coherent domains)
- Temperature-dependent phase transition
""")

# =============================================================================
# TEST 7: Excluding ferroelectrics
# =============================================================================
print("\n" + "="*60)
print("TEST 7: Excluding ferroelectric anomalies")
print("="*60)

# Exclude extreme ferroelectrics
normal_mask = eps_r < 50  # Exclude TiO2, SrTiO3, BaTiO3
names_normal = [names[i] for i in range(len(names)) if normal_mask[i]]
eps_normal = eps_r[normal_mask]
Eg_normal = E_g[normal_mask]
gamma_normal = gamma_optical[normal_mask]

r_normal, _ = stats.pearsonr(1/Eg_normal**2, eps_normal)
print(f"Normal dielectrics (ε < 50): {len(eps_normal)} materials")
print(f"ε vs 1/E_g²: r = {r_normal:.3f}")

r_gamma_normal, _ = stats.pearsonr(gamma_normal, eps_normal)
print(f"ε vs γ_optical: r = {r_gamma_normal:.3f}")

# Power law fit
slope_n, _, r_pow_n, _, _ = stats.linregress(np.log(Eg_normal), np.log(eps_normal))
print(f"ε ∝ E_g^{slope_n:.2f} (r = {r_pow_n:.3f})")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(14, 9))

# Plot 1: ε vs n²
ax = axes[0, 0]
ax.scatter(eps_optical, eps_r, c='blue', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (eps_optical[i], eps_r[i]), fontsize=6, alpha=0.8)
max_val = max(eps_optical.max(), eps_r.max())
ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='ε_s = n²')
ax.set_xlabel('ε_optical = n²')
ax.set_ylabel('ε_static')
ax.set_title(f'Static vs Optical ε (r = {r_eps_opt:.3f})')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: ε vs E_g
ax = axes[0, 1]
ax.scatter(E_g, eps_r, c='green', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (E_g[i], eps_r[i]), fontsize=6, alpha=0.8)
ax.set_xlabel('E_g (eV)')
ax.set_ylabel('ε_static')
ax.set_yscale('log')
ax.set_title(f'Dielectric vs Band Gap (r = {r_log:.3f})')
ax.grid(True, alpha=0.3)

# Plot 3: ε vs γ_optical (excluding ferroelectrics)
ax = axes[0, 2]
ax.scatter(gamma_normal, eps_normal, c='red', s=60, alpha=0.7)
for i, name in enumerate(names_normal):
    idx = names.index(name)
    ax.annotate(name, (gamma_normal[i], eps_normal[i]), fontsize=6, alpha=0.8)
ax.set_xlabel('γ_optical = 2×13.6/E_g')
ax.set_ylabel('ε_static')
ax.set_title(f'Dielectric vs γ_optical (r = {r_gamma_normal:.3f})')
ax.grid(True, alpha=0.3)

# Plot 4: Log-log plot
ax = axes[1, 0]
ax.scatter(E_g[normal_mask], eps_r[normal_mask], c='purple', s=60, alpha=0.7)
for i, name in enumerate(names_normal):
    idx = names.index(name)
    ax.annotate(name, (Eg_normal[i], eps_normal[i]), fontsize=6, alpha=0.8)
# Fit line
x_fit = np.linspace(E_g.min(), E_g.max(), 100)
y_fit = np.exp(intercept) * x_fit**slope
ax.plot(x_fit, y_fit, 'r--', label=f'ε ∝ E_g^{slope:.2f}')
ax.set_xlabel('E_g (eV)')
ax.set_ylabel('ε_static')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f'Log-log: ε ∝ E_g^{slope_n:.2f}')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 5: Clausius-Mossotti
ax = axes[1, 1]
ax.scatter(gamma_optical**2, CM, c='orange', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_optical[i]**2, CM[i]), fontsize=6, alpha=0.8)
ax.set_xlabel('γ_optical²')
ax.set_ylabel('(ε-1)/(ε+2) [Clausius-Mossotti]')
ax.set_title(f'CM function vs γ² (r = {r_CM_gamma:.3f})')
ax.grid(True, alpha=0.3)

# Plot 6: Material classes
ax = axes[1, 2]
class_colors = {'Semiconductors': 'blue', 'II-VI': 'green', 'Lead salts': 'red',
                'Ionic': 'purple', 'Perovskites': 'orange'}
for cls_name, mats in classes.items():
    idx = [names.index(m) for m in mats if m in names]
    if len(idx) > 0:
        ax.scatter(E_g[idx], eps_r[idx], c=class_colors.get(cls_name, 'gray'),
                   s=80, alpha=0.7, label=cls_name)
ax.set_xlabel('E_g (eV)')
ax.set_ylabel('ε_static')
ax.set_yscale('log')
ax.set_title('Dielectric by Material Class')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/dielectric_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: Dielectric Constant and Coherence")
print("="*60)

print("""
Key Results:
1. ε_static vs n²: r = {:.3f} (optical contribution)
2. ε vs 1/E_g²: r = {:.3f} (Penn model)
3. ε ∝ E_g^{:.2f}: Power law exponent
4. ε vs γ_optical (normal): r = {:.3f}

Coherence Framework Interpretation:

1. **Penn model from coherence**
   - ε = 1 + (ω_p/E_g)² where ω_p = plasma frequency
   - E_g sets the optical coherence scale
   - γ_optical = 2×E_ref/E_g
   - So ε ∝ γ_optical² (approximately)

2. **Why ε ∝ 1/E_g²?**
   - Small E_g → easy electronic polarization
   - Small E_g → large γ_optical → loose electrons
   - Easy to polarize = high dielectric constant

3. **Ferroelectric anomaly**
   - BaTiO3, SrTiO3: ε >> Penn prediction
   - These involve COLLECTIVE coherence
   - Soft phonon modes create coherent dipole domains
   - γ_collective < γ_single-electron
   - Can't be explained by simple γ_optical

4. **Ionic contribution**
   - ε_static > n² due to ionic polarization
   - Ionic ε ∝ (Ω_TO)^(-2) from oscillator model
   - Lower Ω_TO (soft phonons) → higher ionic ε

5. **Material class variations**
   - Lead salts (PbS, PbSe, PbTe): Anomalously high ε
   - Due to soft phonon modes and high-Z atoms
   - Ionic crystals: ε dominated by ionic contribution
""".format(r_eps_opt, r_eps_Eg, slope, r_gamma_normal))

print("\nValidation summary:")
print(f"  Normal dielectrics: ε ∝ γ_optical with r = {r_gamma_normal:.3f}")
print(f"  Penn model: ε ∝ 1/E_g² with r = {r_normal:.3f}")
print(f"  Power law: ε ∝ E_g^{slope_n:.2f}")
print(f"\n  Ferroelectrics require collective coherence model")

print(f"\nPlot saved to: dielectric_coherence.png")
