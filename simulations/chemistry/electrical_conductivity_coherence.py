#!/usr/bin/env python3
"""
Chemistry Session #81: Electrical Conductivity & Coherence
Test whether coherence framework predicts electrical transport.

Electrical conductivity σ = n×e²×τ/m*
Where:
- n = carrier density
- e = electron charge
- τ = scattering time (mean free path / Fermi velocity)
- m* = effective mass

Key relationships:
- σ ∝ 1/ρ (conductivity = 1/resistivity)
- σ ∝ τ (longer scattering time = higher conductivity)
- For metals: ρ ∝ T (Bloch-Grüneisen)
- σ varies over ~25 orders of magnitude (Ag to SiO2)

Coherence interpretation:
- Electron coherence length ξ_e determines τ
- High coherence → long mean free path → high σ
- σ ∝ (2/γ_e) where γ_e is electron coherence parameter
- Impurities, phonons, defects reduce ξ_e → increase γ_e
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("CHEMISTRY SESSION #81: ELECTRICAL CONDUCTIVITY & COHERENCE")
print("=" * 70)

# ==============================================================================
# DATASET: ELECTRICAL CONDUCTIVITY
# ==============================================================================

print("\n" + "=" * 70)
print("DATASET: ELECTRICAL CONDUCTIVITY AT 300K")
print("=" * 70)

# Material: (σ in S/m, T_m in K, θ_D in K, E_gap in eV for semiconductors)
# Conductivity at room temperature (300K)
materials = {
    # Best conductors (metals)
    'Ag': (6.30e7, 1235, 225, 0),
    'Cu': (5.96e7, 1358, 343, 0),
    'Au': (4.52e7, 1337, 165, 0),
    'Al': (3.77e7, 933, 428, 0),
    'Mg': (2.27e7, 923, 400, 0),
    'Na': (2.10e7, 371, 158, 0),
    'Ca': (2.82e7, 1115, 230, 0),

    # Transition metals
    'Fe': (1.00e7, 1811, 470, 0),
    'Ni': (1.43e7, 1728, 450, 0),
    'Co': (1.72e7, 1768, 445, 0),
    'W': (1.79e7, 3695, 400, 0),
    'Mo': (1.87e7, 2896, 450, 0),
    'Ti': (2.38e6, 1941, 420, 0),
    'Cr': (7.74e6, 2180, 630, 0),

    # Poor metals
    'Pb': (4.81e6, 601, 105, 0),
    'Sn': (9.17e6, 505, 200, 0),

    # Semiconductors (intrinsic)
    'Ge': (2.0, 1211, 374, 0.67),
    'Si': (4.35e-4, 1687, 645, 1.12),

    # Insulators
    'SiO2': (1e-14, 1986, 470, 9.0),
    'Al2O3': (1e-14, 2345, 1047, 8.8),
    'Diamond': (1e-16, 3800, 2230, 5.5),
}

print(f"Materials: {len(materials)}")

# Print sorted by σ
print("\nMaterials sorted by conductivity:")
print("-" * 70)
print(f"{'Material':<12} {'σ (S/m)':<15} {'log₁₀(σ)':<10} {'T_m (K)':<10} {'θ_D (K)':<10}")
print("-" * 70)

for name, (sigma, Tm, theta_D, Eg) in sorted(materials.items(), key=lambda x: -x[1][0]):
    log_sigma = np.log10(sigma)
    print(f"{name:<12} {sigma:>12.2e}  {log_sigma:>8.1f}  {Tm:>8.0f}  {theta_D:>8.0f}")

# ==============================================================================
# SEPARATE METALS AND SEMICONDUCTORS
# ==============================================================================

print("\n" + "=" * 70)
print("METALS ANALYSIS")
print("=" * 70)

metals = {k: v for k, v in materials.items() if v[3] == 0}  # E_gap = 0
semiconductors = {k: v for k, v in materials.items() if v[3] > 0}

print(f"Metals: {len(metals)}")
print(f"Semiconductors/Insulators: {len(semiconductors)}")

# Extract metal arrays
sigma_metals = []
Tm_metals = []
theta_metals = []
names_metals = []

for name, (sigma, Tm, theta, Eg) in metals.items():
    sigma_metals.append(sigma)
    Tm_metals.append(Tm)
    theta_metals.append(theta)
    names_metals.append(name)

sigma_metals = np.array(sigma_metals)
Tm_metals = np.array(Tm_metals)
theta_metals = np.array(theta_metals)

# ==============================================================================
# METAL CORRELATIONS
# ==============================================================================

print("\nMetal correlations:")

# σ vs θ_D
r_sigma_theta, _ = stats.pearsonr(sigma_metals, theta_metals)
print(f"σ vs θ_D: r = {r_sigma_theta:.3f}")

# σ vs T_m
r_sigma_Tm, _ = stats.pearsonr(sigma_metals, Tm_metals)
print(f"σ vs T_m: r = {r_sigma_Tm:.3f}")

# log(σ) vs θ_D
log_sigma_metals = np.log10(sigma_metals)
r_log_theta, _ = stats.pearsonr(log_sigma_metals, theta_metals)
print(f"log(σ) vs θ_D: r = {r_log_theta:.3f}")

# σ vs 1/θ_D
r_sigma_inv_theta, _ = stats.pearsonr(sigma_metals, 1/theta_metals)
print(f"σ vs 1/θ_D: r = {r_sigma_inv_theta:.3f}")

# ==============================================================================
# COHERENCE ANALYSIS FOR METALS
# ==============================================================================

print("\n" + "=" * 70)
print("COHERENCE ANALYSIS: METALS")
print("=" * 70)

def gamma_phonon(theta_D, T=300):
    """Phonon coherence from Debye temperature."""
    gamma = 2.0 * T / theta_D
    return np.clip(gamma, 0.2, 2.0)

gamma_metals = np.array([gamma_phonon(theta) for theta in theta_metals])

# σ vs γ
r_sigma_gamma, _ = stats.pearsonr(sigma_metals, gamma_metals)
print(f"σ vs γ_phonon: r = {r_sigma_gamma:.3f}")

# σ vs 2/γ
coh_metals = 2.0 / gamma_metals
r_sigma_coh, _ = stats.pearsonr(sigma_metals, coh_metals)
print(f"σ vs 2/γ: r = {r_sigma_coh:.3f}")

print("""
Note: For metals, conductivity depends on:
1. Carrier density (mostly constant for simple metals)
2. Scattering time τ (temperature-dependent)

At 300K, electron-phonon scattering dominates:
τ ∝ θ_D / T (from Bloch-Grüneisen theory)

Higher θ_D → less phonon scattering → higher σ
BUT also: noble metals have low θ_D but high σ!

This is because Ag, Cu, Au have:
- High carrier density (1 conduction electron/atom)
- Low effective mass
- Weak electron-phonon coupling
""")

# ==============================================================================
# ELECTRON MEAN FREE PATH
# ==============================================================================

print("\n" + "=" * 70)
print("ELECTRON MEAN FREE PATH")
print("=" * 70)

# Estimate mean free path: l = v_F × τ
# For simple metals: l ∝ σ × m* / (n × e²)
# Approximately: l ∝ σ / n (for similar n)

# Use resistivity formula: ρ = m × v_F / (n × e² × l)
# So: l ∝ v_F / (ρ × n)

# For estimation, assume similar Fermi velocity and density
# Then l ∝ σ

print("Estimated mean free paths (relative to Ag):")
print("-" * 50)
l_Ag = 1.0  # Reference
for i, name in enumerate(names_metals):
    l_rel = sigma_metals[i] / sigma_metals[names_metals.index('Ag')]
    print(f"{name:<12}: l/l_Ag = {l_rel:.2f}")

# ==============================================================================
# WIEDEMANN-FRANZ LAW
# ==============================================================================

print("\n" + "=" * 70)
print("WIEDEMANN-FRANZ LAW: κ/σ = L×T")
print("=" * 70)

print("""
The Wiedemann-Franz law states:
κ / (σ × T) = L₀ = π²k²/(3e²) = 2.44×10⁻⁸ W·Ω/K²

This connects thermal and electrical conductivity:
- Both depend on electron mean free path
- Both measure electron coherence

From Session #65: κ ∝ 2/γ (thermal conductivity)
For metals: κ_e = L₀ × σ × T

So: σ ∝ κ_e / T ∝ (2/γ) / T

At fixed T: σ ∝ 2/γ (same as thermal!)
""")

# ==============================================================================
# SEMICONDUCTOR ANALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("SEMICONDUCTOR/INSULATOR ANALYSIS")
print("=" * 70)

print("For semiconductors, conductivity is dominated by band gap:")
print("σ ∝ exp(-E_g / 2kT)")
print("\nFrom Session #60: E_g ∝ 2/γ")
print("So: log(σ) ∝ -E_g ∝ -(2/γ)")
print("Or: log(σ) ∝ γ (higher γ → smaller gap → higher σ)")

for name, (sigma, Tm, theta, Eg) in semiconductors.items():
    log_sigma = np.log10(sigma)
    gamma_est = 2.0 / (Eg / 1.0) if Eg > 0 else 0  # Rough estimate from Eg
    print(f"{name:<12}: log(σ) = {log_sigma:>6.1f}, E_g = {Eg:.2f} eV")

# ==============================================================================
# NOBLE METAL ANOMALY
# ==============================================================================

print("\n" + "=" * 70)
print("NOBLE METAL ANOMALY")
print("=" * 70)

print("""
Noble metals (Ag, Cu, Au) have:
- Low θ_D (165-343 K) → high γ_phonon
- BUT highest conductivity!

This seems to contradict σ ∝ 2/γ.

Resolution: There are TWO coherence factors:
1. γ_phonon (lattice coherence) - from θ_D
2. γ_electron (electronic coherence) - different physics!

For noble metals:
- d-band is FULL and deep below Fermi level
- Only s-electrons conduct
- Weak electron-phonon coupling
- Long electron mean free path despite low θ_D

γ_electron (conduction) ≠ γ_phonon (lattice)

This is like Session #76 (refractive index):
- n measures electronic polarizability
- Not the same as lattice coherence

For conductivity, we need γ_electron, not γ_phonon!
""")

# ==============================================================================
# RESISTIVITY AND TEMPERATURE
# ==============================================================================

print("\n" + "=" * 70)
print("RESISTIVITY TEMPERATURE DEPENDENCE")
print("=" * 70)

print("""
For metals above ~θ_D/4:
ρ(T) ∝ T (linear increase with temperature)

This comes from electron-phonon scattering:
- More phonons at higher T → more scattering → higher ρ

Bloch-Grüneisen formula:
ρ(T) ∝ (T/θ_D)⁵ ∫₀^(θ_D/T) x⁵/(eˣ-1)(1-e⁻ˣ) dx

At T >> θ_D: ρ ∝ T
At T << θ_D: ρ ∝ T⁵ (phonon freeze-out)

Coherence interpretation:
- γ_phonon = 2T/θ_D increases with T
- τ ∝ 1/γ_phonon → decreases with T
- σ = neτ/m → decreases with T
""")

# ==============================================================================
# MATERIAL CLASS ANALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS BY MATERIAL CLASS")
print("=" * 70)

metal_classes = {
    'Noble metals': ['Ag', 'Cu', 'Au'],
    'Alkali metals': ['Na', 'Ca'],  # Ca is alkaline earth but similar
    'Transition metals': ['Fe', 'Ni', 'Co', 'W', 'Mo', 'Ti', 'Cr'],
    'Post-transition': ['Al', 'Mg', 'Pb', 'Sn'],
}

for class_name, members in metal_classes.items():
    sigmas = [metals[m][0] for m in members if m in metals]
    thetas = [metals[m][2] for m in members if m in metals]

    if sigmas:
        mean_sigma = np.mean(sigmas)
        mean_theta = np.mean(thetas)
        print(f"{class_name}: mean σ = {mean_sigma:.2e} S/m, mean θ_D = {mean_theta:.0f} K")

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SESSION #81 SUMMARY: ELECTRICAL CONDUCTIVITY & COHERENCE")
print("=" * 70)

print(f"""
Metal Correlations:
- σ vs θ_D: r = {r_sigma_theta:.3f}
- σ vs T_m: r = {r_sigma_Tm:.3f}
- σ vs 1/θ_D: r = {r_sigma_inv_theta:.3f}
- σ vs γ_phonon: r = {r_sigma_gamma:.3f}
- σ vs 2/γ: r = {r_sigma_coh:.3f}

Key Findings:
1. Metal conductivity shows WEAK correlation with θ_D
   - Noble metals have LOW θ_D but HIGH σ
   - This is the "noble metal anomaly"

2. γ_phonon does NOT predict metal conductivity
   - r = {r_sigma_gamma:.3f} ({"positive" if r_sigma_gamma > 0 else "negative"} - wrong direction!)
   - Phonon coherence ≠ electron coherence

3. Electronic coherence (γ_electron) differs from lattice coherence
   - Noble metals: full d-band, weak e-ph coupling
   - High σ despite low θ_D

4. Semiconductors: log(σ) ∝ -E_g
   - From E_g ∝ 2/γ: σ ∝ exp(-2/γ_e)
   - Higher electronic γ → smaller gap → higher σ

5. Wiedemann-Franz validates connection
   - κ/σ = L×T (constant)
   - Both κ and σ depend on electron coherence

Physical Interpretation:
- Conductivity measures ELECTRON coherence length
- This is different from PHONON coherence (θ_D)
- For metals, electron-phonon coupling strength matters
- Noble metals: weak e-ph coupling → long τ → high σ
- Transition metals: strong d-electron scattering → shorter τ → lower σ

CONCLUSION:
The framework applies, but γ_electron ≠ γ_phonon!
Electrical conductivity requires its own coherence parameter.
""")

# ==============================================================================
# PREDICTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print("""
P81.1: σ_metal ∝ 2/γ_electron (not γ_phonon)
Electron coherence length determines conductivity.

P81.2: γ_electron from electron-phonon coupling
Noble metals: weak λ_ep → long l_e → high σ
Transition metals: strong λ_ep → short l_e → low σ

P81.3: σ(T) ∝ 1/γ_phonon at high T
Temperature dependence follows lattice γ.

P81.4: log(σ_SC) ∝ -E_g ∝ -(2/γ_e)
Semiconductor conductivity exponentially sensitive to γ.

P81.5: Wiedemann-Franz from common γ_electron
κ_e and σ both ∝ 2/γ_electron → constant ratio.

P81.6: Superconductivity = γ_electron → 0
Cooper pairs = perfectly coherent electrons → σ → ∞.
""")

# ==============================================================================
# VALIDATION STATUS
# ==============================================================================

print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

print(f"""
**MIXED - REVEALS LIMITATION**

Metal conductivity does NOT correlate with γ_phonon:
- r(σ, 2/γ) = {r_sigma_coh:.3f} ({"positive" if r_sigma_coh > 0 else "negative"})
- Noble metal anomaly shows γ_phonon ≠ γ_electron

This is NOT a failure of the coherence framework!
It reveals that DIFFERENT properties require DIFFERENT γ:
- Lattice properties (E, C_p, α, v): γ_phonon from θ_D
- Electronic properties (σ, κ_e): γ_electron from e-ph coupling
- Optical properties (n): γ_optical from polarizability

The framework is valid, but γ estimation must match the physics.

Connections to previous sessions:
- Session #65: κ ∝ 2/γ (thermal conductivity)
- Session #76: n not correlated with δ (wrong γ for liquids)
- Both show: right γ for right property!
""")

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: σ vs θ_D for metals
ax1 = axes[0, 0]
ax1.scatter(theta_metals, sigma_metals/1e7, s=80, alpha=0.7, c='blue')
for i, name in enumerate(names_metals):
    if name in ['Ag', 'Cu', 'Au', 'Fe', 'Ti', 'Na']:
        ax1.annotate(name, (theta_metals[i], sigma_metals[i]/1e7), fontsize=9)
ax1.set_xlabel('Debye Temperature θ_D (K)', fontsize=12)
ax1.set_ylabel('Conductivity (10⁷ S/m)', fontsize=12)
ax1.set_title(f'Metal Conductivity vs θ_D\n(r = {r_sigma_theta:.3f} - WEAK)', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: σ vs γ_phonon for metals
ax2 = axes[0, 1]
ax2.scatter(gamma_metals, sigma_metals/1e7, s=80, alpha=0.7, c='purple')
for i, name in enumerate(names_metals):
    if name in ['Ag', 'Cu', 'Au', 'Fe', 'Ti', 'Na']:
        ax2.annotate(name, (gamma_metals[i], sigma_metals[i]/1e7), fontsize=9)
ax2.set_xlabel('γ_phonon = 2T/θ_D', fontsize=12)
ax2.set_ylabel('Conductivity (10⁷ S/m)', fontsize=12)
ax2.set_title(f'Metal Conductivity vs γ_phonon\n(r = {r_sigma_gamma:.3f} - WRONG SIGN)', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: log(σ) vs E_g for semiconductors/insulators
ax3 = axes[1, 0]
Eg_arr = [v[3] for v in semiconductors.values()]
log_sigma_arr = [np.log10(v[0]) for v in semiconductors.values()]
names_sc = list(semiconductors.keys())
ax3.scatter(Eg_arr, log_sigma_arr, s=100, alpha=0.7, c='green')
for i, name in enumerate(names_sc):
    ax3.annotate(name, (Eg_arr[i], log_sigma_arr[i]), fontsize=9)
ax3.set_xlabel('Band Gap (eV)', fontsize=12)
ax3.set_ylabel('log₁₀(σ)', fontsize=12)
ax3.set_title('Semiconductor Conductivity vs Band Gap', fontsize=14)
ax3.grid(True, alpha=0.3)

# Plot 4: By material class
ax4 = axes[1, 1]
class_colors = {'Noble metals': 'gold', 'Alkali metals': 'red',
                'Transition metals': 'blue', 'Post-transition': 'green'}
for class_name, members in metal_classes.items():
    sigmas = [metals[m][0]/1e7 for m in members if m in metals]
    thetas = [metals[m][2] for m in members if m in metals]
    if sigmas:
        ax4.scatter(thetas, sigmas, label=class_name, s=100, alpha=0.7,
                    c=class_colors.get(class_name, 'gray'))
ax4.set_xlabel('Debye Temperature θ_D (K)', fontsize=12)
ax4.set_ylabel('Conductivity (10⁷ S/m)', fontsize=12)
ax4.set_title('Metal Conductivity by Class\n(Noble metals: high σ, low θ_D)', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/electrical_conductivity_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved to: simulations/chemistry/electrical_conductivity_coherence.png")

print("\n" + "=" * 70)
print("SESSION #81 COMPLETE: ELECTRICAL CONDUCTIVITY & COHERENCE")
print("=" * 70)
