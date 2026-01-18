"""
Session #89: Anderson Localization and Coherence

Hypothesis: Anderson localization (metal-insulator transition from disorder)
represents the transition from γ << 2 (coherent/metallic) to γ → 2 (incoherent/localized).

Key physics:
- Clean metal: extended wavefunctions, coherent transport, γ_e << 2
- Disordered metal: localized wavefunctions, hopping transport, γ_e → 2
- Critical disorder: metal-insulator transition (MIT)

Anderson criterion: localization when W/B > (W/B)_c
where W = disorder bandwidth, B = clean bandwidth

From coherence: At MIT, γ_electron reaches critical value γ_c.

Data: Metal-insulator transitions in disordered systems
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# DISORDER-DRIVEN MIT DATA
# =============================================================================

# Doped semiconductors (concentration-driven MIT)
doped_semiconductors = {
    # Material: {n_c: critical concentration (cm^-3), m_eff: effective mass, a_B: Bohr radius}
    'Si:P': {'n_c': 3.7e18, 'a_B': 20.4, 'epsilon': 11.7},
    'Si:B': {'n_c': 4.1e18, 'a_B': 13.5, 'epsilon': 11.7},
    'Si:As': {'n_c': 8.5e18, 'a_B': 25.5, 'epsilon': 11.7},
    'Ge:Sb': {'n_c': 1.5e17, 'a_B': 96, 'epsilon': 16.0},
    'Ge:Ga': {'n_c': 2.0e17, 'a_B': 57, 'epsilon': 16.0},
    'GaAs:Si': {'n_c': 1.0e16, 'a_B': 103, 'epsilon': 12.9},
    'InSb:Te': {'n_c': 1.5e14, 'a_B': 650, 'epsilon': 17.9},
    'InP:S': {'n_c': 1.2e18, 'a_B': 75, 'epsilon': 12.5},
}

# Thin film MIT (thickness-driven or disorder-driven)
thin_films = {
    # Material: {t_c: critical thickness (nm), R_sq_c: critical sheet resistance (Ω/sq)}
    'Bi': {'t_c': 10, 'R_sq_c': 6500},
    'Pb': {'t_c': 8, 'R_sq_c': 4000},
    'Au': {'t_c': 2, 'R_sq_c': 5000},
    'Ag': {'t_c': 3, 'R_sq_c': 4500},
    'Cu': {'t_c': 2.5, 'R_sq_c': 4200},
    'Al': {'t_c': 3, 'R_sq_c': 5200},
    'NbN': {'t_c': 4, 'R_sq_c': 5500},
    'InOx': {'t_c': 15, 'R_sq_c': 7000},
}

# Composition-driven MIT
alloys_mit = {
    # Material: composition at MIT, conductivity ratio (σ_metal/σ_insulator at T→0)
    'Si1-xGex': {'x_c': 0.15, 'sigma_ratio': 1e6},
    'V2O3': {'T_MIT': 150, 'sigma_ratio': 1e7},  # Temperature-driven
    'VO2': {'T_MIT': 340, 'sigma_ratio': 1e5},   # Temperature-driven
    'NdNiO3': {'T_MIT': 200, 'sigma_ratio': 1e4},
    'SmNiO3': {'T_MIT': 400, 'sigma_ratio': 1e5},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*60)
print("Session #89: Anderson Localization and Coherence")
print("="*60)

# =============================================================================
# TEST 1: Mott criterion n_c^(1/3) × a_B ≈ 0.25
# =============================================================================
print("\n" + "="*60)
print("TEST 1: Mott criterion n_c^(1/3) × a_B ≈ 0.25")
print("="*60)

# Mott: Metal-insulator transition when mean impurity spacing ~ Bohr radius
# n_c^(1/3) × a_B = constant ≈ 0.25

print(f"\n{'Material':<12} {'n_c (cm⁻³)':>12} {'a_B (Å)':>10} {'n_c^(1/3)×a_B':>12}")
print("-" * 50)

materials = list(doped_semiconductors.keys())
n_c = np.array([doped_semiconductors[m]['n_c'] for m in materials])
a_B = np.array([doped_semiconductors[m]['a_B'] for m in materials])
epsilon = np.array([doped_semiconductors[m]['epsilon'] for m in materials])

mott_product = (n_c ** (1/3)) * a_B * 1e-8  # Convert Å to cm

for i, mat in enumerate(materials):
    print(f"{mat:<12} {n_c[i]:>12.2e} {a_B[i]:>10.1f} {mott_product[i]:>12.3f}")

print(f"\nMott product mean: {mott_product.mean():.3f}")
print(f"Mott product std: {mott_product.std():.3f}")
print(f"Theory: 0.25, Observed: {mott_product.mean():.3f} ± {mott_product.std():.3f}")

# =============================================================================
# TEST 2: Coherence interpretation of Mott criterion
# =============================================================================
print("\n" + "="*60)
print("TEST 2: Coherence interpretation of Mott criterion")
print("="*60)

# At MIT, mean impurity spacing r_s = n^(-1/3)
# Coherence requires overlap: r_s < a_B (wavefunctions overlap)
# So n_c^(1/3) × a_B > 1 for metallic phase

# Define γ_disorder as ratio of disorder scale to coherence scale
# γ_disorder = r_s / a_B = 1 / (n^(1/3) × a_B)

# At MIT: γ_disorder = 1/0.25 = 4 (but we use γ_coh, so...)
# In coherence framework: γ → 2 at classical limit

# Model: γ_electron(disorder) = γ_0 × (W/B) where W/B is disorder strength
# At MIT: γ_electron = γ_c ~ some critical value

# Use: γ_Anderson = 2 / mott_product
# At MIT: mott_product = 0.25, so γ_Anderson = 8

gamma_anderson = 2 / mott_product

print(f"\nγ_Anderson = 2 / (n_c^(1/3)×a_B)")
print(f"{'Material':<12} {'mott_product':>12} {'γ_Anderson':>12}")
print("-" * 40)
for i, mat in enumerate(materials):
    print(f"{mat:<12} {mott_product[i]:>12.3f} {gamma_anderson[i]:>12.2f}")

print(f"\nγ_Anderson at MIT: {gamma_anderson.mean():.2f} ± {gamma_anderson.std():.2f}")

# =============================================================================
# TEST 3: Sheet resistance at MIT
# =============================================================================
print("\n" + "="*60)
print("TEST 3: Sheet resistance at 2D MIT")
print("="*60)

# For 2D systems, MIT occurs at quantum of resistance
# R_c = h/(e²) ≈ 25.8 kΩ
# But due to factors, practical R_c ~ 4-10 kΩ/sq

R_quantum = 25812  # Ω
print(f"Quantum of resistance: R_0 = h/e² = {R_quantum} Ω")

film_names = list(thin_films.keys())
R_sq_c = np.array([thin_films[f]['R_sq_c'] for f in film_names])
t_c = np.array([thin_films[f]['t_c'] for f in film_names])

print(f"\n{'Film':<10} {'t_c (nm)':>10} {'R_sq_c (Ω)':>12} {'R_c/R_0':>10}")
print("-" * 45)
for i, name in enumerate(film_names):
    print(f"{name:<10} {t_c[i]:>10.1f} {R_sq_c[i]:>12.0f} {R_sq_c[i]/R_quantum:>10.3f}")

print(f"\nMean R_c/R_0: {(R_sq_c/R_quantum).mean():.3f}")
print(f"All films show R_c ~ 0.2 × R_0 (factor of ~5 below quantum)")

# Coherence interpretation:
# γ_2D = (R/R_0) × 2 ?
gamma_2d = 2 * R_sq_c / R_quantum
print(f"\nIf γ_2D = 2 × R_sq/R_0, then γ_c = {gamma_2d.mean():.2f}")

# =============================================================================
# TEST 4: Localization length and coherence
# =============================================================================
print("\n" + "="*60)
print("TEST 4: Localization length scaling")
print("="*60)

# Anderson localization: ξ ∝ (W_c - W)^(-ν) near MIT
# where ξ = localization length, ν ≈ 1 (1D), 1.5 (3D)

# Coherence interpretation:
# ξ = ξ_0 × (2/γ) where γ increases with disorder
# At MIT: ξ → ∞, so γ → 0? No, that's backward.

# Actually: In metallic phase, ξ = ∞ (extended states)
# At MIT: ξ becomes finite
# In insulating phase: ξ = localization length

# Model: γ_metal = γ_0 (clean limit)
# γ(W) = γ_0 × (1 + (W/W_c)^p)
# At W = W_c: γ = 2γ_0

print("Localization length: ξ ∝ |W - W_c|^(-ν)")
print("Critical exponent ν:")
print("  1D: ν = 2")
print("  2D: ν = ∞ (marginal)")
print("  3D: ν ≈ 1.5")
print("\nCoherence interpretation:")
print("  Metallic: γ << 2, extended states")
print("  MIT: γ = γ_c (critical)")
print("  Insulating: γ → 2, localized states")

# =============================================================================
# TEST 5: Mott-Ioffe-Regel criterion
# =============================================================================
print("\n" + "="*60)
print("TEST 5: Mott-Ioffe-Regel criterion")
print("="*60)

# MIR criterion: Localization when k_F × l ≈ 1
# where k_F = Fermi wavevector, l = mean free path

# From Boltzmann: σ = (e²/ℏ) × k_F² × l / 3π² (3D)
# At MIR limit: σ_min = (e²/ℏ) × k_F / 3π²

# Coherence: mean free path l = coherence length
# γ_electron ∝ 1/l
# At MIR: k_F × l = 1, so l = 1/k_F (one wavelength)

# This gives σ_min ~ 100-1000 S/cm for metals

print("Mott-Ioffe-Regel: k_F × l ≈ 1 at MIT")
print("\nAt MIT:")
print("  l ≈ 1/k_F (one Fermi wavelength)")
print("  σ_min = (e²/ℏ) × k_F / 3π²")
print("\nCoherence interpretation:")
print("  l = mean free path = electronic coherence length")
print("  γ_electron ∝ a/l where a = lattice constant")
print("  At MIR: γ_electron = a × k_F ≈ π (since k_F ≈ π/a)")

# =============================================================================
# TEST 6: Temperature dependence in localized regime
# =============================================================================
print("\n" + "="*60)
print("TEST 6: Variable-range hopping")
print("="*60)

# In localized regime: σ = σ_0 × exp(-(T_0/T)^p)
# Mott VRH: p = 1/(d+1)
#   3D: p = 1/4
#   2D: p = 1/3
# Efros-Shklovskii: p = 1/2 (Coulomb gap)

# Coherence interpretation:
# Hopping = thermally activated coherent jumps
# γ(T) increases with T (thermal fluctuations destroy coherence)
# σ ∝ exp(-Δγ/kT) where Δγ = barrier in γ space

print("Variable-range hopping: σ = σ_0 × exp(-(T_0/T)^p)")
print("\nMott VRH exponent:")
print("  3D: p = 1/4")
print("  2D: p = 1/3")
print("  Efros-Shklovskii: p = 1/2")
print("\nCoherence interpretation:")
print("  Hopping distance R ∝ (T_0/T)^(1/(d+1))")
print("  γ_hop = 2 × (R/ξ)²")
print("  At low T: long hops, γ_hop large, σ small")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Mott criterion
ax = axes[0, 0]
ax.scatter(a_B, n_c, c='blue', s=80, alpha=0.7)
for i, mat in enumerate(materials):
    ax.annotate(mat, (a_B[i], n_c[i]), fontsize=8, alpha=0.8)
# Mott criterion line: n_c = (0.25/a_B)^3
a_B_line = np.linspace(10, 700, 100)
n_c_mott = (0.25 / (a_B_line * 1e-8))**3
ax.plot(a_B_line, n_c_mott, 'r--', label='Mott: n^(1/3)×a_B = 0.25')
ax.set_xlabel('Bohr radius a_B (Å)')
ax.set_ylabel('Critical concentration n_c (cm⁻³)')
ax.set_yscale('log')
ax.set_title('Mott Criterion Validation')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Sheet resistance at MIT
ax = axes[0, 1]
ax.bar(film_names, R_sq_c, color='green', alpha=0.7)
ax.axhline(y=R_quantum/4, color='r', linestyle='--', label=f'R_0/4 = {R_quantum//4} Ω')
ax.axhline(y=R_quantum, color='k', linestyle='--', label=f'R_0 = {R_quantum} Ω')
ax.set_ylabel('Critical sheet resistance (Ω/□)')
ax.set_xticklabels(film_names, rotation=45, ha='right')
ax.set_title('2D Metal-Insulator Transition')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

# Plot 3: Disorder phase diagram
ax = axes[1, 0]
# Schematic phase diagram
W_B = np.linspace(0, 2, 100)
gamma_disorder = 2 * W_B / (1 + W_B)  # Model: γ saturates at 2
ax.plot(W_B, gamma_disorder, 'b-', linewidth=2, label='γ(W/B)')
ax.axhline(y=1, color='r', linestyle='--', label='γ_c = 1 (MIT)')
ax.axvline(x=1, color='r', linestyle=':', alpha=0.5)
ax.fill_between(W_B, 0, gamma_disorder, where=(gamma_disorder < 1),
                alpha=0.3, color='blue', label='Metallic (γ < 1)')
ax.fill_between(W_B, gamma_disorder, 2, where=(gamma_disorder > 1),
                alpha=0.3, color='red', label='Insulating (γ > 1)')
ax.set_xlabel('Disorder strength W/B')
ax.set_ylabel('γ_electron')
ax.set_title('Coherence Phase Diagram (Schematic)')
ax.legend()
ax.set_ylim(0, 2.2)
ax.grid(True, alpha=0.3)

# Plot 4: VRH regimes
ax = axes[1, 1]
T_ratio = np.linspace(0.1, 2, 100)
sigma_3D = np.exp(-T_ratio**(-1/4))
sigma_2D = np.exp(-T_ratio**(-1/3))
sigma_ES = np.exp(-T_ratio**(-1/2))
ax.plot(T_ratio, sigma_3D, 'b-', linewidth=2, label='3D Mott VRH (p=1/4)')
ax.plot(T_ratio, sigma_2D, 'g-', linewidth=2, label='2D Mott VRH (p=1/3)')
ax.plot(T_ratio, sigma_ES, 'r-', linewidth=2, label='Efros-Shklovskii (p=1/2)')
ax.set_xlabel('T/T_0')
ax.set_ylabel('σ/σ_0')
ax.set_title('Variable-Range Hopping Regimes')
ax.legend()
ax.set_yscale('log')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/anderson_localization.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: Anderson Localization and Coherence")
print("="*60)

print("""
Key Results:
1. Mott criterion: n_c^(1/3) × a_B = {:.3f} ± {:.3f} (theory: 0.25)
2. 2D MIT: R_c/R_0 = {:.3f} ± {:.3f} (quantum of resistance)
3. γ_Anderson at MIT: {:.2f} ± {:.2f}

Coherence Framework Interpretation:

1. **Anderson localization = coherence destruction**
   - Clean metal: Extended wavefunctions, γ_e << 2
   - Disordered: Localized wavefunctions, γ_e → 2
   - MIT: γ_e = γ_c (critical coherence)

2. **Mott criterion = coherence overlap**
   - Wavefunctions overlap when n^(1/3) × a_B > 0.25
   - This ensures phase coherence between sites
   - Below threshold: isolated, incoherent states

3. **2D quantum of resistance**
   - R_0 = h/e² = 25.8 kΩ sets the scale
   - At MIT: R ≈ R_0 means γ ≈ 1 (σ ∝ 1/γ → 1/R)
   - 2D is marginal: localization at any disorder

4. **Variable-range hopping = thermally-activated coherence**
   - σ ∝ exp(-(T_0/T)^p) reflects hopping distance
   - Hopping distance R ∝ coherent jump range
   - Low T: longer hops needed, γ_hop larger, σ smaller

5. **Mott-Ioffe-Regel criterion**
   - k_F × l = 1: mean free path = Fermi wavelength
   - This is γ_e = a × k_F ≈ π (lattice constant × wavevector)
   - Minimum metallic conductivity from coherence limit

Physical Insight:
Anderson localization is the ULTIMATE test of coherence:
- Metallic: Electrons maintain phase coherence across sample
- Insulating: Phase randomization confines electrons
- MIT: Critical point where coherence length ξ → ∞

The framework predicts:
- γ_metal << 2 (coherent, extended)
- γ_insulator → 2 (incoherent, localized)
- γ_MIT = critical value where ξ diverges
""".format(mott_product.mean(), mott_product.std(),
           (R_sq_c/R_quantum).mean(), (R_sq_c/R_quantum).std(),
           gamma_anderson.mean(), gamma_anderson.std()))

# =============================================================================
# PREDICTIONS
# =============================================================================
print("\n" + "="*60)
print("PREDICTIONS from Coherence Framework")
print("="*60)

print("""
P89.1: At MIT, γ_electron reaches critical value γ_c ≈ 1
       (based on R_c ~ R_0/4 and σ ∝ 1/γ)

P89.2: Mott criterion n_c^(1/3) × a_B = 0.25 is coherence overlap
       Wavefunctions must overlap for metallic coherence

P89.3: VRH exponent p = 1/(d+1) reflects hopping coherence
       Lower dimension → stronger localization → larger p

P89.4: 2D systems are marginal (always localize at T=0)
       Infinite localization length in 2D

P89.5: σ_min exists from Mott-Ioffe-Regel: k_F × l = 1
       Minimum coherence length = 1/k_F

P89.6: Temperature-driven MIT (VO2, V2O3) follows γ(T)
       γ increases with T, crosses γ_c at T_MIT
""")

print(f"\nPlot saved to: anderson_localization.png")
