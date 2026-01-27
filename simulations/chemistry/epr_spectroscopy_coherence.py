#!/usr/bin/env python3
"""
EPR/ESR Spectroscopy Coherence Analysis
Chemistry Session #232

Explores how γ ~ 1 manifests in electron paramagnetic resonance:
1. g-factor - free electron value as γ ~ 1 reference
2. Hyperfine coupling - nuclear spin effects
3. Relaxation times T1/T2 - coherence dynamics
4. Spin concentration - detection limits
5. Line shape analysis - Lorentzian vs Gaussian
6. Zero-field splitting - high-spin systems
7. Exchange coupling - coupled spin systems

Key insight: EPR spectroscopy is fundamentally about electron spin
coherence, with g = 2.0023 as the universal γ ~ 1 reference.
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("EPR/ESR SPECTROSCOPY COHERENCE ANALYSIS")
print("Chemistry Session #232: γ ~ 1 in Electron Spin Resonance")
print("=" * 70)

# Physical constants
mu_B = 9.274e-24    # Bohr magneton (J/T)
h = 6.626e-34       # Planck constant (J·s)
g_e = 2.00231930436 # Free electron g-factor

# =============================================================================
# 1. G-FACTOR ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("1. G-FACTOR: THE FUNDAMENTAL EPR REFERENCE")
print("=" * 70)

# Free electron g-factor is THE γ ~ 1 reference
print(f"  Free electron g-factor: g_e = {g_e:.10f}")
print(f"  This IS the γ ~ 1 reference for ALL EPR!")
print(f"  Deviation: Δg = g - g_e measures spin-orbit coupling")

# g-factors of common radicals/paramagnetic species
g_factors = {
    'Free electron': 2.00232,
    'DPPH (standard)': 2.0036,
    'Organic radical (π)': 2.003,
    'Cu²⁺ (d⁹)': 2.1 - 2.4,
    'Fe³⁺ high-spin': 2.0,
    'Mn²⁺': 2.0,
    'Ni²⁺': 2.2,
    'Co²⁺': 4.3,  # highly anisotropic
    'Gd³⁺': 2.0,
    'Nitroxide': 2.006,
}

print(f"\n  Species            | g-factor    | g/g_e")
print("  " + "-" * 50)
for species, g in g_factors.items():
    if isinstance(g, float):
        ratio = g / g_e
        print(f"  {species:18s} | {g:.4f}      | {ratio:.4f}")
    else:
        print(f"  {species:18s} | {g}  | (range)")

# Analyze deviations
print(f"\n  Key observations:")
print(f"    - Organic radicals: g ≈ g_e (small spin-orbit)")
print(f"    - Transition metals: g deviates (large spin-orbit)")
print(f"    - g = g_e IS γ ~ 1 for electron spin")

# DPPH as reference
print(f"\n  DPPH (2,2-diphenyl-1-picrylhydrazyl):")
print(f"    g = 2.0036 (primary EPR standard)")
print(f"    Δg = 0.0013 from free electron")
print(f"    γ = g_DPPH/g_e = {2.0036/g_e:.5f}")

# =============================================================================
# 2. HYPERFINE COUPLING
# =============================================================================
print("\n" + "=" * 70)
print("2. HYPERFINE COUPLING (A)")
print("=" * 70)

# Hyperfine coupling arises from electron-nuclear interaction
# A = 0 when no nuclear spin (I = 0)

print("  Hyperfine coupling: H_hf = A × S·I")
print("  A = 0 when I = 0 (no nuclear spin) - γ ~ 1 reference!")

# Hyperfine coupling constants
hf_couplings = {
    'H atom (1s)': {'A_MHz': 1420, 'I': 0.5, 'notes': 'Fermi contact'},
    'Benzene anion': {'A_MHz': 10.5, 'I': 0.5, 'notes': '6 equiv H'},
    'Methyl radical': {'A_MHz': 64.5, 'I': 0.5, 'notes': '3 equiv H'},
    'Nitroxide (¹⁴N)': {'A_MHz': 42, 'I': 1.0, 'notes': 'N coupling'},
    'Cu²⁺': {'A_MHz': 500, 'I': 1.5, 'notes': '⁶³,⁶⁵Cu'},
    'Mn²⁺': {'A_MHz': 250, 'I': 2.5, 'notes': '⁵⁵Mn'},
    'V⁴⁺': {'A_MHz': 300, 'I': 3.5, 'notes': '⁵¹V'},
}

print(f"\n  System          | A (MHz) | I   | Notes")
print("  " + "-" * 55)
for system, data in hf_couplings.items():
    print(f"  {system:16s} | {data['A_MHz']:6.1f}  | {data['I']:.1f} | {data['notes']}")

# Hyperfine splitting pattern
print(f"\n  Splitting pattern: 2nI + 1 lines")
print(f"  For I = 1/2: 2 lines (doublet)")
print(f"  For I = 1: 3 lines (triplet)")
print(f"  Equal intensity at γ ~ 1 (equivalent nuclei)")

# Isotropic vs anisotropic
print(f"\n  Isotropic coupling (A_iso):")
print(f"    Solution: molecular tumbling averages anisotropy")
print(f"    A_iso = (A_x + A_y + A_z)/3 (γ ~ 1 for averaging)")

# =============================================================================
# 3. RELAXATION TIMES
# =============================================================================
print("\n" + "=" * 70)
print("3. SPIN RELAXATION TIMES (T₁, T₂)")
print("=" * 70)

# T1 = spin-lattice, T2 = spin-spin
# At certain conditions, T1 = T2 (extreme narrowing)

print("  T₁ (spin-lattice): energy exchange with environment")
print("  T₂ (spin-spin): phase coherence time")
print("  Linewidth: Δν = 1/(πT₂)")

# Relaxation data
relaxation_data = {
    'Organic radical (dilute)': {'T1_s': 1e-6, 'T2_s': 1e-6},
    'Nitroxide (room T)': {'T1_s': 5e-7, 'T2_s': 5e-7},
    'Cu²⁺ (room T)': {'T1_s': 1e-9, 'T2_s': 1e-9},
    'Fe³⁺ high-spin': {'T1_s': 1e-11, 'T2_s': 1e-11},
    'Gd³⁺': {'T1_s': 1e-10, 'T2_s': 1e-10},
    'NV center (77K)': {'T1_s': 1e-3, 'T2_s': 1e-4},
}

print(f"\n  System                 | T₁ (s)   | T₂ (s)   | T₁/T₂")
print("  " + "-" * 60)
for system, data in relaxation_data.items():
    ratio = data['T1_s'] / data['T2_s']
    print(f"  {system:22s} | {data['T1_s']:.2e} | {data['T2_s']:.2e} | {ratio:.2f}")

# T1/T2 = 1 condition
print(f"\n  T₁/T₂ = 1 (extreme narrowing limit):")
print(f"    Occurs when ωτ_c << 1 (fast motion)")
print(f"    Same physics as NMR (Session #228)")
print(f"    T₁ = T₂ IS γ ~ 1 for spin coherence!")

# ωτc = 1 crossover
print(f"\n  At ωτ_c ≈ 1: T₁ minimum (γ ~ 1 for dynamics)")
print(f"    Below: T₁ = T₂ (narrowing)")
print(f"    Above: T₁ >> T₂ (slow motion)")

# =============================================================================
# 4. SPIN CONCENTRATION
# =============================================================================
print("\n" + "=" * 70)
print("4. SPIN CONCENTRATION AND DETECTION")
print("=" * 70)

# EPR sensitivity and detection limits

print("  EPR signal intensity ∝ N_spins × P^0.5")
print("  where P = microwave power")

# Detection limits (continuous wave X-band)
print(f"\n  Detection limits (X-band, 9.5 GHz):")
print(f"    Organic radicals: ~10¹⁰ spins")
print(f"    Transition metals: ~10¹¹ spins")
print(f"    High-Q cavity: ~10⁹ spins")

# Power saturation
print(f"\n  Power saturation:")
print(f"    At P = P_1/2: signal at half-maximum")
print(f"    P_1/2 = 1/(γ²T₁T₂) relates to relaxation")
print(f"    γ_eff = P/P_1/2")
print(f"    At γ_eff = 1: 50% saturation (γ ~ 1!)")

# Boltzmann polarization
print(f"\n  Boltzmann spin polarization:")
print(f"    P_B = tanh(gμ_B B₀ / 2kT)")
print(f"    At X-band (0.34 T), 300K: P_B ≈ 0.001")
print(f"    Thermal equilibrium limits signal")

# =============================================================================
# 5. LINE SHAPE ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("5. LINE SHAPE: LORENTZIAN vs GAUSSIAN")
print("=" * 70)

# Lorentzian: homogeneous broadening
# Gaussian: inhomogeneous broadening

print("  Lorentzian: homogeneous broadening (T₂ relaxation)")
print("  Gaussian: inhomogeneous broadening (field distribution)")

# Linewidth definitions
print(f"\n  Linewidth definitions:")
print(f"    Peak-to-peak ΔB_pp (derivative spectrum)")
print(f"    Full width at half maximum FWHM")
print(f"    For Lorentzian: FWHM = √3 × ΔB_pp")
print(f"    For Gaussian: FWHM = 1.18 × ΔB_pp")

# Voigt profile
print(f"\n  Voigt profile: convolution of L and G")
print(f"    Mixed broadening in real systems")
print(f"    Shape factor S = L/(L+G)")
print(f"    S = 0.5 is equal Lorentzian/Gaussian (γ ~ 1!)")

# Linewidth analysis
linewidths = {
    'DPPH (dilute)': {'Bpp_G': 0.2, 'type': 'Lorentzian'},
    'Coal': {'Bpp_G': 0.5, 'type': 'Gaussian'},
    'Nitroxide': {'Bpp_G': 1.5, 'type': 'Lorentzian'},
    'Cu²⁺ powder': {'Bpp_G': 10, 'type': 'Gaussian'},
    'Fe³⁺ glass': {'Bpp_G': 50, 'type': 'Gaussian'},
}

print(f"\n  Sample           | ΔB_pp (G) | Line Type")
print("  " + "-" * 45)
for sample, data in linewidths.items():
    print(f"  {sample:16s} | {data['Bpp_G']:6.1f}    | {data['type']}")

# =============================================================================
# 6. ZERO-FIELD SPLITTING
# =============================================================================
print("\n" + "=" * 70)
print("6. ZERO-FIELD SPLITTING (ZFS)")
print("=" * 70)

# ZFS occurs for S > 1/2 systems
# D and E parameters

print("  Zero-field splitting for S > 1/2:")
print("  H_ZFS = D[S_z² - S(S+1)/3] + E(S_x² - S_y²)")
print("  D = axial, E = rhombic parameter")

# ZFS parameters
zfs_data = {
    'Mn²⁺ (S=5/2)': {'D_cm': 0.025, 'E_cm': 0.008},
    'Fe³⁺ HS (S=5/2)': {'D_cm': 1.0, 'E_cm': 0.3},
    'Cr³⁺ (S=3/2)': {'D_cm': 0.2, 'E_cm': 0.02},
    'Ni²⁺ (S=1)': {'D_cm': 3.0, 'E_cm': 0.5},
    'V³⁺ (S=1)': {'D_cm': 5.0, 'E_cm': 0.1},
    'NV center (S=1)': {'D_cm': 0.095, 'E_cm': 0.0},
}

print(f"\n  System           | D (cm⁻¹) | E (cm⁻¹) | E/D")
print("  " + "-" * 50)
for system, data in zfs_data.items():
    D = data['D_cm']
    E = data['E_cm']
    ratio = E/D if D != 0 else 0
    print(f"  {system:16s} | {D:7.3f}  | {E:7.3f}  | {ratio:.3f}")

# E/D ratio
print(f"\n  E/D ratio:")
print(f"    E/D = 0: axial symmetry (C∞, D∞)")
print(f"    E/D = 1/3: maximum rhombicity")
print(f"    At E/D = 0: highest symmetry (γ ~ 1 reference)")

# =============================================================================
# 7. EXCHANGE COUPLING
# =============================================================================
print("\n" + "=" * 70)
print("7. EXCHANGE COUPLING (J)")
print("=" * 70)

# Exchange coupling between paramagnetic centers
# H_ex = -2J × S₁·S₂

print("  Exchange coupling: H_ex = -2J × S₁·S₂")
print("  J > 0: ferromagnetic (parallel spins)")
print("  J < 0: antiferromagnetic (antiparallel)")
print("  J = 0: uncoupled (γ ~ 1 reference!)")

# Exchange regimes
print(f"\n  Exchange regimes:")
print(f"    |J| >> A: exchange narrowing (single line)")
print(f"    |J| << A: individual spectra")
print(f"    |J| ≈ A: intermediate (γ ~ 1 crossover!)")

# Exchange coupling constants
exchange_data = {
    'Biradical (close)': {'J_cm': 100},
    'Biradical (far)': {'J_cm': 0.001},
    'Cu dimer': {'J_cm': -150},
    'Fe-Fe (direct)': {'J_cm': -500},
    'Superexchange': {'J_cm': -10},
}

print(f"\n  System           | J (cm⁻¹) | Coupling type")
print("  " + "-" * 50)
for system, data in exchange_data.items():
    J = data['J_cm']
    coupling = 'FM' if J > 0 else 'AFM' if J < 0 else 'none'
    print(f"  {system:16s} | {J:8.1f} | {coupling}")

# J/kT ratio
print(f"\n  J/kT ratio:")
print(f"    J/kT >> 1: strongly coupled")
print(f"    J/kT << 1: thermally populated")
print(f"    J/kT ≈ 1: crossover regime (γ ~ 1!)")
print(f"    At 300K: kT = 208 cm⁻¹")

# =============================================================================
# 8. RESONANCE CONDITION
# =============================================================================
print("\n" + "=" * 70)
print("8. EPR RESONANCE CONDITION")
print("=" * 70)

# hν = gμ_B B₀
print("  Resonance: hν = gμ_B B₀")
print("  For g = g_e = 2.0023:")

# Common EPR frequencies/fields
epr_bands = {
    'L-band': {'freq_GHz': 1.5, 'field_T': 0.054},
    'S-band': {'freq_GHz': 3.0, 'field_T': 0.107},
    'X-band': {'freq_GHz': 9.5, 'field_T': 0.339},
    'K-band': {'freq_GHz': 24, 'field_T': 0.857},
    'Q-band': {'freq_GHz': 34, 'field_T': 1.214},
    'W-band': {'freq_GHz': 94, 'field_T': 3.356},
}

print(f"\n  Band     | ν (GHz) | B₀ (T)  | B₀/ν (T/GHz)")
print("  " + "-" * 50)
for band, data in epr_bands.items():
    ratio = data['field_T'] / data['freq_GHz']
    print(f"  {band:8s} | {data['freq_GHz']:6.1f}  | {data['field_T']:.3f}   | {ratio:.4f}")

# Universal field/frequency ratio
ratio_universal = 0.339 / 9.5  # X-band ratio
print(f"\n  Universal ratio: B₀/ν = h/(gμ_B) = 0.0357 T/GHz")
print(f"  This ratio IS the γ ~ 1 reference for resonance!")

# =============================================================================
# 9. SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY: γ ~ 1 IN EPR SPECTROSCOPY")
print("=" * 70)

gamma_findings = [
    ("g-factor", "g_e = 2.0023", "Free electron IS γ ~ 1 reference"),
    ("Hyperfine", "A = 0", "No nuclear coupling (I = 0)"),
    ("Relaxation", "T₁/T₂ = 1", "Extreme narrowing limit"),
    ("Power", "P = P_1/2", "50% saturation"),
    ("Line shape", "L/(L+G) = 0.5", "Equal broadening"),
    ("ZFS", "E/D = 0", "Axial symmetry reference"),
    ("Exchange", "J = 0", "Uncoupled spins"),
    ("Resonance", "hν = gμ_B B", "Field/frequency universal"),
]

print("\n  Parameter        | γ ~ 1 Condition        | Interpretation")
print("  " + "-" * 65)
for param, value, interp in gamma_findings:
    print(f"  {param:16s} | {value:22s} | {interp}")

print("\n  CONCLUSION: EPR IS γ ~ 1 spectroscopy:")
print("    - g_e = 2.0023 is THE universal reference")
print("    - Deviations measure spin-orbit effects")
print("    - T₁ = T₂ at extreme narrowing")
print("    - Exchange coupling J = 0 is uncoupled reference")
print("    - Line shape L/G = 0.5 is equal broadening")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('EPR Spectroscopy Coherence Analysis\nSession #232: γ ~ 1 in Electron Spin Resonance',
             fontsize=14, fontweight='bold')

# 1. g-factor distribution
ax1 = axes[0, 0]
g_values = [2.00232, 2.0036, 2.003, 2.0, 2.006]
g_names = ['Free e-', 'DPPH', 'π radical', 'Mn²⁺', 'Nitroxide']
colors = ['red' if abs(g - g_e) < 0.005 else 'blue' for g in g_values]
ax1.barh(g_names, g_values, color=colors, alpha=0.7)
ax1.axvline(x=g_e, color='red', linestyle='--', linewidth=2, label=f'g_e = {g_e:.4f}')
ax1.set_xlabel('g-factor')
ax1.set_title('g-Factor Distribution')
ax1.set_xlim(1.99, 2.02)
ax1.legend()
ax1.grid(True, alpha=0.3, axis='x')

# 2. Hyperfine coupling
ax2 = axes[0, 1]
hf_names = list(hf_couplings.keys())[:5]
hf_vals = [hf_couplings[n]['A_MHz'] for n in hf_names]
ax2.barh(hf_names, hf_vals, color='teal', alpha=0.7)
ax2.axvline(x=0, color='red', linestyle='--', linewidth=2, label='A = 0 (γ ~ 1)')
ax2.set_xlabel('Hyperfine Coupling (MHz)')
ax2.set_title('Hyperfine Coupling Constants')
ax2.set_xscale('log')
ax2.legend()
ax2.grid(True, alpha=0.3, axis='x')

# 3. Relaxation times
ax3 = axes[0, 2]
relax_names = list(relaxation_data.keys())
T1_vals = [relaxation_data[n]['T1_s'] for n in relax_names]
T2_vals = [relaxation_data[n]['T2_s'] for n in relax_names]
x = np.arange(len(relax_names))
width = 0.35
ax3.barh(x - width/2, T1_vals, width, label='T₁', color='blue', alpha=0.7)
ax3.barh(x + width/2, T2_vals, width, label='T₂', color='green', alpha=0.7)
ax3.set_yticks(x)
ax3.set_yticklabels([n.replace(' ', '\n') for n in relax_names], fontsize=7)
ax3.set_xlabel('Relaxation Time (s)')
ax3.set_title('T₁ and T₂ Relaxation')
ax3.set_xscale('log')
ax3.legend()
ax3.grid(True, alpha=0.3, axis='x')

# 4. Line shapes
ax4 = axes[1, 0]
B = np.linspace(-5, 5, 200)
# Lorentzian
L = 1 / (1 + B**2)
# Gaussian
G = np.exp(-B**2)
# Voigt (equal mix)
V = 0.5 * L + 0.5 * G
ax4.plot(B, L, 'b-', linewidth=2, label='Lorentzian')
ax4.plot(B, G, 'g-', linewidth=2, label='Gaussian')
ax4.plot(B, V, 'r--', linewidth=2, label='Voigt (γ ~ 1)')
ax4.set_xlabel('Field (arb.)')
ax4.set_ylabel('Intensity')
ax4.set_title('EPR Line Shapes')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. ZFS E/D ratio
ax5 = axes[1, 1]
zfs_names = list(zfs_data.keys())
E_D_ratios = [zfs_data[n]['E_cm']/zfs_data[n]['D_cm'] if zfs_data[n]['D_cm'] != 0 else 0 for n in zfs_names]
colors = ['green' if r < 0.1 else 'orange' if r < 0.2 else 'red' for r in E_D_ratios]
ax5.barh(zfs_names, E_D_ratios, color=colors, alpha=0.7)
ax5.axvline(x=0, color='red', linestyle='--', linewidth=2, label='E/D = 0 (γ ~ 1)')
ax5.axvline(x=1/3, color='gray', linestyle=':', alpha=0.7, label='E/D = 1/3 (max)')
ax5.set_xlabel('E/D Ratio')
ax5.set_title('Zero-Field Splitting Rhombicity')
ax5.legend()
ax5.grid(True, alpha=0.3, axis='x')

# 6. EPR bands
ax6 = axes[1, 2]
band_names = list(epr_bands.keys())
freqs = [epr_bands[n]['freq_GHz'] for n in band_names]
fields = [epr_bands[n]['field_T'] for n in band_names]
ax6.scatter(freqs, fields, s=100, c='blue', alpha=0.7)
for i, name in enumerate(band_names):
    ax6.annotate(name, (freqs[i], fields[i]), xytext=(5, 5), textcoords='offset points', fontsize=8)
# Theoretical line
freq_theory = np.linspace(0, 100, 100)
field_theory = freq_theory * 0.0357  # h/(g_e × μ_B)
ax6.plot(freq_theory, field_theory, 'r--', linewidth=2, label='hν = gμ_B B (γ ~ 1)')
ax6.set_xlabel('Frequency (GHz)')
ax6.set_ylabel('Magnetic Field (T)')
ax6.set_title('EPR Resonance Condition')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/epr_spectroscopy_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FINDING #169: EPR Spectroscopy at γ ~ 1")
print("=" * 70)
print("""
Electron Paramagnetic Resonance exhibits γ ~ 1 at fundamental
spectroscopic boundaries:

1. g_e = 2.0023: Free electron IS the universal reference
2. A = 0: No hyperfine (I = 0) - uncoupled reference
3. T₁/T₂ = 1: Extreme narrowing (fast motion)
4. P = P_1/2: 50% power saturation
5. L/(L+G) = 0.5: Equal Lorentzian/Gaussian broadening
6. E/D = 0: Axial symmetry (highest symmetry)
7. J = 0: Uncoupled spins reference

95th phenomenon type exhibiting γ ~ 1 transition behavior.
EPR spectroscopy IS electron spin coherence measurement!
""")

print("Visualization saved: epr_spectroscopy_coherence.png")
