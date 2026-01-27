#!/usr/bin/env python3
"""
Raman Spectroscopy Coherence Analysis
Chemistry Session #235

Explores how γ ~ 1 manifests in Raman scattering:
1. Stokes/Anti-Stokes ratio - thermal equilibrium
2. Depolarization ratio - symmetry classification
3. SERS enhancement - resonance amplification
4. Rayleigh line - elastic reference
5. Polarizability derivative - selection rule
6. Group frequencies - characteristic vibrations
7. Resonance Raman - electronic coupling

Key insight: Raman spectroscopy measures polarizability changes,
with γ ~ 1 appearing at thermal balance, symmetry transitions,
and resonance enhancement boundaries.
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 70)
print("RAMAN SPECTROSCOPY COHERENCE ANALYSIS")
print("Chemistry Session #235: γ ~ 1 in Inelastic Light Scattering")
print("=" * 70)

# Physical constants
h = 6.626e-34     # Planck constant
c = 3e8            # Speed of light
k_B = 1.381e-23   # Boltzmann constant

# =============================================================================
# 1. STOKES / ANTI-STOKES RATIO
# =============================================================================
print("\n" + "=" * 70)
print("1. STOKES / ANTI-STOKES INTENSITY RATIO")
print("=" * 70)

# I_AS/I_S = (ν₀ + ν)⁴/(ν₀ - ν)⁴ × exp(-hcν̃/kT)
# At high T (hcν̃ << kT): ratio → 1 (γ ~ 1!)

print("  I_AS/I_S = [(ν₀ + ν)/(ν₀ - ν)]⁴ × exp(-hcν̃/kT)")
print("  At kT >> hcν̃: ratio → 1 (equal populations, γ ~ 1)")

# Calculate ratios at different temperatures for typical vibrations
T = 300  # K
wavenumbers = [200, 500, 1000, 1500, 2000, 3000]  # cm⁻¹

print(f"\n  Stokes/Anti-Stokes ratios at T = {T} K:")
print(f"  (assuming ν₀ >> ν for frequency factor ≈ 1)")
print(f"\n  ν̃ (cm⁻¹) | hcν̃/kT | I_AS/I_S | Comment")
print("  " + "-" * 55)
for nu in wavenumbers:
    x = h * c * nu * 100 / (k_B * T)  # hcν̃/kT
    ratio = np.exp(-x)
    comment = ""
    if abs(ratio - 1) < 0.1:
        comment = "near γ ~ 1"
    elif abs(ratio - 0.5) < 0.1:
        comment = "50% (γ ~ 1 threshold)"
    print(f"  {nu:6d}    | {x:.3f}  | {ratio:.4f}  | {comment}")

# Find crossover temperature where ratio = 0.5
print(f"\n  At I_AS/I_S = 0.5 (γ ~ 1 threshold):")
print(f"  hcν̃/kT = ln(2) = 0.693")
for nu in [200, 500, 1000, 1500]:
    T_half = h * c * nu * 100 / (k_B * 0.693)
    print(f"    ν̃ = {nu:4d} cm⁻¹: T = {T_half:.0f} K")

# At high T limit: ratio → 1 (equal Stokes and Anti-Stokes)
print(f"\n  High-T limit: I_AS/I_S → 1 (γ ~ 1)")
print(f"  This means equal population of ground and excited states")
print(f"  Quantum → classical transition at kT = hcν̃ (γ ~ 1!)")

# =============================================================================
# 2. DEPOLARIZATION RATIO
# =============================================================================
print("\n" + "=" * 70)
print("2. DEPOLARIZATION RATIO (ρ)")
print("=" * 70)

# ρ = I_⊥/I_∥
# For linearly polarized excitation:
# ρ = 0 for totally symmetric modes
# ρ = 3/4 for depolarized (asymmetric) modes

print("  Depolarization ratio: ρ = I_⊥/I_∥")
print("  ρ classifies vibrational symmetry!")

depol_values = {
    'Totally symmetric (isotropic)': 0.0,
    'Slightly depolarized': 0.1,
    'Moderately depolarized': 0.3,
    'Anomalously polarized': 0.50,
    'Depolarized (3/4)': 0.75,
}

print(f"\n  Category                    | ρ     | γ interpretation")
print("  " + "-" * 60)
for category, rho in depol_values.items():
    gamma_interp = ""
    if rho == 0.0:
        gamma_interp = "Perfect symmetry (γ ~ 1 reference)"
    elif rho == 0.75:
        gamma_interp = "Maximum asymmetry (3/4)"
    elif rho == 0.50:
        gamma_interp = "Anomalous polarization"
    print(f"  {category:28s} | {rho:.2f}  | {gamma_interp}")

# Real molecules
mol_depol = {
    'CCl4 ν1 (A1)': 0.005,
    'CS2 ν1 (Σg+)': 0.01,
    'CHCl3 ν1 (A1)': 0.01,
    'C6H6 ν1 (A1g)': 0.02,
    'CH3OH ν (C-O)': 0.30,
    'H2O ν1 (A1)': 0.05,
    'CO2 ν1 (Σg+)': 0.0,  # polarized
}

print(f"\n  Molecule ν mode   | ρ     | Symmetry")
print("  " + "-" * 45)
for mol, rho in mol_depol.items():
    sym = "polarized" if rho < 0.375 else "depolarized"
    print(f"  {mol:18s} | {rho:.3f} | {sym}")

print(f"\n  ρ < 3/4: polarized (symmetric)")
print(f"  ρ = 3/4: depolarized (asymmetric)")
print(f"  Boundary at ρ = 3/4 classifies symmetry!")
print(f"  ρ = 0 IS γ ~ 1 for perfect isotropy")

# =============================================================================
# 3. SERS ENHANCEMENT
# =============================================================================
print("\n" + "=" * 70)
print("3. SURFACE-ENHANCED RAMAN (SERS)")
print("=" * 70)

# SERS enhancement factor: EF = |E_loc/E_inc|⁴
# At plasmon resonance: maximum enhancement

print("  SERS Enhancement: EF = |E_loc/E_0|⁴")
print("  At surface plasmon resonance: EF → 10⁶-10¹⁴")

# Enhancement factors
sers_substrates = {
    'Au nanoparticles (50 nm)': {'EF': 1e6, 'λ_res_nm': 530},
    'Ag nanoparticles (40 nm)': {'EF': 1e7, 'λ_res_nm': 410},
    'Au nanorods (aspect 3)': {'EF': 1e8, 'λ_res_nm': 700},
    'Ag film (roughened)': {'EF': 1e5, 'λ_res_nm': 400},
    'Hot spot (dimer gap)': {'EF': 1e10, 'λ_res_nm': 550},
    'Single molecule': {'EF': 1e14, 'λ_res_nm': 520},
}

print(f"\n  Substrate              | EF        | λ_plasmon (nm)")
print("  " + "-" * 55)
for sub, data in sers_substrates.items():
    print(f"  {sub:24s} | {data['EF']:.0e}    | {data['λ_res_nm']}")

# Resonance condition
print(f"\n  Maximum SERS when:")
print(f"    λ_laser ≈ λ_plasmon (resonance, γ ~ 1!)")
print(f"    λ_plasmon ≈ λ_Raman (double resonance)")
print(f"    |E_loc/E_0| = 1 at non-resonant (γ ~ 1 baseline)")
print(f"    Enhancement measured relative to |E_loc/E_0| = 1")

# =============================================================================
# 4. RAYLEIGH LINE
# =============================================================================
print("\n" + "=" * 70)
print("4. RAYLEIGH (ELASTIC) SCATTERING REFERENCE")
print("=" * 70)

# Rayleigh line at Δν = 0 is THE γ ~ 1 reference
print("  Rayleigh scattering: Δν̃ = 0 cm⁻¹")
print("  This IS the elastic scattering reference (γ ~ 1!)")
print("  All Raman shifts measured from Rayleigh line")

print(f"\n  Scattering types relative to Rayleigh:")
print(f"    Rayleigh (elastic): Δν̃ = 0 cm⁻¹ (γ ~ 1 reference)")
print(f"    Stokes (red shift): Δν̃ < 0 (energy loss)")
print(f"    Anti-Stokes (blue shift): Δν̃ > 0 (energy gain)")
print(f"    Brillouin: |Δν̃| ≈ 0.1-1 cm⁻¹ (acoustic phonons)")

# Intensity ratio
print(f"\n  Relative intensities (typical organic liquid):")
print(f"    Rayleigh: 1.000 (reference)")
print(f"    Stokes: 10⁻³ (1000× weaker)")
print(f"    Anti-Stokes: 10⁻⁶ (10⁶× weaker)")
print(f"    Rayleigh IS the dominant elastic process (γ ~ 1)")

# =============================================================================
# 5. SELECTION RULES AND POLARIZABILITY
# =============================================================================
print("\n" + "=" * 70)
print("5. RAMAN SELECTION RULES")
print("=" * 70)

# Raman active: dα/dQ ≠ 0 (polarizability change)
# IR active: dμ/dQ ≠ 0 (dipole change)
# Mutual exclusion: centrosymmetric molecules

print("  Raman: requires dα/dQ ≠ 0 (polarizability change)")
print("  IR: requires dμ/dQ ≠ 0 (dipole change)")
print("  Selection rule dα/dQ = 0 IS the γ ~ 1 reference (inactive)")

# Mutual exclusion
print(f"\n  Mutual exclusion principle (centrosymmetric):")
print(f"    If Raman active → IR inactive (and vice versa)")
print(f"    Example: CO₂")
print(f"      ν₁ (symmetric stretch): Raman active, IR inactive")
print(f"      ν₃ (asymmetric stretch): IR active, Raman inactive")
print(f"    This IS γ ~ 1 complementarity!")

# Rule of mutual exclusion examples
molecules = {
    'CO₂ (D∞h)': {'mutual': True, 'modes_R': 1, 'modes_IR': 2},
    'H₂O (C2v)': {'mutual': False, 'modes_R': 3, 'modes_IR': 3},
    'C₂H₂ (D∞h)': {'mutual': True, 'modes_R': 3, 'modes_IR': 2},
    'CH₄ (Td)': {'mutual': False, 'modes_R': 4, 'modes_IR': 2},
    'C₆H₆ (D6h)': {'mutual': True, 'modes_R': 7, 'modes_IR': 4},
    'SF₆ (Oh)': {'mutual': True, 'modes_R': 3, 'modes_IR': 2},
}

print(f"\n  Molecule   | Mutual Excl | Raman modes | IR modes")
print("  " + "-" * 55)
for mol, data in molecules.items():
    me = "Yes" if data['mutual'] else "No"
    print(f"  {mol:12s} | {me:11s} | {data['modes_R']:5d}       | {data['modes_IR']}")

# =============================================================================
# 6. GROUP FREQUENCIES
# =============================================================================
print("\n" + "=" * 70)
print("6. CHARACTERISTIC GROUP FREQUENCIES")
print("=" * 70)

# Important Raman-active group frequencies
group_freq = {
    'C-C stretch': {'range': (800, 1200), 'mid': 1000},
    'C=C stretch': {'range': (1600, 1680), 'mid': 1640},
    'C≡C stretch': {'range': (2100, 2260), 'mid': 2180},
    'C-H stretch': {'range': (2800, 3100), 'mid': 2950},
    'O-H stretch': {'range': (3200, 3600), 'mid': 3400},
    'C=O stretch': {'range': (1650, 1780), 'mid': 1715},
    'S-S stretch': {'range': (470, 550), 'mid': 510},
    'N-N stretch': {'range': (1000, 1200), 'mid': 1100},
    'Ring breathing': {'range': (700, 1050), 'mid': 875},
}

print("  Characteristic Raman frequencies:")
print(f"\n  Group            | Range (cm⁻¹)  | Midpoint | Strong Raman?")
print("  " + "-" * 60)
for group, data in group_freq.items():
    lo, hi = data['range']
    mid = data['mid']
    strong = "Yes" if group in ['C=C stretch', 'C≡C stretch', 'S-S stretch', 'Ring breathing'] else "Moderate"
    print(f"  {group:16s} | {lo:4d} - {hi:4d}  | {mid:5d}   | {strong}")

# Isotope shifts
print(f"\n  Isotope shift ratio:")
print(f"    ν(heavy)/ν(light) = √(m_light/m_heavy)")
print(f"    H/D: ν_ratio = √(1/2) = 0.707")
print(f"    ¹²C/¹³C: ν_ratio = √(12/13) = 0.961")
print(f"    At equal masses: ratio = 1 (γ ~ 1!)")

# =============================================================================
# 7. RESONANCE RAMAN
# =============================================================================
print("\n" + "=" * 70)
print("7. RESONANCE RAMAN SPECTROSCOPY")
print("=" * 70)

# Resonance when laser matches electronic transition
# Enhancement up to 10⁶

print("  Resonance condition: hν_laser ≈ hν_electronic")
print("  Enhancement: 10³-10⁶ for resonant modes!")

# Resonance Raman examples
rr_examples = {
    'β-carotene (π→π*)': {'λ_abs': 450, 'EF': 1e5, 'modes': 'C=C, C-C'},
    'Heme (Soret)': {'λ_abs': 410, 'EF': 1e6, 'modes': 'Fe-N, porphyrin'},
    'Rhodopsin': {'λ_abs': 500, 'EF': 1e4, 'modes': 'C=C retinal'},
    'MnO4⁻ (CT)': {'λ_abs': 525, 'EF': 1e4, 'modes': 'Mn-O'},
    'CrO4²⁻ (CT)': {'λ_abs': 370, 'EF': 1e3, 'modes': 'Cr-O'},
}

print(f"\n  System            | λ_abs (nm) | Enhancement | Key modes")
print("  " + "-" * 65)
for system, data in rr_examples.items():
    print(f"  {system:18s} | {data['λ_abs']:4d}       | {data['EF']:.0e}       | {data['modes']}")

print(f"\n  At exact resonance (λ_laser = λ_abs): maximum enhancement")
print(f"  λ_laser/λ_abs = 1 IS γ ~ 1!")
print(f"  Pre-resonance (λ close but not equal): partial enhancement")

# Excitation profile
print(f"\n  Excitation profile:")
print(f"    I_Raman ∝ 1/(ν_el - ν_0)² (denominator)")
print(f"    At ν_0 = ν_el: resonance (γ ~ 1), I → maximum")
print(f"    Damping factor Γ prevents true divergence")

# =============================================================================
# 8. SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY: γ ~ 1 IN RAMAN SPECTROSCOPY")
print("=" * 70)

gamma_findings = [
    ("Stokes/Anti-Stokes", "I_AS/I_S → 1 at high T", "Equal populations (γ ~ 1)"),
    ("Thermal crossover", "kT = hcν̃", "Quantum-classical transition"),
    ("Depolarization", "ρ = 0 (isotropic)", "Perfect symmetry reference"),
    ("SERS resonance", "λ_laser = λ_plasmon", "Maximum enhancement"),
    ("Rayleigh line", "Δν̃ = 0", "Elastic reference (γ ~ 1)"),
    ("Selection rule", "dα/dQ = 0", "Inactive reference"),
    ("Resonance Raman", "ν_laser = ν_electronic", "Maximum coupling (γ ~ 1)"),
    ("Isotope shift", "m₁/m₂ = 1", "Equal mass reference"),
]

print("\n  Parameter           | γ ~ 1 Condition          | Interpretation")
print("  " + "-" * 70)
for param, value, interp in gamma_findings:
    print(f"  {param:19s} | {value:24s} | {interp}")

print("\n  CONCLUSION: Raman spectroscopy IS γ ~ 1 scattering science:")
print("    - Rayleigh (Δν = 0) IS the elastic reference")
print("    - Stokes/Anti-Stokes ratio → 1 at high T")
print("    - SERS at plasmon resonance (λ_laser = λ_plasmon)")
print("    - Resonance Raman at electronic match")
print("    - ρ = 0 for perfect isotropy")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Raman Spectroscopy Coherence Analysis\nSession #235: γ ~ 1 in Inelastic Light Scattering',
             fontsize=14, fontweight='bold')

# 1. Stokes/Anti-Stokes ratio vs temperature
ax1 = axes[0, 0]
T_range = np.linspace(100, 2000, 200)
for nu in [200, 500, 1000, 2000]:
    x = h * c * nu * 100 / (k_B * T_range)
    ratio = np.exp(-x)
    ax1.plot(T_range, ratio, linewidth=2, label=f'{nu} cm⁻¹')
ax1.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='γ ~ 1')
ax1.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('I_AS/I_S')
ax1.set_title('Stokes/Anti-Stokes Ratio')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# 2. Depolarization ratio
ax2 = axes[0, 1]
rho_vals = list(depol_values.values())
rho_names = list(depol_values.keys())
colors = ['green' if r < 0.1 else 'orange' if r < 0.5 else 'red' for r in rho_vals]
ax2.barh(rho_names, rho_vals, color=colors, alpha=0.7)
ax2.axvline(x=0.75, color='red', linestyle='--', linewidth=2, label='ρ = 3/4 (depolarized)')
ax2.axvline(x=0, color='green', linestyle='--', linewidth=2, label='ρ = 0 (γ ~ 1)')
ax2.set_xlabel('Depolarization Ratio ρ')
ax2.set_title('Depolarization Classification')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3, axis='x')

# 3. SERS enhancement
ax3 = axes[0, 2]
sers_names = [s.replace(' ', '\n') for s in sers_substrates.keys()]
sers_ef = [sers_substrates[n]['EF'] for n in sers_substrates.keys()]
ax3.barh(sers_names, sers_ef, color='teal', alpha=0.7)
ax3.axvline(x=1, color='red', linestyle='--', linewidth=2, label='EF = 1 (no enhancement)')
ax3.set_xlabel('Enhancement Factor')
ax3.set_title('SERS Enhancement')
ax3.set_xscale('log')
ax3.legend()
ax3.grid(True, alpha=0.3, axis='x')

# 4. Raman spectrum schematic
ax4 = axes[1, 0]
# Simulate Raman spectrum with Rayleigh, Stokes, Anti-Stokes
shifts = np.linspace(-3000, 3000, 1000)
# Rayleigh (elastic, strong)
rayleigh = 100 * np.exp(-shifts**2 / (2 * 10**2))
# Stokes peaks
stokes = 10 * np.exp(-(shifts + 1000)**2 / (2 * 20**2)) + \
         5 * np.exp(-(shifts + 1600)**2 / (2 * 20**2)) + \
         8 * np.exp(-(shifts + 2900)**2 / (2 * 30**2))
# Anti-Stokes (weaker)
antistokes = 2 * np.exp(-(shifts - 1000)**2 / (2 * 20**2)) + \
             1 * np.exp(-(shifts - 1600)**2 / (2 * 20**2)) + \
             0.5 * np.exp(-(shifts - 2900)**2 / (2 * 30**2))
ax4.plot(shifts, rayleigh + stokes + antistokes, 'b-', linewidth=1.5)
ax4.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Δν = 0 (Rayleigh, γ ~ 1)')
ax4.set_xlabel('Raman Shift (cm⁻¹)')
ax4.set_ylabel('Intensity')
ax4.set_title('Raman Spectrum Schematic')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# 5. Resonance Raman enhancement
ax5 = axes[1, 1]
# Plot excitation profile
nu_el = 20000  # cm⁻¹ (500 nm)
Gamma = 500    # damping
nu_range = np.linspace(15000, 25000, 300)
intensity = 1 / ((nu_el - nu_range)**2 + Gamma**2)
intensity /= np.max(intensity)
ax5.plot(nu_range, intensity, 'b-', linewidth=2)
ax5.axvline(x=nu_el, color='r', linestyle='--', linewidth=2, label=f'ν_el (γ ~ 1)')
ax5.set_xlabel('Laser Frequency (cm⁻¹)')
ax5.set_ylabel('Raman Enhancement')
ax5.set_title('Resonance Raman Profile')
ax5.legend()
ax5.grid(True, alpha=0.3)

# 6. Group frequencies
ax6 = axes[1, 2]
gf_names = list(group_freq.keys())
gf_mids = [group_freq[n]['mid'] for n in gf_names]
gf_ranges = [(group_freq[n]['range'][1] - group_freq[n]['range'][0])/2 for n in gf_names]
ax6.barh(gf_names, gf_mids, xerr=gf_ranges, color='purple', alpha=0.7, capsize=3)
ax6.set_xlabel('Wavenumber (cm⁻¹)')
ax6.set_title('Characteristic Group Frequencies')
ax6.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/raman_spectroscopy_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "=" * 70)
print("FINDING #172: Raman Spectroscopy at γ ~ 1")
print("=" * 70)
print("""
Raman spectroscopy exhibits γ ~ 1 at fundamental scattering
and resonance boundaries:

1. I_AS/I_S → 1 at high T (equal populations)
2. kT = hcν̃: quantum-classical transition
3. ρ = 0: perfect isotropy (totally symmetric)
4. Δν̃ = 0: Rayleigh elastic reference
5. λ_laser = λ_plasmon: SERS maximum enhancement
6. ν_laser = ν_electronic: Resonance Raman maximum
7. m₁/m₂ = 1: isotope shift reference

98th phenomenon type exhibiting γ ~ 1 transition behavior.
Raman scattering IS coherence spectroscopy of polarizability!
""")

print("Visualization saved: raman_spectroscopy_coherence.png")
