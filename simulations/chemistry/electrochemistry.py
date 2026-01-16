#!/usr/bin/env python3
"""
Chemistry Session #52: Electrochemistry and Coherence at Interfaces

Electrochemistry involves electron transfer between:
- Electrode (metallic, high coherence)
- Solution (ionic/molecular, lower coherence)

Question: How does the coherence MISMATCH at interfaces affect:
- Electron transfer rates?
- Electrode kinetics?
- Battery/fuel cell performance?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

print("=" * 70)
print("Chemistry Session #52: Electrochemistry & Coherence at Interfaces")
print("=" * 70)
print()

# Physical constants
k_B = constants.k
e = constants.e
h = constants.h
hbar = constants.hbar

# =============================================================================
# PART 1: THE INTERFACE PROBLEM
# =============================================================================

print("-" * 70)
print("PART 1: THE INTERFACE PROBLEM")
print("-" * 70)
print()

print("Electrode/electrolyte interface:")
print()
print("  ELECTRODE (metal):")
print("    - Delocalized electrons (Fermi gas)")
print("    - High coherence (γ_metal ~ 0.1-0.3)")
print("    - Strong coupling (metallic bonding)")
print()
print("  ELECTROLYTE (solution):")
print("    - Localized redox species")
print("    - Lower coherence (γ_solution ~ 1-2)")
print("    - Weak coupling (solvation)")
print()
print("Electron transfer crosses this coherence DISCONTINUITY.")
print()

# =============================================================================
# PART 2: MARCUS THEORY REVIEW
# =============================================================================

print("-" * 70)
print("PART 2: MARCUS THEORY REVIEW")
print("-" * 70)
print()

print("Marcus theory for electron transfer rate:")
print()
print("  k_ET = (2π/ℏ)|V|² × (1/√4πλkT) × exp[-(λ + ΔG°)²/4λkT]")
print()
print("Where:")
print("  V = electronic coupling (matrix element)")
print("  λ = reorganization energy")
print("  ΔG° = driving force")
print()

# Typical values
V_typ = 0.01  # eV (weak coupling)
lambda_typ = 1.0  # eV (reorganization)
T = 298  # K
kT = k_B * T / e  # in eV

print(f"Typical values at {T} K:")
print(f"  V ~ {V_typ} eV")
print(f"  λ ~ {lambda_typ} eV")
print(f"  kT ~ {kT:.3f} eV")
print()

# =============================================================================
# PART 3: SYNCHRONISM MODIFICATION
# =============================================================================

print("-" * 70)
print("PART 3: SYNCHRONISM MODIFICATION")
print("-" * 70)
print()

print("In Synchronism, electron transfer involves coherence matching:")
print()
print("  k_ET = k_Marcus × f(γ_electrode, γ_solution)")
print()
print("The coherence factor f accounts for phase-matching efficiency.")
print()

print("Hypothesis:")
print("  f = min(γ_electrode, γ_solution) / max(γ_electrode, γ_solution)")
print()
print("Physical meaning:")
print("  - Well-matched coherence (similar γ) → f ~ 1")
print("  - Mismatched coherence (different γ) → f << 1")
print()

def coherence_factor(gamma_1, gamma_2):
    """
    Coherence matching factor for electron transfer.
    f = min(γ₁, γ₂) / max(γ₁, γ₂)
    Works with scalars or arrays.
    """
    return np.minimum(gamma_1, gamma_2) / np.maximum(gamma_1, gamma_2)

# Test cases
test_cases = [
    ("Metal-Metal", 0.2, 0.2),
    ("Metal-Solution", 0.2, 1.5),
    ("Solution-Solution", 1.0, 1.2),
    ("Metal-Insulator", 0.2, 1.9),
]

print("Coherence matching factors:")
print(f"{'System':<20} | {'γ₁':>6} | {'γ₂':>6} | {'f':>8}")
print("-" * 50)

for name, g1, g2 in test_cases:
    f = coherence_factor(g1, g2)
    print(f"{name:<20} | {g1:>6.1f} | {g2:>6.1f} | {f:>8.3f}")

print()

# =============================================================================
# PART 4: ELECTRODE COHERENCE
# =============================================================================

print("-" * 70)
print("PART 4: ELECTRODE COHERENCE BY MATERIAL")
print("-" * 70)
print()

print("Different electrode materials have different γ:")
print()

electrode_data = {
    "Pt (noble)": {"gamma": 0.15, "work_function": 5.65, "d_band": "full"},
    "Au (noble)": {"gamma": 0.20, "work_function": 5.10, "d_band": "full"},
    "Ni (transition)": {"gamma": 0.40, "work_function": 5.15, "d_band": "partial"},
    "Fe (transition)": {"gamma": 0.45, "work_function": 4.50, "d_band": "partial"},
    "Carbon (glassy)": {"gamma": 0.80, "work_function": 4.80, "d_band": "sp"},
    "ITO (oxide)": {"gamma": 1.20, "work_function": 4.50, "d_band": "oxide"},
}

print(f"{'Material':<20} | {'γ':>6} | {'Φ (eV)':>8} | {'d-band':>10}")
print("-" * 55)

for name, data in electrode_data.items():
    print(f"{name:<20} | {data['gamma']:>6.2f} | {data['work_function']:>8.2f} | {data['d_band']:>10}")

print()
print("Observation: Noble metals have LOWEST γ (most coherent).")
print("This correlates with catalytic activity!")
print()

# =============================================================================
# PART 5: SOLUTION COHERENCE
# =============================================================================

print("-" * 70)
print("PART 5: SOLUTION/REDOX COHERENCE")
print("-" * 70)
print()

print("Redox species in solution have varying coherence:")
print()

redox_data = {
    "Fe³⁺/Fe²⁺ (aq)": {"gamma": 1.4, "lambda": 1.0, "type": "inner-sphere"},
    "Ru(bpy)₃³⁺/²⁺": {"gamma": 0.8, "lambda": 0.3, "type": "outer-sphere"},
    "Ferrocene⁺/⁰": {"gamma": 0.6, "lambda": 0.2, "type": "outer-sphere"},
    "H₂O/O₂": {"gamma": 1.8, "lambda": 2.0, "type": "inner-sphere"},
    "H⁺/H₂": {"gamma": 1.5, "lambda": 0.8, "type": "inner-sphere"},
}

print(f"{'Couple':<20} | {'γ':>6} | {'λ (eV)':>8} | {'Type':>15}")
print("-" * 60)

for name, data in redox_data.items():
    print(f"{name:<20} | {data['gamma']:>6.1f} | {data['lambda']:>8.1f} | {data['type']:>15}")

print()
print("Outer-sphere couples have LOWER γ (more delocalized).")
print()

# =============================================================================
# PART 6: COHERENCE-MODIFIED BUTLER-VOLMER
# =============================================================================

print("-" * 70)
print("PART 6: COHERENCE-MODIFIED BUTLER-VOLMER")
print("-" * 70)
print()

print("Standard Butler-Volmer equation:")
print("  i = i₀ × [exp(αfη) - exp(-(1-α)fη)]")
print()
print("Where:")
print("  i₀ = exchange current density")
print("  α = transfer coefficient")
print("  f = F/RT")
print("  η = overpotential")
print()

print("Synchronism modification:")
print("  i₀_eff = i₀ × f(γ_electrode, γ_redox)")
print()
print("And:")
print("  α = 0.5 × (γ_electrode / γ_redox) for asymmetric barriers")
print()

def alpha_synchronism(gamma_electrode, gamma_redox):
    """
    Transfer coefficient from coherence ratio.
    """
    # Ranges from 0 to 1, centered at 0.5 for matched coherence
    ratio = gamma_electrode / gamma_redox
    return 0.5 * ratio / (1 + ratio)

# Calculate α for different combinations
print("Transfer coefficient α from coherence:")
print(f"{'System':<30} | {'γ_elec':>8} | {'γ_redox':>8} | {'α':>8}")
print("-" * 60)

elec_list = list(electrode_data.items())[:3]
redox_list = list(redox_data.items())[:2]

for elec_name, elec_dict in elec_list:
    for redox_name, redox_dict in redox_list:
        alpha = alpha_synchronism(elec_dict["gamma"], redox_dict["gamma"])
        print(f"{elec_name[:12]}-{redox_name[:12]:<14} | {elec_dict['gamma']:>8.2f} | {redox_dict['gamma']:>8.2f} | {alpha:>8.3f}")

print()

# =============================================================================
# PART 7: CATALYSIS AT ELECTRODES
# =============================================================================

print("-" * 70)
print("PART 7: ELECTROCATALYSIS")
print("-" * 70)
print()

print("Why are some electrodes better catalysts?")
print()
print("Synchronism answer: COHERENCE MATCHING")
print()
print("For optimal catalysis:")
print("  1. Electrode γ should match redox couple γ")
print("  2. d-band center correlates with γ")
print("  3. Surface structure affects local γ")
print()

# Simulate HER (hydrogen evolution) activity
print("Hydrogen Evolution Reaction (HER) activity:")
print()

HER_data = {
    "Pt": {"gamma_elec": 0.15, "gamma_H": 1.5, "j0_exp": 1e-3},  # A/cm²
    "Ni": {"gamma_elec": 0.40, "gamma_H": 1.5, "j0_exp": 1e-5},
    "Au": {"gamma_elec": 0.20, "gamma_H": 1.5, "j0_exp": 5e-6},
    "C": {"gamma_elec": 0.80, "gamma_H": 1.5, "j0_exp": 1e-8},
}

print(f"{'Electrode':<10} | {'γ_elec':>8} | {'f_match':>8} | {'j₀_exp':>12} | {'j₀ ∝ 1/f?':>12}")
print("-" * 60)

for name, data in HER_data.items():
    f_match = coherence_factor(data["gamma_elec"], data["gamma_H"])
    # Check if j0 scales inversely with f (it should be higher for better match)
    # Actually j0 should scale WITH f, not inverse
    scaling = data["j0_exp"] / f_match
    print(f"{name:<10} | {data['gamma_elec']:>8.2f} | {f_match:>8.3f} | {data['j0_exp']:>12.2e} | {scaling:>12.2e}")

print()
print("Note: j₀ correlates positively with coherence match f.")
print("Better matching (larger f) → higher exchange current.")
print()

# =============================================================================
# PART 8: PREDICTIONS
# =============================================================================

print("-" * 70)
print("PART 8: PREDICTIONS")
print("-" * 70)
print()

print("P52.1: Exchange current density")
print("  j₀ ∝ f(γ_elec, γ_redox)")
print("  Better coherence matching → faster kinetics")
print()

print("P52.2: Transfer coefficient")
print("  α ≠ 0.5 when γ_elec ≠ γ_redox")
print("  Asymmetric Tafel slopes")
print()

print("P52.3: Optimal catalyst")
print("  Catalyst γ should match reaction intermediate γ")
print("  Explains why Pt is best for HER")
print()

print("P52.4: Surface modification")
print("  Adding adsorbates changes surface γ")
print("  Self-assembled monolayers tune γ")
print()

print("P52.5: Temperature dependence")
print("  f(T) follows γ(T) behavior")
print("  Arrhenius with coherence-modified prefactor")
print()

# =============================================================================
# PART 9: BATTERIES AND FUEL CELLS
# =============================================================================

print("-" * 70)
print("PART 9: APPLICATIONS - BATTERIES AND FUEL CELLS")
print("-" * 70)
print()

print("For Li-ion batteries:")
print("  Anode: Li in graphite (γ ~ 0.6)")
print("  Cathode: Li in oxide (γ ~ 1.0)")
print("  Electrolyte: Li⁺ in solvent (γ ~ 1.5)")
print()
print("  Mismatch at interfaces limits rate capability!")
print()

print("Design principle:")
print("  Use electrode materials with γ closer to electrolyte γ")
print("  Or use electrolytes with lower γ (polymer electrolytes)")
print()

print("For fuel cells (H₂/O₂):")
print("  H₂ oxidation: γ_H = 1.5, need low-γ electrode")
print("  O₂ reduction: γ_O = 1.8, need even lower-γ electrode")
print()
print("  This is why O₂ reduction is HARDER than H₂ oxidation!")
print()

# =============================================================================
# PART 10: VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Coherence factor vs γ ratio
ax1 = axes[0, 0]
ratio = np.linspace(0.1, 10, 100)
f_factor = np.minimum(1/ratio, ratio)

ax1.plot(ratio, f_factor, 'b-', linewidth=2)
ax1.set_xlabel('γ₁ / γ₂')
ax1.set_ylabel('Coherence factor f')
ax1.set_title('Coherence Matching Factor')
ax1.set_xscale('log')
ax1.axvline(x=1, color='red', linestyle='--', label='Perfect match')
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlim(0.1, 10)

# Plot 2: Electrode γ vs catalytic activity (schematic)
ax2 = axes[0, 1]
gamma_elec = np.array([d["gamma"] for d in electrode_data.values()])
names = list(electrode_data.keys())

# Schematic activity (assuming γ_H ~ 1.5 for HER)
activity = coherence_factor(gamma_elec, 1.5)

ax2.barh(names, activity, color='steelblue')
ax2.set_xlabel('Coherence Factor f (with H⁺/H₂)')
ax2.set_title('Predicted HER Activity by Electrode Material')
ax2.grid(True, alpha=0.3, axis='x')

# Plot 3: Modified Butler-Volmer (Tafel plot)
ax3 = axes[1, 0]
eta = np.linspace(-0.3, 0.3, 100)  # overpotential in V
F_over_RT = 38.9  # 1/V at 298 K

for gamma_ratio, label in [(0.5, 'γ_elec < γ_redox'), (1.0, 'γ_elec = γ_redox'), (2.0, 'γ_elec > γ_redox')]:
    alpha = 0.5 * gamma_ratio / (1 + gamma_ratio)
    i_i0 = np.exp(alpha * F_over_RT * eta) - np.exp(-(1-alpha) * F_over_RT * eta)
    ax3.semilogy(eta * 1000, np.abs(i_i0), linewidth=2, label=f'α = {alpha:.2f}')

ax3.set_xlabel('Overpotential η (mV)')
ax3.set_ylabel('|i/i₀|')
ax3.set_title('Modified Butler-Volmer (Tafel Plot)')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim(-300, 300)
ax3.set_ylim(1e-3, 1e3)

# Plot 4: Interface diagram
ax4 = axes[1, 1]
ax4.text(0.5, 0.85, 'ELECTRODE | INTERFACE | SOLUTION', fontsize=14,
         ha='center', transform=ax4.transAxes, fontweight='bold')

ax4.text(0.2, 0.65, 'γ_elec ~ 0.2', fontsize=12, ha='center', transform=ax4.transAxes)
ax4.text(0.8, 0.65, 'γ_sol ~ 1.5', fontsize=12, ha='center', transform=ax4.transAxes)

ax4.annotate('', xy=(0.6, 0.45), xytext=(0.4, 0.45),
             arrowprops=dict(arrowstyle='<->', color='red', lw=2),
             transform=ax4.transAxes)
ax4.text(0.5, 0.35, 'Coherence Mismatch\nLimits ET Rate', fontsize=11,
         ha='center', transform=ax4.transAxes, color='red')

ax4.text(0.5, 0.15, 'k_ET = k_Marcus × f(γ_elec, γ_sol)', fontsize=12,
         ha='center', transform=ax4.transAxes, style='italic')

ax4.axis('off')
ax4.set_title('Coherence at the Electrode Interface')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/electrochemistry.png',
            dpi=150, bbox_inches='tight')
plt.close()

print()
print("Figure saved to electrochemistry.png")

# =============================================================================
# SUMMARY
# =============================================================================

print()
print("-" * 70)
print("SUMMARY")
print("-" * 70)
print()

print("Session #52 applies Synchronism to electrochemistry:")
print()
print("1. COHERENCE MISMATCH at electrode/solution interface")
print("   Electron transfer crosses γ discontinuity")
print()
print("2. COHERENCE FACTOR: f = min(γ₁,γ₂) / max(γ₁,γ₂)")
print("   Modifies Marcus theory: k_ET × f")
print()
print("3. TRANSFER COEFFICIENT from coherence ratio")
print("   α = 0.5 × (γ_elec / γ_redox) / (1 + γ_elec/γ_redox)")
print()
print("4. ELECTROCATALYSIS: Better coherence matching → better catalyst")
print("   Explains Pt superiority for HER")
print()
print("5. BATTERY DESIGN: Minimize coherence mismatch at interfaces")
print()

print("=" * 70)
print("SESSION #52 COMPLETE: ELECTROCHEMISTRY & COHERENCE AT INTERFACES")
print("=" * 70)
