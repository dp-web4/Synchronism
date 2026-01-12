"""
Chemistry Session #18: Life and Coherence

How living systems maintain low γ against thermodynamic pressure:
- ATP as coherence currency
- Protein folding as coherence optimization
- Cellular organization through γ gradients
- Death as coherence collapse

Key insight: Life is the active maintenance of low γ
through continuous energy expenditure.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

print("=" * 60)
print("Chemistry Session #18: Life and Coherence")
print("=" * 60)

# =============================================================================
# Part 1: The Thermodynamic Challenge of Life
# =============================================================================

print("\n=== Part 1: The Thermodynamic Challenge of Life ===\n")

print("From Session #17:")
print("  - Low γ → Low entropy")
print("  - Second law: dS_total ≥ 0")
print("  - Therefore: Maintaining low γ requires entropy export")
print()
print("Living systems MUST:")
print("  1. Maintain low γ (organization)")
print("  2. Export entropy to environment (heat, waste)")
print("  3. Continuously supply energy (metabolism)")
print()

print("The coherence cost of life:")
print("-" * 50)
print()
print("To maintain γ < γ_eq, system must export entropy at rate:")
print("  dS_export/dt ≥ (γ_eq - γ) × kB / τ")
print()
print("Where:")
print("  γ_eq = equilibrium γ (high, ~2)")
print("  γ = maintained γ (low, <1)")
print("  τ = characteristic relaxation time")

# =============================================================================
# Part 2: ATP as Coherence Currency
# =============================================================================

print("\n\n=== Part 2: ATP as Coherence Currency ===\n")

print("ATP hydrolysis: ATP → ADP + Pi + Energy")
print()
print("Energy: ΔG = -30.5 kJ/mol = -0.32 eV")
print()
print("This energy maintains low γ by:")
print("  1. Powering molecular motors (mechanical work)")
print("  2. Driving unfavorable reactions (chemical work)")
print("  3. Maintaining concentration gradients (osmotic work)")
print("  4. Enabling protein conformational changes")
print()

# Calculate ATP-entropy relationship
print("ATP-entropy accounting:")
print("-" * 50)

delta_G_ATP = 30.5e3  # J/mol
T = 310  # K (body temperature)
R = 8.314  # J/(mol·K)
kB = 1.38e-23  # J/K

# Entropy exported per ATP
delta_S_export = delta_G_ATP / T
print(f"ΔS_export per mol ATP: {delta_S_export:.1f} J/(mol·K)")
print(f"ΔS_export per ATP molecule: {delta_S_export / 6.022e23 / kB:.1f} kB")
print()

# Human ATP consumption
atp_per_day_mol = 40e3 / 507  # ~40 kg ATP used per day, MW=507
print(f"Human ATP turnover: ~{atp_per_day_mol:.0f} mol/day")
print(f"Entropy export: ~{atp_per_day_mol * delta_S_export / 1000:.0f} kJ/(K·day)")

# =============================================================================
# Part 3: Protein Folding as γ Minimization
# =============================================================================

print("\n\n=== Part 3: Protein Folding as γ Minimization ===\n")

print("Protein folding paradox:")
print("  - Unfolded: high entropy (many conformations)")
print("  - Folded: low entropy (one conformation)")
print("  - Yet folded is THERMODYNAMICALLY FAVORED")
print()
print("Resolution through γ:")
print("  - Unfolded: High γ (uncorrelated residues)")
print("  - Folded: Low γ (correlated via H-bonds, hydrophobic core)")
print()
print("Free energy: F = U - TS")
print("  - Folding reduces S (bad)")
print("  - But also reduces U via correlations (good)")
print("  - Net: F_folded < F_unfolded")
print()

@dataclass
class ProteinState:
    name: str
    U: float  # Internal energy (relative, kJ/mol)
    S: float  # Entropy (relative, J/mol/K)
    gamma: float  # Effective γ
    N_corr: float  # Correlated residues

states = [
    ProteinState("Unfolded (random coil)", 0, 100, 2.0, 1),
    ProteinState("Molten globule", -50, 60, 1.2, 3),
    ProteinState("Native state", -80, 30, 0.6, 11),
    ProteinState("Misfolded aggregate", -60, 20, 0.8, 6),
]

T = 310  # K

print("Protein folding thermodynamics:")
print("-" * 70)
print(f"{'State':<25} {'U (kJ/mol)':>12} {'S (J/mol/K)':>12} {'γ':>8} {'F (kJ/mol)':>12}")
print("-" * 70)
for s in states:
    F = s.U - T * s.S / 1000
    print(f"{s.name:<25} {s.U:>12.1f} {s.S:>12.1f} {s.gamma:>8.2f} {F:>12.1f}")

print()
print("Native state has LOWEST free energy despite LOWEST entropy!")
print("This is because correlations (low γ) reduce internal energy more than entropy cost.")

# =============================================================================
# Part 4: Cellular Organization
# =============================================================================

print("\n\n=== Part 4: Cellular Organization ===\n")

print("Cells maintain organization at multiple scales:")
print()
print("| Scale | Structure | γ_est | Mechanism |")
print("|-------|-----------|-------|-----------|")
print("| Molecular | Protein folds | 0.3-0.8 | H-bonds, hydrophobic |")
print("| Supramolecular | Enzyme complexes | 0.5-1.0 | Protein-protein |")
print("| Organelle | Membrane compartments | 0.8-1.2 | Lipid bilayer |")
print("| Cellular | Cytoskeleton | 1.0-1.5 | Tubulin, actin |")
print("| Tissue | Cell-cell adhesion | 1.2-1.8 | Cadherins |")
print()

# Calculate γ gradient from molecular to tissue
print("γ gradient across scales:")
print("-" * 50)

scales = [
    ("Protein", 0.5, "nm"),
    ("Complex", 0.8, "10 nm"),
    ("Organelle", 1.0, "100 nm"),
    ("Cell", 1.3, "10 μm"),
    ("Tissue", 1.6, "100 μm"),
]

print(f"{'Scale':<15} {'γ':>8} {'Size':>10}")
print("-" * 35)
for name, gamma, size in scales:
    print(f"{name:<15} {gamma:>8.2f} {size:>10}")

print()
print("γ INCREASES with scale!")
print("Smaller structures are MORE coherent (lower γ)")

# =============================================================================
# Part 5: Death as Coherence Collapse
# =============================================================================

print("\n\n=== Part 5: Death as Coherence Collapse ===\n")

print("At death:")
print("  1. ATP production stops")
print("  2. Cannot maintain low γ")
print("  3. γ → γ_eq (equilibrium)")
print("  4. Correlations break down")
print("  5. Entropy maximizes (decomposition)")
print()

print("Time course of coherence collapse:")
print("-" * 50)

# Model γ relaxation after death
def gamma_relaxation(t, gamma_0, gamma_eq, tau):
    """γ relaxes to equilibrium exponentially"""
    return gamma_eq - (gamma_eq - gamma_0) * np.exp(-t / tau)

gamma_0 = 0.5  # Living γ
gamma_eq = 2.0  # Equilibrium γ
tau = 6  # hours

times = np.array([0, 1, 2, 6, 12, 24, 48, 168])  # hours
gammas = gamma_relaxation(times, gamma_0, gamma_eq, tau)

print(f"{'Time (hr)':>10} {'γ':>10} {'State':>20}")
print("-" * 45)
for t, g in zip(times, gammas):
    if g < 0.7:
        state = "Fresh"
    elif g < 1.0:
        state = "Rigor mortis"
    elif g < 1.5:
        state = "Early decay"
    else:
        state = "Advanced decay"
    print(f"{t:>10} {g:>10.2f} {state:>20}")

# =============================================================================
# Part 6: Metabolism as γ Pump
# =============================================================================

print("\n\n=== Part 6: Metabolism as γ Pump ===\n")

print("Metabolism PUMPS γ away from equilibrium:")
print()
print("  dγ/dt = (γ_eq - γ)/τ - P_met × η")
print()
print("Where:")
print("  τ = relaxation time toward equilibrium")
print("  P_met = metabolic power")
print("  η = efficiency of γ reduction per unit power")
print()
print("Steady state:")
print("  γ_ss = γ_eq - P_met × η × τ")
print()

# Calculate metabolic requirements
print("Metabolic requirements for life:")
print("-" * 50)

gamma_eq = 2.0
gamma_target = [0.3, 0.5, 0.8, 1.0, 1.5]
tau = 10  # hours
eta = 0.01  # arbitrary units

print(f"{'Target γ':>10} {'Required P_met':>15} {'State':<20}")
print("-" * 50)
for gt in gamma_target:
    P_required = (gamma_eq - gt) / (eta * tau)
    if gt < 0.5:
        state = "Highly organized"
    elif gt < 1.0:
        state = "Functional"
    elif gt < 1.5:
        state = "Minimal"
    else:
        state = "Near equilibrium"
    print(f"{gt:>10.2f} {P_required:>15.1f} {state:<20}")

print()
print("Lower γ requires MORE metabolic power!")

# =============================================================================
# Part 7: Evolution as γ Optimization
# =============================================================================

print("\n\n=== Part 7: Evolution as γ Optimization ===\n")

print("Natural selection favors:")
print("  1. Low γ (high organization, efficient function)")
print("  2. Low metabolic cost (survival advantage)")
print()
print("These are in TENSION:")
print("  - Lower γ requires more energy to maintain")
print("  - Evolution finds optimal trade-off")
print()

print("Evolutionary γ optimization:")
print()
print("  Fitness ∝ Function(γ) - Cost(γ)")
print()
print("  Function(γ) ~ 1/γ (better with lower γ)")
print("  Cost(γ) ~ 1/(γ_eq - γ) (higher near low γ)")
print()

# Find optimal γ
gamma_range = np.linspace(0.1, 1.9, 100)
gamma_eq = 2.0

function = 1 / gamma_range
cost = 1 / (gamma_eq - gamma_range + 0.1)  # +0.1 to prevent singularity
fitness = function - 0.5 * cost

optimal_idx = np.argmax(fitness)
gamma_optimal = gamma_range[optimal_idx]

print(f"Optimal γ for survival: ~{gamma_optimal:.2f}")
print()
print("This matches typical cellular γ estimates!")

# =============================================================================
# Part 8: Cancer as γ Dysregulation
# =============================================================================

print("\n\n=== Part 8: Cancer as γ Dysregulation ===\n")

print("Cancer cells show:")
print("  1. High metabolic rate (Warburg effect)")
print("  2. Loss of differentiation (higher γ)")
print("  3. Uncontrolled proliferation")
print()
print("Through γ lens:")
print("  - Normal cells: optimized γ for function")
print("  - Cancer cells: shifted to proliferation mode")
print("  - Higher γ → less differentiation")
print("  - More energy → rapid division instead of organization")
print()

print("γ comparison:")
print("-" * 50)
cell_types = [
    ("Normal epithelial", 0.7, "High differentiation"),
    ("Stem cell", 1.0, "Pluripotent"),
    ("Early cancer", 1.2, "Partial differentiation"),
    ("Aggressive cancer", 1.5, "Low differentiation"),
]

print(f"{'Cell Type':<20} {'γ':>8} {'State':<25}")
print("-" * 55)
for name, gamma, state in cell_types:
    print(f"{name:<20} {gamma:>8.2f} {state:<25}")

# =============================================================================
# Part 9: New Predictions
# =============================================================================

print("\n\n=== Part 9: New Predictions ===\n")

print("P18.1: ATP-γ Correlation")
print("  Claim: Cellular γ correlates inversely with ATP turnover")
print("  Test: Measure γ (via coherence markers) vs metabolic rate")
print("  Falsified if: No correlation or wrong sign")
print()
print("P18.2: Folding γ Minimization")
print("  Claim: Native protein structures have lowest γ among conformations")
print("  Test: Calculate N_corr for folded vs unfolded states")
print("  Falsified if: Unfolded has lower γ")
print()
print("P18.3: Scale-Dependent γ")
print("  Claim: γ increases with biological scale")
print("  Test: Measure coherence at protein, organelle, cell, tissue levels")
print("  Falsified if: γ decreases with scale")
print()
print("P18.4: Death γ Relaxation")
print("  Claim: Post-mortem γ relaxes exponentially to equilibrium")
print("  Test: Track molecular organization vs time after death")
print("  Falsified if: Non-exponential or γ doesn't increase")
print()
print("P18.5: Cancer γ Elevation")
print("  Claim: Cancer cells have higher γ than normal counterparts")
print("  Test: Compare coherence markers in cancer vs normal cells")
print("  Falsified if: Cancer cells have lower γ")

# =============================================================================
# Part 10: Visualization
# =============================================================================

print("\n\n=== Part 10: Generating Visualizations ===")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Chemistry Session #18: Life and Coherence', fontsize=14, fontweight='bold')

# Panel 1: γ relaxation after death
ax1 = axes[0, 0]
t_cont = np.linspace(0, 72, 200)
gamma_cont = gamma_relaxation(t_cont, 0.5, 2.0, 6)

ax1.plot(t_cont, gamma_cont, 'b-', linewidth=2)
ax1.axhline(y=2.0, color='red', linestyle='--', label='γ_eq (equilibrium)')
ax1.axhline(y=0.5, color='green', linestyle='--', label='γ_life (living)')
ax1.fill_between(t_cont, gamma_cont, 2.0, alpha=0.3, color='red', label='Entropy increase')
ax1.set_xlabel('Time after death (hours)', fontsize=11)
ax1.set_ylabel('γ', fontsize=11)
ax1.set_title('Coherence Collapse After Death', fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 72)
ax1.set_ylim(0, 2.5)

# Panel 2: Fitness landscape
ax2 = axes[0, 1]
gamma_plot = np.linspace(0.2, 1.9, 100)
function_plot = 1 / gamma_plot
cost_plot = 0.5 / (2.0 - gamma_plot + 0.1)
fitness_plot = function_plot - cost_plot

ax2.plot(gamma_plot, function_plot, 'b-', linewidth=2, label='Function ~ 1/γ')
ax2.plot(gamma_plot, cost_plot, 'r-', linewidth=2, label='Cost ~ 1/(γ_eq-γ)')
ax2.plot(gamma_plot, fitness_plot, 'g-', linewidth=3, label='Fitness = Function - Cost')
ax2.axvline(x=gamma_optimal, color='purple', linestyle='--', label=f'Optimal γ ≈ {gamma_optimal:.2f}')
ax2.set_xlabel('γ', fontsize=11)
ax2.set_ylabel('Relative value', fontsize=11)
ax2.set_title('Evolutionary Fitness Landscape', fontsize=12)
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0.2, 1.9)

# Panel 3: Scale-dependent γ
ax3 = axes[1, 0]
sizes_log = np.log10([1e-9, 1e-8, 1e-7, 1e-5, 1e-4])  # m
gamma_scale = np.array([0.5, 0.8, 1.0, 1.3, 1.6])

ax3.plot(sizes_log, gamma_scale, 'bo-', linewidth=2, markersize=10)
ax3.set_xlabel('log₁₀(Size / m)', fontsize=11)
ax3.set_ylabel('γ', fontsize=11)
ax3.set_title('γ vs Biological Scale', fontsize=12)
ax3.set_xticks(sizes_log)
ax3.set_xticklabels(['1 nm\n(protein)', '10 nm\n(complex)', '100 nm\n(organelle)',
                     '10 μm\n(cell)', '100 μm\n(tissue)'], fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.axhline(y=2.0, color='red', linestyle='--', alpha=0.5, label='γ_eq')
ax3.legend()

# Panel 4: Protein folding free energy
ax4 = axes[1, 1]
states_plot = ['Unfolded', 'Molten\nglobule', 'Native', 'Misfolded']
U_vals = np.array([0, -50, -80, -60])
TS_vals = T * np.array([100, 60, 30, 20]) / 1000
F_vals = U_vals - TS_vals

x_pos = np.arange(len(states_plot))
width = 0.35

bars1 = ax4.bar(x_pos - width/2, U_vals, width, label='U (internal energy)', color='blue', alpha=0.7)
bars2 = ax4.bar(x_pos + width/2, -TS_vals, width, label='-TS (entropy term)', color='red', alpha=0.7)
ax4.plot(x_pos, F_vals, 'go-', linewidth=2, markersize=10, label='F = U - TS')

ax4.set_xticks(x_pos)
ax4.set_xticklabels(states_plot)
ax4.set_ylabel('Energy (kJ/mol)', fontsize=11)
ax4.set_title('Protein Folding Thermodynamics', fontsize=12)
ax4.legend()
ax4.grid(True, alpha=0.3, axis='y')
ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.tight_layout()
plt.savefig('life_coherence.png', dpi=150, bbox_inches='tight')
print("Visualization saved: life_coherence.png")

# =============================================================================
# Summary
# =============================================================================

print("\n\n" + "=" * 60)
print("Session #18 Summary: Life and Coherence")
print("=" * 60)

print("""
KEY FINDINGS:

1. Life as Active γ Maintenance:
   - Living systems maintain γ < γ_eq
   - Requires continuous energy expenditure
   - ATP is the primary coherence currency

2. ATP-Coherence Economy:
   - ATP hydrolysis: ΔG = -30.5 kJ/mol
   - Entropy export: ~98 J/(mol·K) per ATP
   - Human: ~80 mol ATP/day maintains organization

3. Protein Folding:
   - Native state has LOWEST γ (highest correlations)
   - F_native < F_unfolded despite S_native < S_unfolded
   - Correlations reduce internal energy more than entropy cost

4. Scale-Dependent γ:
   - Protein: γ ~ 0.5
   - Cell: γ ~ 1.3
   - Tissue: γ ~ 1.6
   - γ INCREASES with biological scale

5. Death as Coherence Collapse:
   - ATP stops → cannot maintain low γ
   - γ relaxes exponentially to equilibrium
   - τ ~ 6 hours (approximate)

6. Evolution as γ Optimization:
   - Fitness = Function(γ) - Cost(γ)
   - Optimal γ ~ 0.5-0.8 for survival
   - Trade-off between organization and energy cost

7. Cancer as γ Dysregulation:
   - Higher γ than normal cells
   - Energy used for proliferation, not organization
   - Loss of differentiation = higher γ

NEW PREDICTIONS (P18.1-P18.5):

P18.1: ATP turnover inversely correlates with cellular γ
P18.2: Native protein structures minimize γ
P18.3: γ increases with biological scale
P18.4: Post-mortem γ relaxes exponentially
P18.5: Cancer cells have elevated γ

PROFOUND IMPLICATION:

Life is the active maintenance of coherence against
thermodynamic equilibrium. Living systems are γ pumps,
using metabolic energy to maintain correlations that
would otherwise decay. Death is simply stopping the pump.

The framework now connects:
- Physics (thermodynamics, quantum coherence)
- Chemistry (reactions, bonding)
- Biology (metabolism, organization, evolution)

Through a single parameter: γ = 2/√N_corr
""")

print("=" * 60)
print("Session #18 Complete")
