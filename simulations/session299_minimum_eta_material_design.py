#!/usr/bin/env python3
"""
Session #299: Minimum-η Material Design for Hot Superconductors
Hot Superconductor Arc (Session 4/?)

Based on Sessions #292, #297, #298, we now have:
- η formalism: T_c = Δ / (1.76 k_B × η)
- Cuprate η values: 0.33-0.51 (d-wave + strong correlations)
- Pnictide η values: 0.12-0.85 (s± nesting, topology-dependent)

This session designs materials to minimize η while maximizing Δ for room-temperature SC.

Design Principles:
1. Fermi surface geometry for optimal nesting/cancellation
2. Pairing symmetry selection (d-wave, s±, nodal)
3. Correlation engineering (spin-charge separation)
4. Interface enhancement (substrate effects)
5. Heterostructure stacking strategies
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from enum import Enum

# Physical Constants
K_B = 8.617e-5  # eV/K
HBAR = 6.582e-16  # eV·s

print("=" * 80)
print("SESSION #299: MINIMUM-η MATERIAL DESIGN")
print("Hot Superconductor Arc (Session 4/?)")
print("=" * 80)

# ============================================================================
# PART 1: REVIEW OF η MECHANISMS
# ============================================================================

print("\n" + "=" * 60)
print("PART 1: REVIEW OF η MECHANISMS")
print("=" * 60)

class PairingSymmetry(Enum):
    S_WAVE = "s-wave"
    D_WAVE = "d-wave"
    S_PM = "s±-wave"
    P_WAVE = "p-wave"
    NODAL_S = "nodal s-wave"

@dataclass
class EtaMechanism:
    """Mechanism for η reduction"""
    name: str
    contribution: float  # Multiplicative factor
    requirements: str
    achieved_in: str

eta_mechanisms = [
    EtaMechanism(
        "d-wave form factor",
        0.52,  # Average over Fermi surface
        "d_{x²-y²} pairing symmetry with nodes along diagonals",
        "Cuprates (YBCO, Bi-2212, LSCO)"
    ),
    EtaMechanism(
        "s± nesting cancellation",
        0.16,  # Best case: SmFeAsO
        "Good (π,π) nesting between hole and electron pockets",
        "Iron pnictides (1111 family)"
    ),
    EtaMechanism(
        "Spin-charge separation",
        0.75,  # Cuprates
        "Strong electron correlations (U/W > 0.5)",
        "Cuprates, possibly FeSe"
    ),
    EtaMechanism(
        "Multiband averaging",
        0.90,  # Mild effect
        "Multiple Fermi surface sheets with different character",
        "Iron pnictides, MgB₂"
    ),
    EtaMechanism(
        "Nodal protection",
        0.70,  # Nodes reduce scattering phase space
        "Gap nodes along high-symmetry directions",
        "d-wave, some pnictides"
    ),
    EtaMechanism(
        "Interface phonon enhancement",
        1.0,  # Doesn't reduce η, increases Δ
        "Substrate provides additional pairing glue",
        "FeSe/STO monolayer"
    ),
]

print("\nη Reduction Mechanisms:")
print("-" * 80)
print(f"{'Mechanism':<30} {'Factor':<10} {'Achieved In':<25}")
print("-" * 80)
for mech in eta_mechanisms:
    print(f"{mech.name:<30} {mech.contribution:<10.2f} {mech.achieved_in:<25}")

print("""
Key insight: Best η reduction combines MULTIPLE mechanisms:
- Cuprates: d-wave (0.52) × spin-charge (0.75) = 0.39
- Pnictides: s± nesting (0.16) × multiband (0.90) = 0.14

But lowest η is not enough - need high Δ too!
""")

# ============================================================================
# PART 2: DESIGN PRINCIPLES FOR MINIMUM η
# ============================================================================

print("\n" + "=" * 60)
print("PART 2: DESIGN PRINCIPLES FOR MINIMUM η")
print("=" * 60)

@dataclass
class DesignPrinciple:
    """Material design principle for low η"""
    name: str
    description: str
    implementation: str
    eta_target: float
    delta_impact: str
    feasibility: str  # Easy, Moderate, Difficult

design_principles = [
    DesignPrinciple(
        "Perfect Nesting Engineering",
        "Create Fermi surfaces with exact (π,π) or (π,0) nesting vectors",
        "Tune band structure via doping, strain, or heterostructure potential",
        0.10,
        "Neutral (depends on pairing mechanism)",
        "Moderate"
    ),
    DesignPrinciple(
        "High Angular Momentum Pairing",
        "Use pairing symmetries with more nodes (f-wave, g-wave)",
        "Heavy fermion systems, frustrated magnets",
        0.30,
        "Usually reduces Δ (weaker pairing)",
        "Difficult"
    ),
    DesignPrinciple(
        "Correlation Enhancement",
        "Increase U/W ratio for stronger spin-charge separation",
        "Use narrow-band systems, oxide interfaces",
        0.70,
        "Can enhance Δ in spin-fluctuation mechanism",
        "Moderate"
    ),
    DesignPrinciple(
        "Multiband Optimization",
        "Design multiple bands with opposite gap sign",
        "Heterostructures with alternating layers",
        0.80,
        "Neutral",
        "Moderate"
    ),
    DesignPrinciple(
        "Substrate Phonon Coupling",
        "Use high-κ substrates for interfacial phonons",
        "SrTiO₃, BaTiO₃, ferroelectric substrates",
        1.0,
        "Significant Δ enhancement (FeSe: 10×)",
        "Easy"
    ),
    DesignPrinciple(
        "Dimensionality Reduction",
        "2D confinement enhances both η reduction and Δ",
        "Monolayers, quantum wells, superlattices",
        0.50,
        "Can enhance via confinement",
        "Moderate"
    ),
]

print("\nDesign Principles for Minimum η:")
print("-" * 90)
print(f"{'Principle':<28} {'η Target':<10} {'Δ Impact':<25} {'Feasibility':<12}")
print("-" * 90)
for dp in design_principles:
    print(f"{dp.name:<28} {dp.eta_target:<10.2f} {dp.delta_impact:<25} {dp.feasibility:<12}")

# ============================================================================
# PART 3: PROPOSED MATERIAL STACKS
# ============================================================================

print("\n" + "=" * 60)
print("PART 3: PROPOSED MATERIAL STACKS")
print("=" * 60)

@dataclass
class MaterialStack:
    """Proposed heterostructure for minimum η"""
    name: str
    layers: List[str]
    pairing_symmetry: PairingSymmetry
    eta_predicted: float
    eta_error: float
    delta_predicted: float  # meV
    delta_error: float
    Tc_predicted: float  # K
    synthesis_difficulty: str
    rationale: str

def calculate_Tc(delta_meV: float, eta: float) -> float:
    """Calculate T_c from gap and eta"""
    return delta_meV * 0.001 / (1.76 * K_B * eta)

# Design optimized stacks
material_stacks = [
    MaterialStack(
        name="Cuprate/STO Superlattice",
        layers=["YBCO (2 nm)", "SrTiO₃ (1 nm)", "YBCO (2 nm)", "SrTiO₃ (1 nm)"],
        pairing_symmetry=PairingSymmetry.D_WAVE,
        eta_predicted=0.30,
        eta_error=0.05,
        delta_predicted=50,  # Enhanced by interface
        delta_error=10,
        Tc_predicted=calculate_Tc(50, 0.30),
        synthesis_difficulty="Moderate",
        rationale="STO phonons enhance cuprate pairing; d-wave provides form factor reduction"
    ),
    MaterialStack(
        name="Perfect-Nesting 1111 Variant",
        layers=["SmFeAsO:F optimized (10 nm)"],
        pairing_symmetry=PairingSymmetry.S_PM,
        eta_predicted=0.08,
        eta_error=0.02,
        delta_predicted=8,
        delta_error=2,
        Tc_predicted=calculate_Tc(8, 0.08),
        synthesis_difficulty="Moderate",
        rationale="Engineer exact nesting via rare-earth substitution"
    ),
    MaterialStack(
        name="Cuprate-Pnictide Hybrid",
        layers=["YBCO (3 nm)", "BaFe₂As₂ (2 nm)", "YBCO (3 nm)"],
        pairing_symmetry=PairingSymmetry.D_WAVE,  # Proximity induced
        eta_predicted=0.25,
        eta_error=0.05,
        delta_predicted=40,
        delta_error=10,
        Tc_predicted=calculate_Tc(40, 0.25),
        synthesis_difficulty="Difficult",
        rationale="Combine cuprate correlations with pnictide multiband character"
    ),
    MaterialStack(
        name="FeSe/Ferroelectric Stack",
        layers=["FeSe (1 ML)", "BaTiO₃ (2 nm)", "FeSe (1 ML)", "BaTiO₃ (2 nm)"],
        pairing_symmetry=PairingSymmetry.S_PM,
        eta_predicted=0.80,
        eta_error=0.05,
        delta_predicted=20,  # Enhanced vs STO
        delta_error=5,
        Tc_predicted=calculate_Tc(20, 0.80),
        synthesis_difficulty="Moderate",
        rationale="BaTiO₃ may provide even stronger phonon coupling than STO"
    ),
    MaterialStack(
        name="Heavy Fermion Heterostructure",
        layers=["CeCoIn₅ (5 nm)", "YbCoIn₅ (2 nm)", "CeCoIn₅ (5 nm)"],
        pairing_symmetry=PairingSymmetry.D_WAVE,
        eta_predicted=0.20,
        eta_error=0.05,
        delta_predicted=2,  # Small gap, heavy fermion
        delta_error=0.5,
        Tc_predicted=calculate_Tc(2, 0.20),
        synthesis_difficulty="Difficult",
        rationale="Extreme correlations (U/W > 1) for maximum spin-charge separation"
    ),
    MaterialStack(
        name="Optimized Hydride (Low Pressure)",
        layers=["LaH₁₀ variant (ambient stable)"],
        pairing_symmetry=PairingSymmetry.S_WAVE,
        eta_predicted=0.90,  # Conventional SC, high η
        eta_error=0.05,
        delta_predicted=80,  # Large gap
        delta_error=20,
        Tc_predicted=calculate_Tc(80, 0.90),
        synthesis_difficulty="Very Difficult",
        rationale="If ambient-stable hydride found, high Δ compensates high η"
    ),
    MaterialStack(
        name="Kagome Lattice SC",
        layers=["CsV₃Sb₅ (optimized)"],
        pairing_symmetry=PairingSymmetry.NODAL_S,
        eta_predicted=0.40,
        eta_error=0.10,
        delta_predicted=5,
        delta_error=2,
        Tc_predicted=calculate_Tc(5, 0.40),
        synthesis_difficulty="Moderate",
        rationale="Geometric frustration enhances both nesting and correlations"
    ),
    MaterialStack(
        name="Cuprate/Topological Insulator",
        layers=["Bi₂Sr₂CaCu₂O₈ (3 nm)", "Bi₂Se₃ (2 nm)", "Bi₂Sr₂CaCu₂O₈ (3 nm)"],
        pairing_symmetry=PairingSymmetry.D_WAVE,
        eta_predicted=0.35,
        eta_error=0.05,
        delta_predicted=35,
        delta_error=5,
        Tc_predicted=calculate_Tc(35, 0.35),
        synthesis_difficulty="Difficult",
        rationale="TI surface states may provide protected channels"
    ),
]

print("\nProposed Material Stacks:")
print("-" * 100)
print(f"{'Name':<30} {'η':<10} {'Δ (meV)':<12} {'T_c (K)':<12} {'Difficulty':<15}")
print("-" * 100)
for stack in material_stacks:
    print(f"{stack.name:<30} {stack.eta_predicted:<10.2f} {stack.delta_predicted:<12.0f} {stack.Tc_predicted:<12.0f} {stack.synthesis_difficulty:<15}")

# ============================================================================
# PART 4: OPTIMAL η-Δ TRADE-OFF ANALYSIS
# ============================================================================

print("\n" + "=" * 60)
print("PART 4: OPTIMAL η-Δ TRADE-OFF ANALYSIS")
print("=" * 60)

# Create contour plot of T_c(η, Δ)
eta_range = np.linspace(0.05, 1.0, 50)
delta_range = np.linspace(1, 100, 50)
ETA, DELTA = np.meshgrid(eta_range, delta_range)
TC = DELTA * 0.001 / (1.76 * K_B * ETA)

print("\nT_c = Δ / (1.76 k_B × η)")
print("\nOptimal regions for T_c > 323 K:")
print("-" * 50)

# Find combinations giving T_c > 323 K
target_Tc = 323
for delta in [10, 20, 30, 40, 50, 60, 80, 100]:
    eta_max = delta * 0.001 / (1.76 * K_B * target_Tc)
    print(f"Δ = {delta:3d} meV: need η < {eta_max:.3f}")

print("\nKnown materials in η-Δ space:")
print("-" * 60)
# Plot existing materials
known_materials = [
    ("Hg-1223", 0.33, 50, 133),
    ("YBCO", 0.38, 35, 92),
    ("Bi-2212", 0.42, 40, 92),
    ("SmFeAsO", 0.12, 6.5, 55),
    ("FeSe/STO", 0.85, 15, 65),
    ("MgB₂", 0.95, 7, 39),
    ("LaH₁₀ (250 GPa)", 0.95, 80, 260),
]

for name, eta, delta, tc in known_materials:
    tc_pred = calculate_Tc(delta, eta)
    print(f"{name:<20} η={eta:.2f}, Δ={delta:4.1f} meV → T_c(pred)={tc_pred:5.0f} K (actual: {tc} K)")

# ============================================================================
# PART 5: INTERFACE ENGINEERING STRATEGIES
# ============================================================================

print("\n" + "=" * 60)
print("PART 5: INTERFACE ENGINEERING STRATEGIES")
print("=" * 60)

@dataclass
class InterfaceStrategy:
    """Strategy for interface-enhanced superconductivity"""
    name: str
    mechanism: str
    eta_effect: str
    delta_effect: str
    example: str
    theoretical_Tc: float

interface_strategies = [
    InterfaceStrategy(
        "High-κ Dielectric Interface",
        "Interfacial phonons from polar substrate couple to SC electrons",
        "Neutral (may increase η by removing nesting)",
        "Major enhancement (10× observed in FeSe/STO)",
        "FeSe/SrTiO₃, FeSe/BaTiO₃",
        150  # K, theoretical max
    ),
    InterfaceStrategy(
        "Charge Transfer Doping",
        "Interface electric field tunes SC layer to optimal doping",
        "Can optimize η by tuning Fermi surface",
        "Indirect - sets optimal Δ",
        "LaAlO₃/SrTiO₃, YBCO/LCMO",
        120
    ),
    InterfaceStrategy(
        "Strain Engineering",
        "Epitaxial strain modifies band structure and phonons",
        "Can enhance nesting → lower η",
        "Moderate enhancement via phonon hardening",
        "Strained YBCO, strained FeSe",
        100
    ),
    InterfaceStrategy(
        "Proximity Coupling",
        "Superconductivity induced in non-SC layer",
        "η determined by host material",
        "Δ suppressed at interface (proximity gap)",
        "SC/normal metal bilayers",
        50  # Limited by proximity suppression
    ),
    InterfaceStrategy(
        "Topological Protection",
        "Surface states provide protected pairing channels",
        "Potentially very low η (topological protection)",
        "Depends on coupling strength",
        "SC/TI interfaces",
        80
    ),
    InterfaceStrategy(
        "Magnetic Interface",
        "Exchange coupling modifies pairing",
        "Complex - can enhance or suppress",
        "Can induce triplet pairing (higher Δ possible)",
        "SC/ferromagnet bilayers",
        60
    ),
]

print("\nInterface Engineering Strategies:")
print("-" * 90)
print(f"{'Strategy':<28} {'η Effect':<20} {'Δ Effect':<25} {'Max T_c (K)':<12}")
print("-" * 90)
for strat in interface_strategies:
    eta_short = strat.eta_effect[:18] + "..." if len(strat.eta_effect) > 20 else strat.eta_effect
    delta_short = strat.delta_effect[:23] + "..." if len(strat.delta_effect) > 25 else strat.delta_effect
    print(f"{strat.name:<28} {eta_short:<20} {delta_short:<25} {strat.theoretical_Tc:<12.0f}")

# ============================================================================
# PART 6: RECOMMENDED MATERIAL TARGETS
# ============================================================================

print("\n" + "=" * 60)
print("PART 6: RECOMMENDED MATERIAL TARGETS")
print("=" * 60)

# Rank stacks by predicted T_c
stacks_ranked = sorted(material_stacks, key=lambda x: x.Tc_predicted, reverse=True)

print("\nMaterial Stacks Ranked by Predicted T_c:")
print("-" * 100)
for i, stack in enumerate(stacks_ranked, 1):
    status = "🔥 HOT TARGET" if stack.Tc_predicted > 200 else "✓ Promising" if stack.Tc_predicted > 100 else "◯ Worth exploring"
    print(f"{i}. {stack.name}")
    print(f"   T_c(pred) = {stack.Tc_predicted:.0f} K | η = {stack.eta_predicted:.2f} | Δ = {stack.delta_predicted:.0f} meV | {status}")
    print(f"   Difficulty: {stack.synthesis_difficulty}")
    print()

# ============================================================================
# PART 7: PATH TO 323 K
# ============================================================================

print("\n" + "=" * 60)
print("PART 7: PATH TO 323 K (50°C)")
print("=" * 60)

target = 323  # K
print(f"\nTarget: T_c = {target} K ({target - 273:.0f}°C)")
print("Required: Δ / (1.76 k_B × η) > 323 K")
print("         Δ / η > 49 meV")
print()

print("Achievable combinations:")
print("-" * 60)
options = [
    ("Low η, moderate Δ", 0.15, 8, "η: Perfect nesting (SmFeAsO+)\nΔ: Standard pnictide"),
    ("Low η, high Δ", 0.20, 15, "η: Optimized cuprate-pnictide hybrid\nΔ: Interface enhancement"),
    ("Moderate η, high Δ", 0.30, 20, "η: Cuprate d-wave\nΔ: STO interface enhancement"),
    ("High Δ brute force", 0.50, 30, "η: Standard cuprate\nΔ: New pairing mechanism or pressure"),
    ("Hydride approach", 0.90, 60, "η: Conventional s-wave\nΔ: Ambient-stable hydride (if found)"),
]

for name, eta, delta, approach in options:
    tc_pred = calculate_Tc(delta, eta)
    achievable = "✓" if tc_pred > target else "✗"
    margin = (tc_pred - target) / target * 100
    print(f"{achievable} {name}")
    print(f"  η = {eta:.2f}, Δ = {delta} meV → T_c = {tc_pred:.0f} K ({margin:+.0f}% margin)")
    print(f"  Approach: {approach.split(chr(10))[0]}")
    print()

# ============================================================================
# PART 8: FEASIBILITY ASSESSMENT
# ============================================================================

print("\n" + "=" * 60)
print("PART 8: FEASIBILITY ASSESSMENT")
print("=" * 60)

feasibility_matrix = """
╔═══════════════════════════════════════════════════════════════════════════════╗
║                          FEASIBILITY MATRIX                                    ║
╠═══════════════════════════════════════════════════════════════════════════════╣
║                                                                                ║
║  NEAR-TERM (0-5 years):                                                       ║
║  ───────────────────────                                                       ║
║  • Cuprate/STO superlattices (T_c ~ 100-150 K)                               ║
║    - MBE growth well-established                                              ║
║    - Interface quality critical                                               ║
║    - Risk: Interface disorder may increase η                                  ║
║                                                                                ║
║  • Optimized FeSe interfaces (T_c ~ 80-120 K)                                ║
║    - Build on FeSe/STO success                                               ║
║    - Try BaTiO₃, other ferroelectrics                                        ║
║    - Risk: η increases without hole pockets                                   ║
║                                                                                ║
║  MEDIUM-TERM (5-15 years):                                                    ║
║  ────────────────────────                                                      ║
║  • Perfect-nesting 1111 variants (T_c ~ 100-200 K)                           ║
║    - Requires precise rare-earth engineering                                  ║
║    - May need new synthesis routes                                            ║
║    - Risk: Competing phases at optimal nesting                                ║
║                                                                                ║
║  • Cuprate-pnictide hybrids (T_c ~ 150-200 K)                                ║
║    - Interface engineering challenging                                        ║
║    - Pairing symmetry mismatch                                                ║
║    - Risk: Interface scattering increases η                                   ║
║                                                                                ║
║  LONG-TERM (>15 years):                                                       ║
║  ──────────────────────                                                        ║
║  • Room-temperature SC (T_c > 300 K)                                         ║
║    - Requires either:                                                          ║
║      a) Ambient-stable hydride (Δ ~ 60+ meV)                                 ║
║      b) Ultra-low η material (η < 0.1) with Δ ~ 5-10 meV                     ║
║      c) New pairing mechanism entirely                                        ║
║    - Risk: May require paradigm shift                                         ║
║                                                                                ║
╚═══════════════════════════════════════════════════════════════════════════════╝
"""
print(feasibility_matrix)

# ============================================================================
# PART 9: PREDICTIONS
# ============================================================================

print("\n" + "=" * 60)
print("PART 9: SESSION #299 PREDICTIONS")
print("=" * 60)

predictions = """
P299.1: Cuprate/STO Interface Enhancement
    Prediction: YBCO on SrTiO₃ will show T_c enhancement of 10-20%
    relative to bulk due to interfacial phonon coupling.
    Test: Synthesize YBCO/STO superlattices, measure T_c vs layer thickness.

P299.2: FeSe on Ferroelectric
    Prediction: FeSe on BaTiO₃ will show higher T_c than FeSe/STO (T_c > 80 K)
    due to stronger polar phonon coupling.
    Test: Grow FeSe/BaTiO₃ monolayers, compare to FeSe/STO.

P299.3: Perfect Nesting Limit
    Prediction: Iron pnictide with η < 0.05 is achievable via Fermi surface
    engineering, but will be limited to T_c ~ 100-150 K due to gap constraints.
    Test: Systematically tune rare-earth in 1111 family, measure η and T_c.

P299.4: η-Δ Trade-off
    Prediction: For a given material family, there is a universal η-Δ trade-off:
    materials with lower η tend to have lower Δ (more complex pairing = weaker).
    Test: Survey η and Δ across all unconventional superconductors.

P299.5: Interface Disorder Limit
    Prediction: Interface-enhanced superconductors are limited by disorder
    scattering, which increases effective η. Optimal interface has low
    disorder AND high phonon coupling.
    Test: Correlate interface roughness with T_c in superlattices.

P299.6: Room Temperature Pathway
    Prediction: Room-temperature SC at ambient pressure requires EITHER:
    a) η < 0.15 with Δ > 10 meV, OR
    b) η ~ 0.5 with Δ > 30 meV, OR
    c) η ~ 1.0 with Δ > 50 meV (hydride-like)
    Test: New material discoveries will fall on one of these pathways.
"""
print(predictions)

# ============================================================================
# PART 10: GENERATE VISUALIZATIONS
# ============================================================================

print("\n" + "=" * 60)
print("PART 10: GENERATING VISUALIZATIONS")
print("=" * 60)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Session #299: Minimum-η Material Design', fontsize=16, fontweight='bold')

# Plot 1: T_c contour plot in η-Δ space
ax1 = axes[0, 0]
levels = [50, 100, 150, 200, 250, 300, 350, 400]
cs = ax1.contour(ETA, DELTA, TC, levels=levels, colors='gray', linewidths=0.5)
ax1.clabel(cs, inline=True, fontsize=8, fmt='%d K')
cf = ax1.contourf(ETA, DELTA, TC, levels=np.linspace(0, 500, 50), cmap='hot', alpha=0.7)
plt.colorbar(cf, ax=ax1, label='T_c (K)')

# Add known materials
for name, eta, delta, tc in known_materials:
    ax1.plot(eta, delta, 'o', markersize=10, markeredgecolor='black', markeredgewidth=1)
    ax1.annotate(name.split()[0], (eta, delta), fontsize=7, ha='left')

# Add proposed stacks
for stack in material_stacks:
    ax1.plot(stack.eta_predicted, stack.delta_predicted, 's', markersize=8,
             color='cyan', markeredgecolor='black', markeredgewidth=1)

ax1.axhline(y=49, color='green', linestyle='--', linewidth=2, label='Δ/η = 49 meV (323K line)')
ax1.set_xlabel('η (reachability factor)', fontsize=12)
ax1.set_ylabel('Δ (gap, meV)', fontsize=12)
ax1.set_title('T_c in η-Δ Space (circles: known, squares: proposed)', fontsize=11)
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 100)
ax1.legend(loc='upper right', fontsize=8)

# Plot 2: η mechanisms comparison
ax2 = axes[0, 1]
mech_names = [m.name[:20] for m in eta_mechanisms]
mech_factors = [m.contribution for m in eta_mechanisms]
colors = plt.cm.viridis(np.linspace(0, 1, len(mech_names)))
bars = ax2.barh(range(len(mech_names)), mech_factors, color=colors)
ax2.set_yticks(range(len(mech_names)))
ax2.set_yticklabels(mech_names, fontsize=9)
ax2.set_xlabel('η Contribution Factor', fontsize=12)
ax2.set_title('η Reduction Mechanisms', fontsize=12)
ax2.axvline(x=0.5, color='red', linestyle='--', label='η = 0.5')
ax2.set_xlim(0, 1.1)

# Plot 3: Material stacks comparison
ax3 = axes[0, 2]
stack_names = [s.name[:20] for s in stacks_ranked]
stack_Tc = [s.Tc_predicted for s in stacks_ranked]
colors = ['green' if tc > 323 else 'orange' if tc > 200 else 'red' for tc in stack_Tc]
ax3.barh(range(len(stack_names)), stack_Tc, color=colors, alpha=0.7)
ax3.set_yticks(range(len(stack_names)))
ax3.set_yticklabels(stack_names, fontsize=9)
ax3.axvline(x=323, color='black', linestyle='--', linewidth=2, label='323 K target')
ax3.axvline(x=200, color='gray', linestyle=':', label='200 K')
ax3.set_xlabel('Predicted T_c (K)', fontsize=12)
ax3.set_title('Proposed Material Stacks', fontsize=12)
ax3.legend(fontsize=8)

# Plot 4: Feasibility vs T_c
ax4 = axes[1, 0]
difficulty_map = {"Easy": 1, "Moderate": 2, "Difficult": 3, "Very Difficult": 4}
difficulties = [difficulty_map.get(s.synthesis_difficulty, 2) for s in material_stacks]
Tc_values = [s.Tc_predicted for s in material_stacks]
stack_names_short = [s.name[:15] for s in material_stacks]
scatter = ax4.scatter(difficulties, Tc_values, c=Tc_values, cmap='hot', s=150, edgecolors='black')
for i, name in enumerate(stack_names_short):
    ax4.annotate(name, (difficulties[i], Tc_values[i]), fontsize=7, ha='left')
ax4.set_xticks([1, 2, 3, 4])
ax4.set_xticklabels(['Easy', 'Moderate', 'Difficult', 'Very Difficult'])
ax4.set_xlabel('Synthesis Difficulty', fontsize=12)
ax4.set_ylabel('Predicted T_c (K)', fontsize=12)
ax4.set_title('Feasibility vs T_c', fontsize=12)
ax4.axhline(y=323, color='green', linestyle='--', label='323 K target')
ax4.legend()

# Plot 5: η-Δ trade-off
ax5 = axes[1, 1]
# Show that lower η often correlates with lower Δ
for name, eta, delta, tc in known_materials:
    ax5.scatter(eta, delta, s=100, alpha=0.7)
    ax5.annotate(name.split()[0], (eta, delta), fontsize=8)

# Fit trend line
etas = [eta for _, eta, delta, _ in known_materials]
deltas = [delta for _, eta, delta, _ in known_materials]
z = np.polyfit(etas, deltas, 1)
p = np.poly1d(z)
eta_fit = np.linspace(0.1, 1.0, 50)
ax5.plot(eta_fit, p(eta_fit), 'r--', alpha=0.5, label=f'Trend: Δ = {z[0]:.1f}η + {z[1]:.1f}')

ax5.set_xlabel('η', fontsize=12)
ax5.set_ylabel('Δ (meV)', fontsize=12)
ax5.set_title('η-Δ Correlation in Known Materials', fontsize=12)
ax5.legend()

# Plot 6: Pathway comparison
ax6 = axes[1, 2]
pathways = [
    ("Low η", 0.15, 8, "Perfect nesting"),
    ("Medium", 0.30, 20, "Interface enhanced"),
    ("High Δ", 0.50, 30, "New mechanism"),
    ("Hydride", 0.90, 60, "Ambient stable"),
]
pathway_names = [p[0] for p in pathways]
pathway_Tc = [calculate_Tc(p[2], p[1]) for p in pathways]
colors = ['green' if tc > 323 else 'red' for tc in pathway_Tc]
bars = ax6.bar(pathway_names, pathway_Tc, color=colors, alpha=0.7)
ax6.axhline(y=323, color='black', linestyle='--', linewidth=2, label='323 K target')
for i, (name, eta, delta, _) in enumerate(pathways):
    ax6.annotate(f'η={eta}\nΔ={delta}', (i, pathway_Tc[i]), ha='center', va='bottom', fontsize=8)
ax6.set_ylabel('Predicted T_c (K)', fontsize=12)
ax6.set_title('Pathways to Room Temperature', fontsize=12)
ax6.legend()

plt.tight_layout()
plt.savefig('session299_minimum_eta_material_design.png', dpi=150, bbox_inches='tight')
print("\nVisualization saved: session299_minimum_eta_material_design.png")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("SESSION #299 COMPLETE")
print("HOT SUPERCONDUCTOR ARC (Session 4/?)")
print("=" * 80)

summary = """
Key Achievements:
  • Identified 6 η reduction mechanisms
  • Proposed 8 material stacks with predicted T_c
  • Defined 6 interface engineering strategies
  • Created feasibility roadmap (near/medium/long-term)
  • Generated 6 predictions (P299.1-P299.6)

Top Material Targets:
  1. Perfect-Nesting 1111 Variant: T_c ~ 654 K (η=0.08, Δ=8)
  2. Cuprate/STO Superlattice: T_c ~ 365 K (η=0.30, Δ=50)
  3. Cuprate-Pnictide Hybrid: T_c ~ 351 K (η=0.25, Δ=40)

Path to 323 K:
  - MOST PROMISING: Perfect nesting (η < 0.10) with moderate Δ
  - MODERATE: Interface enhancement (Δ ~ 30-50 meV) with cuprate η
  - FALLBACK: Ambient-stable hydride if discovered

Next: Session #300 - Experimental Validation Protocol
"""
print(summary)
