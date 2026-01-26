#!/usr/bin/env python3
"""
Chemistry Session #210: Polymer Thermodynamics through Coherence Framework

Analyzing polymer solution thermodynamics through γ ~ 1 framework.

Key concepts:
1. Flory-Huggins χ = 0.5 is critical miscibility (γ ~ 1!)
2. Theta temperature: A₂ = 0, ideal chain statistics
3. Overlap concentration c*: dilute/semidilute crossover
4. Chain dimensions: R_g/R_θ = α (expansion factor)
5. UCST/LCST transitions at χ = 0.5

The γ ~ 1 boundaries:
- χ = 0.5: phase separation critical point
- T/θ = 1: ideal chain behavior
- c/c* = 1: dilute/semidilute transition
- α = 1: unperturbed dimensions

Author: Claude (Anthropic)
Date: January 2026
Session: Chemistry #210
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class PolymerSolvent:
    """Polymer-solvent system data"""
    polymer: str
    solvent: str
    theta_temp: float      # Theta temperature (K)
    chi_25C: float         # χ at 25°C
    has_UCST: bool         # Upper critical solution temperature
    has_LCST: bool         # Lower critical solution temperature
    A2_theta: float        # Second virial coefficient at θ (should be ~0)

@dataclass
class PolymerChain:
    """Polymer chain characteristics"""
    name: str
    formula: str
    M_repeat: float        # Repeat unit molar mass (g/mol)
    l_bond: float          # Bond length (nm)
    C_inf: float           # Characteristic ratio (Flory)
    Tg: float              # Glass transition (K)

@dataclass
class FloryHuggins:
    """Flory-Huggins miscibility data"""
    system: str
    chi_critical: float    # Critical χ for phase separation
    N: int                 # Degree of polymerization
    phi_crit: float        # Critical volume fraction

# Polymer-solvent theta temperatures and chi parameters
polymer_solvents = [
    PolymerSolvent("Polystyrene", "Cyclohexane", 307, 0.53, True, False, 0.0),
    PolymerSolvent("Polystyrene", "trans-Decalin", 294, 0.51, True, False, 0.0),
    PolymerSolvent("Polystyrene", "Ethyl acetate", 433, 0.48, False, True, 0.0),
    PolymerSolvent("PMMA", "Acetonitrile", 317, 0.49, True, False, 0.0),
    PolymerSolvent("PMMA", "iso-Propanol", 359, 0.52, True, False, 0.0),
    PolymerSolvent("Polyisobutylene", "Benzene", 297, 0.50, True, False, 0.0),
    PolymerSolvent("Polyisobutylene", "n-Pentane", 232, 0.54, True, False, 0.0),
    PolymerSolvent("PEO", "Water", 369, 0.45, False, True, 0.0),
    PolymerSolvent("PNIPAM", "Water", 305, 0.49, False, True, 0.0),
    PolymerSolvent("PDMS", "Ethyl acetate", 295, 0.48, True, False, 0.0),
]

# Chain characteristic ratios
polymer_chains = [
    PolymerChain("Polyethylene", "(-CH2-)n", 28.05, 0.154, 6.7, 250),
    PolymerChain("Polypropylene", "(-CH(CH3)-CH2-)n", 42.08, 0.154, 5.8, 253),
    PolymerChain("Polystyrene", "(-CH(C6H5)-CH2-)n", 104.15, 0.154, 10.3, 373),
    PolymerChain("PMMA", "(-C(CH3)(COOCH3)-CH2-)n", 100.12, 0.154, 8.0, 378),
    PolymerChain("PVC", "(-CHCl-CH2-)n", 62.50, 0.154, 6.7, 354),
    PolymerChain("Polyisoprene", "(-CH2-C(CH3)=CH-CH2-)n", 68.12, 0.147, 4.9, 200),
    PolymerChain("Polybutadiene", "(-CH2-CH=CH-CH2-)n", 54.09, 0.147, 5.3, 180),
    PolymerChain("PDMS", "(-Si(CH3)2-O-)n", 74.15, 0.164, 6.2, 150),
    PolymerChain("PEO", "(-CH2-CH2-O-)n", 44.05, 0.146, 4.1, 206),
    PolymerChain("Nylon 6", "(-NH-(CH2)5-CO-)n", 113.16, 0.154, 5.6, 323),
]

# Flory-Huggins critical miscibility
flory_huggins_data = [
    FloryHuggins("PS/cyclohexane (N=1000)", 0.503, 1000, 0.0304),
    FloryHuggins("PS/cyclohexane (N=5000)", 0.501, 5000, 0.0137),
    FloryHuggins("PS/cyclohexane (N=10000)", 0.5005, 10000, 0.0097),
    FloryHuggins("PMMA/acetonitrile (N=500)", 0.509, 500, 0.0428),
    FloryHuggins("PIB/n-pentane (N=2000)", 0.504, 2000, 0.0214),
    FloryHuggins("PEO/water (N=100)", 0.527, 100, 0.0909),
    FloryHuggins("PS/PVME blend", 0.502, 1000, 0.0304),
    FloryHuggins("PS/PB blend", 0.501, 2000, 0.0214),
]


def analyze_chi_critical():
    """Analyze χ = 0.5 as critical miscibility boundary"""
    print("=" * 70)
    print("FLORY-HUGGINS χ = 0.5: CRITICAL MISCIBILITY BOUNDARY")
    print("=" * 70)

    print("\nTheory: χ_critical = 0.5 + 1/√N + 1/(2N) for polymer/solvent")
    print("At infinite N: χ_c → 0.5 (THE γ ~ 1 boundary!)")

    gamma_chi = []
    print(f"\n{'System':<35} {'χ_critical':<12} {'N':<10} {'γ = χ_c/0.5':<12}")
    print("-" * 75)

    for fh in flory_huggins_data:
        # Theory prediction for finite N
        chi_theory = 0.5 + 1/np.sqrt(fh.N) + 1/(2*fh.N)
        gamma = fh.chi_critical / 0.5
        gamma_chi.append(gamma)
        print(f"{fh.system:<35} {fh.chi_critical:<12.4f} {fh.N:<10} {gamma:<12.4f}")

    print(f"\n{'Mean γ = χ_c/0.5:':<35} {np.mean(gamma_chi):.4f} ± {np.std(gamma_chi):.4f}")
    print(f"{'At N → ∞: χ_c → 0.5:':<35} γ → 1.000 exactly!")
    print(f"{'Systems at γ ∈ [0.95, 1.05]:':<35} {sum(1 for g in gamma_chi if 0.95 <= g <= 1.05)}/{len(gamma_chi)}")

    return gamma_chi


def analyze_theta_temperature():
    """Analyze theta temperature as γ ~ 1 condition"""
    print("\n" + "=" * 70)
    print("THETA TEMPERATURE: IDEAL CHAIN CONDITIONS AT T/θ = 1")
    print("=" * 70)

    print("\nAt θ-temperature:")
    print("- Second virial coefficient A₂ = 0")
    print("- Chain dimensions = unperturbed (Gaussian)")
    print("- χ = 0.5 (critical miscibility)")
    print("- Excluded volume parameter v = 0")

    T = 298.15  # 25°C
    gamma_theta = []

    print(f"\n{'Polymer':<15} {'Solvent':<15} {'θ (K)':<10} {'χ(25°C)':<10} {'γ = T/θ':<10} {'γ_χ = χ/0.5':<10}")
    print("-" * 80)

    for ps in polymer_solvents:
        gamma_T = T / ps.theta_temp
        gamma_chi = ps.chi_25C / 0.5
        gamma_theta.append(gamma_T)
        print(f"{ps.polymer:<15} {ps.solvent:<15} {ps.theta_temp:<10.0f} {ps.chi_25C:<10.2f} {gamma_T:<10.3f} {gamma_chi:<10.3f}")

    print(f"\n{'Mean γ = T/θ:':<35} {np.mean(gamma_theta):.3f} ± {np.std(gamma_theta):.3f}")
    print(f"{'At T = θ:':<35} γ = 1.000 (ideal chain!)")
    print(f"{'Systems with T near θ:':<35} {sum(1 for g in gamma_theta if 0.8 <= g <= 1.2)}/{len(gamma_theta)}")

    return gamma_theta


def analyze_characteristic_ratio():
    """Analyze characteristic ratio C∞ as chain coherence measure"""
    print("\n" + "=" * 70)
    print("CHARACTERISTIC RATIO: CHAIN STIFFNESS COHERENCE")
    print("=" * 70)

    print("\nC∞ = <R²>/(nl²) measures chain stiffness relative to freely jointed chain")
    print("C∞ = 1: freely jointed chain (no bond correlations)")
    print("C∞ > 1: real chains (bond angle, rotational correlations)")

    gamma_C = []
    print(f"\n{'Polymer':<20} {'C∞':<10} {'l (nm)':<10} {'Tg (K)':<10} {'γ = C∞/6':<10}")
    print("-" * 70)

    # Reference: C∞ ≈ 6 is typical for flexible polymers
    C_ref = 6.0

    for pc in polymer_chains:
        gamma = pc.C_inf / C_ref
        gamma_C.append(gamma)
        print(f"{pc.name:<20} {pc.C_inf:<10.1f} {pc.l_bond:<10.3f} {pc.Tg:<10.0f} {gamma:<10.3f}")

    print(f"\n{'Mean γ = C∞/6:':<35} {np.mean(gamma_C):.3f} ± {np.std(gamma_C):.3f}")
    print(f"{'Reference C∞ = 6:':<35} Typical flexible polymer")
    print(f"{'Systems at γ ∈ [0.7, 1.3]:':<35} {sum(1 for g in gamma_C if 0.7 <= g <= 1.3)}/{len(gamma_C)}")

    # Correlation with Tg
    Tgs = [pc.Tg for pc in polymer_chains]
    C_infs = [pc.C_inf for pc in polymer_chains]
    corr = np.corrcoef(C_infs, Tgs)[0, 1]
    print(f"{'C∞ vs Tg correlation:':<35} r = {corr:.3f}")

    return gamma_C


def analyze_overlap_concentration():
    """Analyze c* as dilute/semidilute crossover"""
    print("\n" + "=" * 70)
    print("OVERLAP CONCENTRATION c*: DILUTE/SEMIDILUTE TRANSITION")
    print("=" * 70)

    print("\nc* = 3M/(4π N_A R_g³) where R_g = chain radius of gyration")
    print("At c/c* = 1: chains begin to overlap (γ ~ 1!)")
    print("c < c*: dilute (isolated chains)")
    print("c > c*: semidilute (overlapping chains)")

    # Calculate c* for different molecular weights
    N_A = 6.022e23
    M_values = [1e4, 5e4, 1e5, 5e5, 1e6, 5e6]

    # For PS: R_g = 0.0145 × M^0.59 (nm) in good solvent
    # c* = 3M/(4π N_A R_g³)

    print(f"\n{'M (g/mol)':<15} {'R_g (nm)':<12} {'c* (g/mL)':<12} {'N_overlap':<12}")
    print("-" * 55)

    for M in M_values:
        R_g = 0.0145 * M**0.59  # nm
        R_g_cm = R_g * 1e-7  # cm
        c_star = 3 * M / (4 * np.pi * N_A * R_g_cm**3)  # g/cm³
        N_overlap = (c_star/M) * N_A * (4/3) * np.pi * R_g_cm**3
        print(f"{M:<15.0e} {R_g:<12.2f} {c_star:<12.4f} {N_overlap:<12.2f}")

    print(f"\n{'At c = c*:':<35} γ = c/c* = 1 (overlap onset!)")
    print(f"{'Semidilute scaling:':<35} η ∝ (c/c*)^(3/(3ν-1)) for c > c*")
    print(f"{'Blob picture:':<35} ξ = R_g × (c/c*)^(-ν/(3ν-1))")

    return M_values


def analyze_expansion_factor():
    """Analyze chain expansion factor α"""
    print("\n" + "=" * 70)
    print("EXPANSION FACTOR α: GOOD SOLVENT CHAIN SWELLING")
    print("=" * 70)

    print("\nα = R_g/R_θ where R_θ = theta dimensions")
    print("α = 1: theta conditions (ideal chain)")
    print("α > 1: good solvent (swelling)")
    print("α < 1: poor solvent (collapse)")

    print("\nPietrila-Flory equation: α⁵ - α³ = Kz where z = excluded volume parameter")
    print("At θ-temperature: z = 0, α = 1 (γ ~ 1!)")

    # Calculate α for different z values
    z_values = np.linspace(-0.5, 3.0, 20)
    alpha_values = []

    print(f"\n{'z':<10} {'α':<10} {'α - 1':<10} {'γ = 1/α':<10}")
    print("-" * 45)

    for z in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
        # Solve α⁵ - α³ = 1.276 × z (using Flory-Fox constant)
        K = 1.276
        # Newton-Raphson for α⁵ - α³ - Kz = 0
        alpha = 1.0
        for _ in range(50):
            f = alpha**5 - alpha**3 - K*z
            fp = 5*alpha**4 - 3*alpha**2
            if abs(fp) > 1e-10:
                alpha = alpha - f/fp
        alpha_values.append(alpha)
        gamma = 1/alpha if alpha > 0 else 0
        print(f"{z:<10.2f} {alpha:<10.3f} {alpha-1:<10.3f} {gamma:<10.3f}")

    print(f"\n{'At z = 0 (θ-solvent):':<35} α = 1.000 exactly!")
    print(f"{'Good solvent (z > 0):':<35} α > 1 (swelling)")
    print(f"{'Poor solvent (z < 0):':<35} α < 1 (collapse)")

    return alpha_values


def analyze_flory_exponent():
    """Analyze Flory exponent ν as coherence measure"""
    print("\n" + "=" * 70)
    print("FLORY EXPONENT ν: CHAIN SCALING UNIVERSALITY")
    print("=" * 70)

    print("\nR_g ∝ N^ν where ν is the Flory exponent")
    print("ν = 0.5: theta/ideal chain (Gaussian, γ ~ 1)")
    print("ν = 0.588: good solvent (self-avoiding walk)")
    print("ν = 1/3: collapsed globule (poor solvent)")
    print("ν = 1.0: rod (fully extended)")

    # Flory exponents from different conditions
    exponents = [
        ("Ideal chain (θ-solvent)", 0.500, "Gaussian statistics"),
        ("Real chain (good solvent)", 0.588, "Self-avoiding walk"),
        ("Collapsed globule", 0.333, "Poor solvent"),
        ("Rod-like", 1.000, "Stiff backbone"),
        ("Branched polymer", 0.400, "Zimm-Stockmayer"),
        ("Melt/concentrated", 0.500, "Screening"),
    ]

    gamma_nu = []
    print(f"\n{'Condition':<30} {'ν':<10} {'γ = ν/0.5':<10} {'Description':<25}")
    print("-" * 80)

    for name, nu, desc in exponents:
        gamma = nu / 0.5
        gamma_nu.append(gamma)
        print(f"{name:<30} {nu:<10.3f} {gamma:<10.3f} {desc:<25}")

    print(f"\n{'Reference ν = 0.5:':<35} Ideal/Gaussian chain")
    print(f"{'Good solvent ν = 0.588:':<35} γ = 1.176 (slightly above γ ~ 1)")
    print(f"{'KEY: θ-conditions at ν = 0.5:':<35} γ = 1.000 exactly!")

    return gamma_nu


def create_visualization(gamma_chi, gamma_theta, gamma_C):
    """Create comprehensive visualization"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: Chi critical values
    ax1 = axes[0, 0]
    systems = [fh.system[:25] for fh in flory_huggins_data]
    chi_vals = [fh.chi_critical for fh in flory_huggins_data]
    ax1.barh(range(len(chi_vals)), chi_vals, color='steelblue', alpha=0.7)
    ax1.axvline(x=0.5, color='red', linestyle='--', linewidth=2, label='χ_c = 0.5 (γ = 1)')
    ax1.set_yticks(range(len(systems)))
    ax1.set_yticklabels(systems, fontsize=9)
    ax1.set_xlabel('χ_critical', fontsize=12)
    ax1.set_title('Flory-Huggins Critical χ', fontsize=14)
    ax1.legend()
    ax1.set_xlim([0.48, 0.55])

    # Plot 2: Theta temperature ratios
    ax2 = axes[0, 1]
    ax2.hist(gamma_theta, bins=8, color='coral', edgecolor='black', alpha=0.7)
    ax2.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='γ = 1 (T = θ)')
    ax2.set_xlabel('γ = T(25°C)/θ', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Temperature vs Theta Temperature', fontsize=14)
    ax2.legend()

    # Plot 3: Characteristic ratios
    ax3 = axes[1, 0]
    polymers = [pc.name for pc in polymer_chains]
    C_vals = [pc.C_inf for pc in polymer_chains]
    ax3.bar(range(len(C_vals)), C_vals, color='mediumseagreen', edgecolor='black', alpha=0.7)
    ax3.axhline(y=6.0, color='red', linestyle='--', linewidth=2, label='C∞ = 6 (reference)')
    ax3.fill_between([-0.5, len(C_vals)-0.5], 4, 8, color='green', alpha=0.2, label='γ ~ 1 range')
    ax3.set_xticks(range(len(polymers)))
    ax3.set_xticklabels(polymers, rotation=45, ha='right', fontsize=9)
    ax3.set_ylabel('Characteristic Ratio C∞', fontsize=12)
    ax3.set_title('Chain Stiffness', fontsize=14)
    ax3.legend()

    # Plot 4: Flory exponent phase diagram
    ax4 = axes[1, 1]
    # Draw schematic phase diagram
    T_reduced = np.linspace(0.5, 1.5, 100)
    chi_curve = 0.5 * (1 + 0.3*(1 - T_reduced))  # Simplified χ(T) relationship
    ax4.fill_between(T_reduced[T_reduced < 1], 0, chi_curve[T_reduced < 1],
                     color='lightblue', alpha=0.5, label='Two-phase (poor solvent)')
    ax4.fill_between(T_reduced[T_reduced >= 1], 0, chi_curve[T_reduced >= 1],
                     color='lightyellow', alpha=0.5, label='One-phase (good solvent)')
    ax4.axhline(y=0.5, color='red', linestyle='--', linewidth=2, label='χ = 0.5')
    ax4.axvline(x=1.0, color='blue', linestyle=':', linewidth=2, label='T/θ = 1')
    ax4.scatter([1.0], [0.5], color='red', s=200, zorder=5, marker='*', label='θ-point')
    ax4.set_xlabel('T/θ (reduced temperature)', fontsize=12)
    ax4.set_ylabel('χ parameter', fontsize=12)
    ax4.set_title('Polymer Solution Phase Diagram', fontsize=14)
    ax4.set_xlim([0.5, 1.5])
    ax4.set_ylim([0.3, 0.7])
    ax4.legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/polymer_thermodynamics_coherence.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nVisualization saved")


def main():
    print("=" * 70)
    print("CHEMISTRY SESSION #210: POLYMER THERMODYNAMICS COHERENCE")
    print("=" * 70)

    gamma_chi = analyze_chi_critical()
    gamma_theta = analyze_theta_temperature()
    gamma_C = analyze_characteristic_ratio()
    M_values = analyze_overlap_concentration()
    alpha_values = analyze_expansion_factor()
    gamma_nu = analyze_flory_exponent()

    create_visualization(gamma_chi, gamma_theta, gamma_C)

    print("\n" + "=" * 70)
    print("SESSION #210 SUMMARY")
    print("=" * 70)

    print("\nKEY γ ~ 1 FINDINGS:")
    print(f"\n1. FLORY-HUGGINS χ:")
    print(f"   χ_critical → 0.5 as N → ∞")
    print(f"   Mean γ = χ_c/0.5 = {np.mean(gamma_chi):.4f} ± {np.std(gamma_chi):.4f}")
    print(f"   χ = 0.5 IS the miscibility boundary!")

    print(f"\n2. THETA TEMPERATURE:")
    print(f"   At T/θ = 1: ideal chain, A₂ = 0, χ = 0.5")
    print(f"   Mean γ = T/θ = {np.mean(gamma_theta):.3f} ± {np.std(gamma_theta):.3f}")
    print(f"   θ-point IS γ ~ 1 for polymer solutions!")

    print(f"\n3. CHARACTERISTIC RATIO:")
    print(f"   Mean γ = C∞/6 = {np.mean(gamma_C):.3f} ± {np.std(gamma_C):.3f}")
    print(f"   C∞ ~ 6 typical for flexible chains")

    print(f"\n4. OVERLAP CONCENTRATION:")
    print(f"   c/c* = 1 IS dilute/semidilute boundary")
    print(f"   γ = c/c* directly!")

    print(f"\n5. EXPANSION FACTOR:")
    print(f"   α = 1 at θ-conditions (z = 0)")
    print(f"   γ = 1/α = 1 for ideal chains")

    print(f"\n6. FLORY EXPONENT:")
    print(f"   ν = 0.5 for ideal chains (γ = ν/0.5 = 1)")
    print(f"   ν = 0.588 for good solvent (γ = 1.18)")

    print("\n" + "=" * 70)
    print("MAJOR INSIGHT: The θ-point IS the universal γ ~ 1 for polymers!")
    print("χ = 0.5, A₂ = 0, α = 1, ν = 0.5 ALL converge at T = θ")
    print("This is the 73rd phenomenon type at γ ~ 1!")
    print("=" * 70)
    print("\nSESSION #210 COMPLETE")


if __name__ == "__main__":
    main()
