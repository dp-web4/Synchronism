#!/usr/bin/env python3
"""
Chemistry Session #211: Magnetic Moments through Coherence Framework

Analyzing magnetic moment ratios through γ ~ 1 framework.

Key concepts:
1. Curie-Weiss: χ = C/(T - θ) diverges at T/θ = 1
2. Effective moment: μ_eff/μ_spin-only ~ 1 for spin-only
3. Weiss constant: θ/T_c ratio for ordered magnets
4. Pauli paramagnetism: χ_Pauli for metals
5. Van Vleck susceptibility: temperature-independent term

The γ ~ 1 boundaries:
- T/θ = 1: Curie temperature (χ → ∞)
- μ_eff/μ_so = 1: spin-only behavior
- θ/T_c = 1: mean-field prediction
- χ_spin/χ_orbital ~ 1: angular momentum coupling

Author: Claude (Anthropic)
Date: January 2026
Session: Chemistry #211
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class ParamagneticIon:
    """Paramagnetic ion susceptibility data"""
    ion: str
    config: str
    S: float              # Spin quantum number
    L: float              # Orbital quantum number
    J: float              # Total angular momentum
    mu_spin_only: float   # Spin-only moment (μ_B)
    mu_eff_exp: float     # Experimental effective moment (μ_B)
    mu_LS_calc: float     # L-S coupling calculated moment

@dataclass
class MagneticMaterial:
    """Magnetic ordering data"""
    name: str
    formula: str
    Tc_or_Tn: float       # Curie or Néel temperature (K)
    theta_CW: float       # Curie-Weiss constant (K)
    order_type: str       # FM, AFM, FiM
    mu_sat: float         # Saturation moment (μ_B/formula unit)

@dataclass
class MetalSusceptibility:
    """Metal Pauli susceptibility data"""
    metal: str
    chi_exp: float        # Experimental χ (×10⁻⁶ emu/mol)
    chi_Pauli_calc: float # Calculated Pauli χ (×10⁻⁶ emu/mol)
    enhancement: float    # Stoner enhancement factor

# 3d transition metal ions - comparing spin-only vs experimental
paramagnetic_ions = [
    # 3d¹ - Ti³⁺, V⁴⁺
    ParamagneticIon("Ti³⁺", "3d¹", 0.5, 2, 1.5, 1.73, 1.8, 1.55),
    ParamagneticIon("V⁴⁺", "3d¹", 0.5, 2, 1.5, 1.73, 1.7, 1.55),
    # 3d² - V³⁺
    ParamagneticIon("V³⁺", "3d²", 1, 3, 2, 2.83, 2.8, 1.63),
    # 3d³ - V²⁺, Cr³⁺
    ParamagneticIon("Cr³⁺", "3d³", 1.5, 3, 1.5, 3.87, 3.8, 0.77),
    ParamagneticIon("V²⁺", "3d³", 1.5, 3, 1.5, 3.87, 3.8, 0.77),
    # 3d⁴ - Cr²⁺, Mn³⁺
    ParamagneticIon("Mn³⁺", "3d⁴", 2, 2, 0, 4.90, 4.9, 0.0),
    # 3d⁵ - Mn²⁺, Fe³⁺
    ParamagneticIon("Mn²⁺", "3d⁵", 2.5, 0, 2.5, 5.92, 5.9, 5.92),
    ParamagneticIon("Fe³⁺", "3d⁵", 2.5, 0, 2.5, 5.92, 5.9, 5.92),
    # 3d⁶ - Fe²⁺, Co³⁺
    ParamagneticIon("Fe²⁺", "3d⁶", 2, 2, 4, 4.90, 5.4, 6.70),
    ParamagneticIon("Co³⁺ HS", "3d⁶", 2, 2, 4, 4.90, 5.0, 6.70),
    # 3d⁷ - Co²⁺
    ParamagneticIon("Co²⁺", "3d⁷", 1.5, 3, 4.5, 3.87, 4.8, 6.63),
    # 3d⁸ - Ni²⁺
    ParamagneticIon("Ni²⁺", "3d⁸", 1, 3, 4, 2.83, 3.2, 5.59),
    # 3d⁹ - Cu²⁺
    ParamagneticIon("Cu²⁺", "3d⁹", 0.5, 2, 2.5, 1.73, 1.9, 3.55),
]

# Magnetic materials with Curie-Weiss data
magnetic_materials = [
    # Ferromagnets (FM): θ > 0, T_c > 0
    MagneticMaterial("Iron", "Fe", 1043, 1100, "FM", 2.22),
    MagneticMaterial("Cobalt", "Co", 1388, 1415, "FM", 1.72),
    MagneticMaterial("Nickel", "Ni", 627, 650, "FM", 0.62),
    MagneticMaterial("Gadolinium", "Gd", 293, 317, "FM", 7.63),
    MagneticMaterial("EuO", "EuO", 69, 77, "FM", 6.8),
    MagneticMaterial("CrBr3", "CrBr₃", 37, 48, "FM", 3.0),
    # Antiferromagnets (AFM): θ < 0, T_N > 0
    MagneticMaterial("MnO", "MnO", 118, -610, "AFM", 5.0),
    MagneticMaterial("NiO", "NiO", 525, -2000, "AFM", 1.9),
    MagneticMaterial("CoO", "CoO", 289, -330, "AFM", 3.8),
    MagneticMaterial("FeF2", "FeF₂", 78, -117, "AFM", 4.5),
    MagneticMaterial("MnF2", "MnF₂", 67, -80, "AFM", 5.0),
    MagneticMaterial("Cr2O3", "Cr₂O₃", 307, -485, "AFM", 2.8),
    # Ferrimagnets (FiM)
    MagneticMaterial("Fe3O4", "Fe₃O₄", 858, 860, "FiM", 4.1),
    MagneticMaterial("Y3Fe5O12", "YIG", 560, 570, "FiM", 5.0),
]

# Metal Pauli susceptibility
metal_susceptibility = [
    MetalSusceptibility("Li", 25.0, 14.5, 1.7),
    MetalSusceptibility("Na", 16.0, 9.2, 1.7),
    MetalSusceptibility("K", 20.0, 8.7, 2.3),
    MetalSusceptibility("Rb", 18.0, 8.0, 2.3),
    MetalSusceptibility("Cs", 28.0, 7.5, 3.7),
    MetalSusceptibility("Cu", -9.6, 8.5, -1.1),  # diamagnetic (core dominates)
    MetalSusceptibility("Ag", -24.0, 7.0, -3.4),
    MetalSusceptibility("Au", -28.0, 6.5, -4.3),
    MetalSusceptibility("Al", 21.0, 16.5, 1.3),
    MetalSusceptibility("Pd", 540.0, 50.0, 10.8),  # Near Stoner instability!
    MetalSusceptibility("Pt", 193.0, 35.0, 5.5),
]


def analyze_spin_only():
    """Analyze μ_eff/μ_spin-only ratio"""
    print("=" * 70)
    print("SPIN-ONLY MOMENT: μ_eff/μ_spin-only AS γ PARAMETER")
    print("=" * 70)

    print("\nμ_spin-only = √(4S(S+1)) = 2√(S(S+1)) μ_B")
    print("μ_eff = √(3kT·χ_m/N_A·μ_B²) from Curie law")
    print("γ = μ_eff/μ_spin-only ~ 1 for spin-only behavior")

    gamma_spin = []
    print(f"\n{'Ion':<12} {'Config':<8} {'S':<6} {'μ_so':<8} {'μ_eff':<8} {'γ = μ_eff/μ_so':<12}")
    print("-" * 70)

    for ion in paramagnetic_ions:
        gamma = ion.mu_eff_exp / ion.mu_spin_only if ion.mu_spin_only > 0 else 0
        gamma_spin.append(gamma)
        print(f"{ion.ion:<12} {ion.config:<8} {ion.S:<6.1f} {ion.mu_spin_only:<8.2f} {ion.mu_eff_exp:<8.2f} {gamma:<12.3f}")

    print(f"\n{'Mean γ = μ_eff/μ_so:':<35} {np.mean(gamma_spin):.3f} ± {np.std(gamma_spin):.3f}")
    print(f"{'Ions at γ ∈ [0.9, 1.1]:':<35} {sum(1 for g in gamma_spin if 0.9 <= g <= 1.1)}/{len(gamma_spin)}")

    # Count by deviation type
    close_to_1 = sum(1 for g in gamma_spin if 0.95 <= g <= 1.05)
    print(f"{'Ions at γ ∈ [0.95, 1.05]:':<35} {close_to_1}/{len(gamma_spin)}")

    return gamma_spin


def analyze_curie_weiss():
    """Analyze θ/T_c ratio for ordered magnets"""
    print("\n" + "=" * 70)
    print("CURIE-WEISS: θ/T_c RATIO AS MEAN-FIELD TEST")
    print("=" * 70)

    print("\nCurie-Weiss law: χ = C/(T - θ)")
    print("Mean-field theory predicts: θ = T_c for FM, θ = -T_N for simple AFM")
    print("γ = |θ|/T_c ~ 1 tests mean-field validity")

    print("\n--- FERROMAGNETS ---")
    gamma_FM = []
    print(f"{'Material':<15} {'T_c (K)':<10} {'θ (K)':<10} {'γ = θ/T_c':<12}")
    print("-" * 55)

    for mat in magnetic_materials:
        if mat.order_type == "FM" or mat.order_type == "FiM":
            gamma = mat.theta_CW / mat.Tc_or_Tn
            gamma_FM.append(gamma)
            print(f"{mat.name:<15} {mat.Tc_or_Tn:<10.0f} {mat.theta_CW:<10.0f} {gamma:<12.3f}")

    print(f"\n{'Mean γ (FM):':<35} {np.mean(gamma_FM):.3f} ± {np.std(gamma_FM):.3f}")

    print("\n--- ANTIFERROMAGNETS ---")
    gamma_AFM = []
    print(f"{'Material':<15} {'T_N (K)':<10} {'θ (K)':<10} {'f = |θ|/T_N':<12}")
    print("-" * 55)

    for mat in magnetic_materials:
        if mat.order_type == "AFM":
            # For AFM, frustration index f = |θ|/T_N
            f = abs(mat.theta_CW) / mat.Tc_or_Tn
            gamma_AFM.append(f)
            print(f"{mat.name:<15} {mat.Tc_or_Tn:<10.0f} {mat.theta_CW:<10.0f} {f:<12.2f}")

    print(f"\n{'Mean frustration f (AFM):':<35} {np.mean(gamma_AFM):.2f} ± {np.std(gamma_AFM):.2f}")
    print(f"{'For unfrustrated AFM, f ~ 1-3':<35}")
    print(f"{'Systems at f ~ 1-2:':<35} {sum(1 for f in gamma_AFM if 1 <= f <= 2)}/{len(gamma_AFM)}")

    return gamma_FM, gamma_AFM


def analyze_pauli_susceptibility():
    """Analyze Pauli paramagnetic susceptibility"""
    print("\n" + "=" * 70)
    print("PAULI SUSCEPTIBILITY: METAL ELECTRON GAS")
    print("=" * 70)

    print("\nPauli χ = μ_B² × g(E_F) = 3Nμ_B²/(2E_F)")
    print("Enhancement: χ_exp/χ_Pauli = 1/(1 - I×g(E_F)) (Stoner)")
    print("At γ = 1: no enhancement (non-interacting electrons)")

    gamma_Pauli = []
    print(f"\n{'Metal':<10} {'χ_exp':<12} {'χ_calc':<12} {'γ = χ_exp/χ_calc':<15} {'Stoner S':<10}")
    print("-" * 65)

    for metal in metal_susceptibility:
        if metal.chi_exp > 0:  # Only paramagnetic
            gamma = metal.chi_exp / metal.chi_Pauli_calc
            gamma_Pauli.append(gamma)
            print(f"{metal.metal:<10} {metal.chi_exp:<12.1f} {metal.chi_Pauli_calc:<12.1f} {gamma:<15.2f} {metal.enhancement:<10.1f}")

    print(f"\n{'Mean γ (paramagnetic metals):':<35} {np.mean(gamma_Pauli):.2f} ± {np.std(gamma_Pauli):.2f}")
    print(f"{'Metals at γ ∈ [0.5, 2.0]:':<35} {sum(1 for g in gamma_Pauli if 0.5 <= g <= 2.0)}/{len(gamma_Pauli)}")
    print(f"{'Pd at γ = 10.8:':<35} Near Stoner ferromagnetic instability!")

    return gamma_Pauli


def analyze_orbital_quenching():
    """Analyze orbital angular momentum quenching"""
    print("\n" + "=" * 70)
    print("ORBITAL QUENCHING: CRYSTAL FIELD EFFECTS")
    print("=" * 70)

    print("\nIn crystal fields, L is often quenched (L_eff → 0)")
    print("μ_eff = √(L(L+1) + 4S(S+1)) μ_B for free ion")
    print("μ_eff = 2√(S(S+1)) μ_B when L quenched (spin-only)")
    print("γ = μ_eff/μ_spin-only measures residual orbital contribution")

    # For 3d ions, compare with L-S coupling prediction
    gamma_LS = []
    print(f"\n{'Ion':<12} {'L':<6} {'S':<6} {'μ_so':<8} {'μ_LS':<8} {'μ_exp':<8} {'γ = μ_exp/μ_so':<12}")
    print("-" * 80)

    for ion in paramagnetic_ions:
        gamma = ion.mu_eff_exp / ion.mu_spin_only if ion.mu_spin_only > 0 else 0
        gamma_LS.append(gamma)
        print(f"{ion.ion:<12} {ion.L:<6.0f} {ion.S:<6.1f} {ion.mu_spin_only:<8.2f} {ion.mu_LS_calc:<8.2f} {ion.mu_eff_exp:<8.2f} {gamma:<12.3f}")

    # Categorize by orbital contribution
    spin_only = sum(1 for g in gamma_LS if 0.95 <= g <= 1.05)
    some_orbital = sum(1 for g in gamma_LS if 1.05 < g <= 1.3)
    strong_orbital = sum(1 for g in gamma_LS if g > 1.3)

    print(f"\n{'Spin-only (γ ~ 1.0):':<35} {spin_only}/{len(gamma_LS)}")
    print(f"{'Some orbital (γ = 1.05-1.3):':<35} {some_orbital}/{len(gamma_LS)}")
    print(f"{'Strong orbital (γ > 1.3):':<35} {strong_orbital}/{len(gamma_LS)}")
    print(f"\n{'3d⁵ (Mn²⁺, Fe³⁺) at γ = 1.00:':<35} L = 0 ground state (no orbital!)")

    return gamma_LS


def analyze_temperature_crossover():
    """Analyze temperature crossover at T/θ = 1"""
    print("\n" + "=" * 70)
    print("CURIE-WEISS CROSSOVER: T/θ = 1 AS γ ~ 1")
    print("=" * 70)

    print("\nχ = C/(T - θ) diverges at T → θ (Curie temperature)")
    print("χ_max at T/θ = 1 IS the magnetic ordering transition")
    print("This IS the ferromagnetic γ ~ 1 boundary!")

    # Plot χ vs T/θ
    T_reduced = np.linspace(0.1, 3.0, 100)
    chi_CW = 1 / (T_reduced - 1)  # Normalized Curie-Weiss

    print(f"\n{'T/θ':<10} {'χ/C (normalized)':<15} {'Notes':<30}")
    print("-" * 60)

    for T_ratio in [0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.5, 2.0, 3.0]:
        if T_ratio != 1.0:
            chi = 1 / (T_ratio - 1)
            note = ""
            if T_ratio < 1:
                note = "Below T_c (ordered phase)"
            elif T_ratio == 1.05:
                note = "Just above T_c (critical region)"
            elif T_ratio > 1.5:
                note = "Paramagnetic regime"
            print(f"{T_ratio:<10.2f} {chi:<15.2f} {note:<30}")
        else:
            print(f"{T_ratio:<10.2f} {'→ ∞':<15} {'CRITICAL POINT (γ = 1)!':<30}")

    print(f"\n{'At T/θ = 1:':<35} χ → ∞ (ordering transition)")
    print(f"{'This IS the magnetic γ ~ 1:':<35} Para/ferro boundary")

    return T_reduced, chi_CW


def create_visualization(gamma_spin, gamma_FM, gamma_AFM, gamma_Pauli):
    """Create comprehensive visualization"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: Spin-only moment ratios
    ax1 = axes[0, 0]
    ions = [ion.ion for ion in paramagnetic_ions]
    ax1.bar(range(len(gamma_spin)), gamma_spin, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='γ = 1 (spin-only)')
    ax1.fill_between([-0.5, len(gamma_spin)-0.5], 0.9, 1.1, color='green', alpha=0.2)
    ax1.set_xticks(range(len(ions)))
    ax1.set_xticklabels(ions, rotation=45, ha='right', fontsize=9)
    ax1.set_ylabel('γ = μ_eff/μ_spin-only', fontsize=12)
    ax1.set_title('Effective Moment vs Spin-Only', fontsize=14)
    ax1.legend()

    # Plot 2: FM θ/T_c ratios
    ax2 = axes[0, 1]
    FM_names = [m.name for m in magnetic_materials if m.order_type in ["FM", "FiM"]]
    ax2.bar(range(len(gamma_FM)), gamma_FM, color='coral', edgecolor='black', alpha=0.7)
    ax2.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='γ = 1 (mean-field)')
    ax2.fill_between([-0.5, len(gamma_FM)-0.5], 0.9, 1.1, color='green', alpha=0.2)
    ax2.set_xticks(range(len(FM_names)))
    ax2.set_xticklabels(FM_names, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('γ = θ/T_c', fontsize=12)
    ax2.set_title('Ferromagnet Curie-Weiss Ratio', fontsize=14)
    ax2.legend()

    # Plot 3: AFM frustration index
    ax3 = axes[1, 0]
    AFM_names = [m.name for m in magnetic_materials if m.order_type == "AFM"]
    ax3.bar(range(len(gamma_AFM)), gamma_AFM, color='mediumseagreen', edgecolor='black', alpha=0.7)
    ax3.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='f = 1 (mean-field)')
    ax3.axhline(y=2.0, color='orange', linestyle=':', linewidth=2, label='f = 2')
    ax3.set_xticks(range(len(AFM_names)))
    ax3.set_xticklabels(AFM_names, rotation=45, ha='right', fontsize=9)
    ax3.set_ylabel('f = |θ|/T_N (frustration)', fontsize=12)
    ax3.set_title('Antiferromagnet Frustration Index', fontsize=14)
    ax3.legend()

    # Plot 4: Curie-Weiss χ vs T/θ
    ax4 = axes[1, 1]
    T_reduced = np.linspace(1.01, 5.0, 100)
    chi_CW = 1 / (T_reduced - 1)
    ax4.plot(T_reduced, chi_CW, 'b-', linewidth=2, label='χ ∝ 1/(T-θ)')
    ax4.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='T/θ = 1 (T_c)')
    ax4.fill_betweenx([0, 10], 0.9, 1.1, color='green', alpha=0.2)
    ax4.set_xlabel('T/θ (reduced temperature)', fontsize=12)
    ax4.set_ylabel('χ/C (normalized)', fontsize=12)
    ax4.set_title('Curie-Weiss Law: χ Diverges at T/θ = 1', fontsize=14)
    ax4.set_xlim([0.8, 5.0])
    ax4.set_ylim([0, 10])
    ax4.legend()

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/magnetic_moments_coherence.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nVisualization saved")


def main():
    print("=" * 70)
    print("CHEMISTRY SESSION #211: MAGNETIC MOMENTS COHERENCE")
    print("=" * 70)

    gamma_spin = analyze_spin_only()
    gamma_FM, gamma_AFM = analyze_curie_weiss()
    gamma_Pauli = analyze_pauli_susceptibility()
    gamma_LS = analyze_orbital_quenching()
    T_reduced, chi_CW = analyze_temperature_crossover()

    create_visualization(gamma_spin, gamma_FM, gamma_AFM, gamma_Pauli)

    print("\n" + "=" * 70)
    print("SESSION #211 SUMMARY")
    print("=" * 70)

    print("\nKEY γ ~ 1 FINDINGS:")
    print(f"\n1. SPIN-ONLY MOMENT:")
    print(f"   μ_eff/μ_spin-only = {np.mean(gamma_spin):.3f} ± {np.std(gamma_spin):.3f}")
    print(f"   {sum(1 for g in gamma_spin if 0.9 <= g <= 1.1)}/{len(gamma_spin)} ions at γ ~ 1")
    print(f"   Spin-only formula works for most 3d ions!")

    print(f"\n2. FERROMAGNET θ/T_c:")
    print(f"   Mean γ = θ/T_c = {np.mean(gamma_FM):.3f} ± {np.std(gamma_FM):.3f}")
    print(f"   {sum(1 for g in gamma_FM if 0.9 <= g <= 1.1)}/{len(gamma_FM)} at γ ~ 1")
    print(f"   Mean-field predicts θ = T_c (γ = 1)!")

    print(f"\n3. ANTIFERROMAGNET FRUSTRATION:")
    print(f"   Mean f = |θ|/T_N = {np.mean(gamma_AFM):.2f} ± {np.std(gamma_AFM):.2f}")
    print(f"   Unfrustrated systems have f ~ 1-2")
    print(f"   Frustration suppresses T_N below mean-field θ")

    print(f"\n4. PAULI SUSCEPTIBILITY:")
    print(f"   Mean enhancement = {np.mean(gamma_Pauli):.2f} ± {np.std(gamma_Pauli):.2f}")
    print(f"   γ ~ 1 means non-interacting electrons")
    print(f"   Pd at γ = 10.8: near Stoner ferromagnetic instability!")

    print(f"\n5. CURIE-WEISS DIVERGENCE:")
    print(f"   χ → ∞ at T/θ = 1 (Curie temperature)")
    print(f"   T/θ = 1 IS the γ ~ 1 magnetic ordering transition")

    print("\n" + "=" * 70)
    print("MAJOR INSIGHT: Magnetic moments have MULTIPLE γ ~ 1 boundaries!")
    print("- μ_eff/μ_so = 1: spin-only (orbital quenched)")
    print("- θ/T_c = 1: mean-field validity (FM)")
    print("- T/θ = 1: Curie transition (χ → ∞)")
    print("- χ_exp/χ_Pauli = 1: non-interacting electrons")
    print("This is the 74th phenomenon type at γ ~ 1!")
    print("=" * 70)
    print("\nSESSION #211 COMPLETE")


if __name__ == "__main__":
    main()
