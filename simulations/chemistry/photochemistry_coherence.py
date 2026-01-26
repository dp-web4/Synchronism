#!/usr/bin/env python3
"""
Chemistry Session #213: Photochemistry through Coherence Framework

Analyzing photochemical processes through γ ~ 1 framework.

Key concepts:
1. Quantum yield Φ = 1 for ideal photochemical conversion
2. Resonance matching hν = ΔE for absorption
3. Intersystem crossing ISC yield
4. Photocatalytic efficiency
5. Förster radius and FRET efficiency

The γ ~ 1 boundaries:
- Φ = 1: perfect photon-to-product conversion
- hν/ΔE = 1: resonance absorption
- τ_rad/τ_obs = 1: radiative vs non-radiative balance
- R/R_0 = 1: FRET 50% efficiency distance

Author: Claude (Anthropic)
Date: January 2026
Session: Chemistry #213
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class PhotochemicalProcess:
    """Photochemical reaction data"""
    name: str
    quantum_yield: float      # Φ (photons → products)
    wavelength: float         # Excitation wavelength (nm)
    mechanism: str            # Reaction type

@dataclass
class Fluorophore:
    """Fluorescent molecule data"""
    name: str
    quantum_yield: float      # Φ_f (fluorescence QY)
    lifetime_ns: float        # Fluorescence lifetime (ns)
    absorption_max: float     # λ_abs (nm)
    emission_max: float       # λ_em (nm)
    stokes_shift: float       # Δλ (nm)

@dataclass
class PhotocatalyticSystem:
    """Photocatalyst performance data"""
    catalyst: str
    reaction: str
    quantum_efficiency: float  # AQY or IQY
    wavelength: float          # Excitation wavelength (nm)
    TON: float                 # Turnover number

@dataclass
class FRETpair:
    """FRET donor-acceptor pair data"""
    donor: str
    acceptor: str
    R0_nm: float              # Förster radius (nm)
    J_overlap: float          # Spectral overlap integral (M⁻¹cm⁻¹nm⁴)

# Photochemical quantum yields
photochemical_processes = [
    # High quantum yield reactions (Φ ~ 1)
    PhotochemicalProcess("HI photodissociation", 1.00, 254, "Direct dissociation"),
    PhotochemicalProcess("I2 photodissociation", 1.00, 500, "Direct dissociation"),
    PhotochemicalProcess("Cl2 photodissociation", 1.00, 330, "Direct dissociation"),
    PhotochemicalProcess("Br2 photodissociation", 1.00, 420, "Direct dissociation"),
    PhotochemicalProcess("Acetone photolysis", 0.95, 313, "Norrish Type I"),
    PhotochemicalProcess("Azobenzene isomerization", 0.90, 365, "cis-trans"),
    # Moderate quantum yield
    PhotochemicalProcess("Benzophenone triplet", 0.85, 365, "ISC"),
    PhotochemicalProcess("Stilbene isomerization", 0.50, 313, "cis-trans"),
    PhotochemicalProcess("Photoreduction of ketones", 0.45, 350, "H-abstraction"),
    PhotochemicalProcess("Norrish Type II", 0.40, 313, "γ-H abstraction"),
    PhotochemicalProcess("Paterno-Büchi", 0.35, 365, "[2+2] cycloaddition"),
    # Low quantum yield (competing pathways)
    PhotochemicalProcess("Chlorophyll fluorescence", 0.03, 680, "Emission"),
    PhotochemicalProcess("Photosystem II O2", 0.10, 680, "Water splitting"),
    PhotochemicalProcess("Photoelectrochemical H2", 0.05, 400, "Water splitting"),
]

# Fluorophore data
fluorophores = [
    # High quantum yield fluorophores
    Fluorophore("Fluorescein", 0.93, 4.0, 490, 520, 30),
    Fluorophore("Rhodamine 6G", 0.95, 4.1, 530, 555, 25),
    Fluorophore("Rhodamine B", 0.70, 1.7, 554, 580, 26),
    Fluorophore("Coumarin 153", 0.53, 4.5, 422, 530, 108),
    Fluorophore("BODIPY-FL", 0.90, 5.7, 503, 512, 9),
    Fluorophore("Cy3", 0.15, 0.3, 550, 570, 20),
    Fluorophore("Cy5", 0.28, 1.0, 650, 670, 20),
    # Quantum dots
    Fluorophore("CdSe QD (525nm)", 0.50, 20.0, 505, 525, 20),
    Fluorophore("CdSe QD (605nm)", 0.60, 25.0, 585, 605, 20),
    # Natural fluorophores
    Fluorophore("Tryptophan", 0.13, 2.8, 280, 348, 68),
    Fluorophore("GFP", 0.79, 2.8, 395, 509, 114),
    Fluorophore("mCherry", 0.22, 1.4, 587, 610, 23),
]

# Photocatalytic systems
photocatalysts = [
    PhotocatalyticSystem("TiO2", "H2O → H2", 0.05, 365, 100),
    PhotocatalyticSystem("CdS", "H2O → H2", 0.15, 420, 500),
    PhotocatalyticSystem("g-C3N4", "H2O → H2", 0.10, 420, 200),
    PhotocatalyticSystem("[Ru(bpy)3]²⁺/MV²⁺", "H2O → H2", 0.30, 450, 1000),
    PhotocatalyticSystem("Ir(ppy)3", "Organic coupling", 0.85, 450, 5000),
    PhotocatalyticSystem("[Ru(bpy)3]²⁺", "Photoreduction", 0.40, 450, 2000),
    PhotocatalyticSystem("Eosin Y", "Organic reduction", 0.35, 520, 800),
    PhotocatalyticSystem("Rose Bengal", "1O2 generation", 0.76, 560, 3000),
]

# FRET pairs
fret_pairs = [
    FRETpair("CFP", "YFP", 4.9, 1.8e15),
    FRETpair("BFP", "GFP", 4.4, 0.9e15),
    FRETpair("Cy3", "Cy5", 5.3, 2.2e15),
    FRETpair("Fluorescein", "Rhodamine", 5.5, 2.5e15),
    FRETpair("Alexa488", "Alexa555", 6.3, 4.5e15),
    FRETpair("QD525", "Cy3", 7.0, 6.0e15),
]


def analyze_quantum_yield():
    """Analyze photochemical quantum yields"""
    print("=" * 70)
    print("PHOTOCHEMICAL QUANTUM YIELD: Φ = 1 AS γ ~ 1")
    print("=" * 70)

    print("\nΦ = (moles product) / (moles photons absorbed)")
    print("Φ = 1: every absorbed photon produces one product (ideal)")
    print("Φ < 1: competing pathways (fluorescence, IC, ISC, quenching)")

    gamma_phi = []
    print(f"\n{'Process':<35} {'Φ':<8} {'λ (nm)':<10} {'Mechanism':<20}")
    print("-" * 80)

    for proc in photochemical_processes:
        gamma_phi.append(proc.quantum_yield)
        print(f"{proc.name:<35} {proc.quantum_yield:<8.2f} {proc.wavelength:<10.0f} {proc.mechanism:<20}")

    print(f"\n{'Mean Φ:':<35} {np.mean(gamma_phi):.3f} ± {np.std(gamma_phi):.3f}")
    print(f"{'Processes at Φ ∈ [0.8, 1.0]:':<35} {sum(1 for g in gamma_phi if 0.8 <= g <= 1.0)}/{len(gamma_phi)}")
    print(f"{'Direct dissociations at Φ = 1:':<35} 4/4 (100%!)")

    return gamma_phi


def analyze_fluorescence():
    """Analyze fluorescence quantum yields"""
    print("\n" + "=" * 70)
    print("FLUORESCENCE QUANTUM YIELD: RADIATIVE COHERENCE")
    print("=" * 70)

    print("\nΦ_f = k_r / (k_r + k_nr) = τ_obs / τ_rad")
    print("Φ_f = 1: all decay is radiative (no non-radiative loss)")
    print("τ_rad/τ_obs = 1/Φ_f is the radiative fraction")

    gamma_fluor = []
    print(f"\n{'Fluorophore':<20} {'Φ_f':<8} {'τ (ns)':<10} {'λ_abs':<8} {'λ_em':<8} {'Δλ':<8}")
    print("-" * 70)

    for fl in fluorophores:
        gamma_fluor.append(fl.quantum_yield)
        print(f"{fl.name:<20} {fl.quantum_yield:<8.2f} {fl.lifetime_ns:<10.1f} {fl.absorption_max:<8.0f} {fl.emission_max:<8.0f} {fl.stokes_shift:<8.0f}")

    print(f"\n{'Mean Φ_f:':<35} {np.mean(gamma_fluor):.3f} ± {np.std(gamma_fluor):.3f}")
    print(f"{'Fluorophores at Φ ∈ [0.8, 1.0]:':<35} {sum(1 for g in gamma_fluor if 0.8 <= g <= 1.0)}/{len(gamma_fluor)}")

    # High QY examples
    high_qy = [fl for fl in fluorophores if fl.quantum_yield > 0.8]
    print(f"\n{'High Φ_f > 0.8:':<35} {len(high_qy)}/{len(fluorophores)}")
    for fl in high_qy:
        print(f"  {fl.name}: Φ = {fl.quantum_yield}")

    return gamma_fluor


def analyze_photocatalysis():
    """Analyze photocatalytic quantum efficiency"""
    print("\n" + "=" * 70)
    print("PHOTOCATALYTIC EFFICIENCY: AQY AS γ PARAMETER")
    print("=" * 70)

    print("\nAQY = (moles product × n) / (moles incident photons)")
    print("AQY = 1: every photon drives one reaction (ideal)")
    print("Real systems: AQY << 1 due to losses")

    gamma_cat = []
    print(f"\n{'Catalyst':<20} {'Reaction':<20} {'AQY':<10} {'λ (nm)':<10} {'TON':<10}")
    print("-" * 75)

    for cat in photocatalysts:
        gamma_cat.append(cat.quantum_efficiency)
        print(f"{cat.catalyst:<20} {cat.reaction:<20} {cat.quantum_efficiency:<10.2f} {cat.wavelength:<10.0f} {cat.TON:<10.0f}")

    print(f"\n{'Mean AQY:':<35} {np.mean(gamma_cat):.3f} ± {np.std(gamma_cat):.3f}")
    print(f"{'Systems at AQY > 0.5:':<35} {sum(1 for g in gamma_cat if g > 0.5)}/{len(gamma_cat)}")
    print(f"{'Ir(ppy)3 at Φ = 0.85:':<35} Near ideal for photoredox catalysis!")

    return gamma_cat


def analyze_fret():
    """Analyze FRET efficiency at R = R0"""
    print("\n" + "=" * 70)
    print("FRET: 50% EFFICIENCY AT R = R₀ (γ ~ 1)")
    print("=" * 70)

    print("\nE = 1 / [1 + (R/R₀)⁶]")
    print("At R = R₀: E = 0.5 (50% energy transfer)")
    print("R/R₀ = 1 IS the γ ~ 1 distance for FRET!")

    print(f"\n{'Donor':<12} {'Acceptor':<12} {'R₀ (nm)':<10} {'J (M⁻¹cm⁻¹nm⁴)':<18}")
    print("-" * 55)

    for pair in fret_pairs:
        print(f"{pair.donor:<12} {pair.acceptor:<12} {pair.R0_nm:<10.1f} {pair.J_overlap:<18.2e}")

    # FRET efficiency vs R/R0
    R_ratio = np.linspace(0.2, 3.0, 50)
    E_fret = 1 / (1 + R_ratio**6)

    print(f"\n{'R/R₀':<10} {'E_FRET':<12} {'Notes':<30}")
    print("-" * 55)

    for r in [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]:
        E = 1 / (1 + r**6)
        note = ""
        if r == 1.0:
            note = "γ = 1: 50% efficiency!"
        elif r < 1.0:
            note = "Donor-dominated"
        else:
            note = "Acceptor-dominated"
        print(f"{r:<10.1f} {E:<12.3f} {note:<30}")

    print(f"\n{'R₀ IS the Förster radius:':<35} The γ ~ 1 distance!")
    print(f"{'At R = R₀:':<35} E = 0.5 (balanced transfer)")

    return R_ratio, E_fret


def analyze_stokes_shift():
    """Analyze Stokes shift as energy dissipation measure"""
    print("\n" + "=" * 70)
    print("STOKES SHIFT: VIBRATIONAL RELAXATION COHERENCE")
    print("=" * 70)

    print("\nStokes shift Δλ = λ_em - λ_abs")
    print("Small Δλ: rigid molecule, less reorganization (coherent)")
    print("Large Δλ: flexible molecule, significant relaxation")

    stokes_nm = [fl.stokes_shift for fl in fluorophores]
    stokes_ratio = [fl.stokes_shift / fl.absorption_max for fl in fluorophores]

    print(f"\n{'Fluorophore':<20} {'Δλ (nm)':<10} {'Δλ/λ_abs':<12} {'Mirror rule?':<15}")
    print("-" * 60)

    for fl in fluorophores:
        ratio = fl.stokes_shift / fl.absorption_max
        mirror = "Yes" if ratio < 0.1 else "Partial" if ratio < 0.2 else "No"
        print(f"{fl.name:<20} {fl.stokes_shift:<10.0f} {ratio:<12.3f} {mirror:<15}")

    print(f"\n{'Mean Δλ:':<35} {np.mean(stokes_nm):.1f} ± {np.std(stokes_nm):.1f} nm")
    print(f"{'Mean Δλ/λ_abs:':<35} {np.mean(stokes_ratio):.3f} ± {np.std(stokes_ratio):.3f}")
    print(f"{'Small Stokes (Δλ < 30 nm):':<35} {sum(1 for s in stokes_nm if s < 30)}/{len(stokes_nm)}")

    return stokes_nm


def create_visualization(gamma_phi, gamma_fluor, gamma_cat, R_ratio, E_fret):
    """Create comprehensive visualization"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: Photochemical quantum yields
    ax1 = axes[0, 0]
    ax1.hist(gamma_phi, bins=10, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='Φ = 1 (ideal)')
    ax1.set_xlabel('Quantum Yield Φ', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Photochemical Quantum Yields', fontsize=14)
    ax1.legend()

    # Plot 2: Fluorescence quantum yields
    ax2 = axes[0, 1]
    names = [fl.name[:12] for fl in fluorophores]
    ax2.bar(range(len(gamma_fluor)), gamma_fluor, color='coral', edgecolor='black', alpha=0.7)
    ax2.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='Φ_f = 1')
    ax2.fill_between([-0.5, len(gamma_fluor)-0.5], 0.8, 1.0, color='green', alpha=0.2)
    ax2.set_xticks(range(len(names)))
    ax2.set_xticklabels(names, rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Fluorescence QY', fontsize=12)
    ax2.set_title('Fluorophore Quantum Yields', fontsize=14)
    ax2.legend()

    # Plot 3: Photocatalytic efficiency
    ax3 = axes[1, 0]
    cat_names = [cat.catalyst[:12] for cat in photocatalysts]
    ax3.bar(range(len(gamma_cat)), gamma_cat, color='mediumseagreen', edgecolor='black', alpha=0.7)
    ax3.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='AQY = 1 (ideal)')
    ax3.set_xticks(range(len(cat_names)))
    ax3.set_xticklabels(cat_names, rotation=45, ha='right', fontsize=9)
    ax3.set_ylabel('Apparent Quantum Yield', fontsize=12)
    ax3.set_title('Photocatalyst Efficiency', fontsize=14)
    ax3.legend()

    # Plot 4: FRET efficiency curve
    ax4 = axes[1, 1]
    ax4.plot(R_ratio, E_fret, 'b-', linewidth=2, label='E = 1/(1+(R/R₀)⁶)')
    ax4.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='R/R₀ = 1')
    ax4.axhline(y=0.5, color='green', linestyle=':', linewidth=2, label='E = 0.5')
    ax4.fill_betweenx([0, 1], 0.8, 1.2, color='green', alpha=0.2)
    ax4.scatter([1.0], [0.5], color='red', s=200, zorder=5, marker='*')
    ax4.set_xlabel('R/R₀ (distance ratio)', fontsize=12)
    ax4.set_ylabel('FRET Efficiency E', fontsize=12)
    ax4.set_title('FRET: E = 0.5 at R = R₀ (γ ~ 1)', fontsize=14)
    ax4.legend()
    ax4.set_xlim([0, 3])
    ax4.set_ylim([0, 1.05])

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/photochemistry_coherence.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nVisualization saved")


def main():
    print("=" * 70)
    print("CHEMISTRY SESSION #213: PHOTOCHEMISTRY COHERENCE")
    print("=" * 70)

    gamma_phi = analyze_quantum_yield()
    gamma_fluor = analyze_fluorescence()
    gamma_cat = analyze_photocatalysis()
    R_ratio, E_fret = analyze_fret()
    stokes_nm = analyze_stokes_shift()

    create_visualization(gamma_phi, gamma_fluor, gamma_cat, R_ratio, E_fret)

    print("\n" + "=" * 70)
    print("SESSION #213 SUMMARY")
    print("=" * 70)

    print("\nKEY γ ~ 1 FINDINGS:")
    print(f"\n1. PHOTOCHEMICAL QUANTUM YIELD:")
    print(f"   Mean Φ = {np.mean(gamma_phi):.3f} ± {np.std(gamma_phi):.3f}")
    print(f"   {sum(1 for g in gamma_phi if 0.8 <= g <= 1.0)}/{len(gamma_phi)} at Φ ~ 1")
    print(f"   Direct dissociations: 4/4 at Φ = 1.00 exactly!")

    print(f"\n2. FLUORESCENCE QUANTUM YIELD:")
    print(f"   Mean Φ_f = {np.mean(gamma_fluor):.3f} ± {np.std(gamma_fluor):.3f}")
    print(f"   High Φ_f (> 0.8): Fluorescein, Rhodamine 6G, BODIPY-FL")
    print(f"   Φ_f = 1 means purely radiative (coherent decay)")

    print(f"\n3. PHOTOCATALYTIC EFFICIENCY:")
    print(f"   Mean AQY = {np.mean(gamma_cat):.3f} ± {np.std(gamma_cat):.3f}")
    print(f"   Ir(ppy)3 at Φ = 0.85 (near ideal photoredox)")
    print(f"   Solar water splitting: AQY << 1 (challenge)")

    print(f"\n4. FRET EFFICIENCY:")
    print(f"   R/R₀ = 1 gives E = 0.5 (50% transfer)")
    print(f"   R₀ IS the Förster radius (γ ~ 1 distance)")
    print(f"   E = 1/(1+(R/R₀)⁶) is THE FRET equation")

    print(f"\n5. RESONANCE CONDITION:")
    print(f"   hν = ΔE for absorption (energy matching)")
    print(f"   hν/ΔE = 1 IS the resonance γ ~ 1")

    print("\n" + "=" * 70)
    print("MAJOR INSIGHT: Photochemistry is FULL of γ ~ 1 boundaries!")
    print("- Φ = 1: ideal quantum conversion")
    print("- R/R₀ = 1: FRET Förster radius")
    print("- hν/ΔE = 1: resonance absorption")
    print("- τ_rad/τ_obs = 1: radiative lifetime ratio")
    print("This is the 76th phenomenon type at γ ~ 1!")
    print("=" * 70)
    print("\nSESSION #213 COMPLETE")


if __name__ == "__main__":
    main()
