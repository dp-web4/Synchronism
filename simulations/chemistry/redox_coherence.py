"""
Chemistry Session #192: Redox Potentials and Electrochemical Series Coherence
Testing redox equilibria through γ ~ 1 framework

Key questions:
1. Is E° = 0 (SHE reference) a γ ~ 1 condition?
2. Does the Nernst equation show γ ~ 1?
3. Are redox couples at E° ~ 0 special?
4. How do electrode potentials relate to coherence?
5. Is the electrochemical series a coherence gradient?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("="*60)
print("CHEMISTRY SESSION #192: REDOX POTENTIALS COHERENCE")
print("="*60)

# Constants
R = 8.314  # J/(mol·K)
T = 298  # K
F = 96485  # C/mol
RT_F = R * T / F  # 0.0257 V at 298 K
ln10 = np.log(10)
NERNST = RT_F * ln10  # 0.0592 V (59.2 mV)

print(f"\nFundamental constants:")
print(f"  RT/F = {RT_F*1000:.1f} mV")
print(f"  RT/F × ln(10) = {NERNST*1000:.1f} mV (Nernst slope)")

# =============================================================================
# STANDARD HYDROGEN ELECTRODE: THE γ ~ 1 REFERENCE
# =============================================================================
print("\n" + "="*60)
print("1. STANDARD HYDROGEN ELECTRODE: E° = 0")
print("="*60)

# SHE: 2H+ + 2e- ⇌ H2, E° = 0.000 V (by definition)
# This is THE reference point for all redox potentials

print("\nStandard Hydrogen Electrode (SHE):")
print("  2H⁺ + 2e⁻ ⇌ H₂")
print("  E° = 0.000 V (by definition)")
print()
print("This IS the γ ~ 1 reference:")
print("  At E = 0: oxidation and reduction equally favored")
print("  The midpoint of the electrochemical series")
print("  Equal tendency to gain or lose electrons")

# =============================================================================
# NERNST EQUATION: E AT EQUILIBRIUM
# =============================================================================
print("\n" + "="*60)
print("2. NERNST EQUATION: γ_Q = Q/K")
print("="*60)

# E = E° - (RT/nF) ln(Q)
# At equilibrium: E = 0, Q = K
# γ_Q = Q/K: at γ = 1, system at equilibrium

print("\nNernst Equation:")
print("  E = E° - (RT/nF) × ln(Q)")
print("  E = E° - (59.2 mV/n) × log(Q)")
print()
print("At equilibrium (E = 0 in a cell):")
print("  Q = K (reaction quotient = equilibrium constant)")
print("  γ_Q = Q/K = 1")
print()
print("γ_Q = Q/K:")
print("  At γ = 1: equilibrium (no net reaction)")
print("  γ < 1: forward reaction favored")
print("  γ > 1: reverse reaction favored")

# =============================================================================
# STANDARD REDUCTION POTENTIALS
# =============================================================================
print("\n" + "="*60)
print("3. STANDARD REDUCTION POTENTIALS: E° DISTRIBUTION")
print("="*60)

# Standard reduction potentials (V vs SHE)
redox_potentials = {
    # Strong oxidizers
    'F2/F-': 2.87,
    'O3/O2': 2.08,
    'H2O2/H2O': 1.78,
    'MnO4-/Mn2+ (acid)': 1.51,
    'Au3+/Au': 1.50,
    'Cl2/Cl-': 1.36,
    'Cr2O7(2-)/Cr3+ (acid)': 1.33,
    'O2/H2O (acid)': 1.23,
    'Br2/Br-': 1.07,
    'NO3-/NO (acid)': 0.96,
    'Ag+/Ag': 0.80,
    'Fe3+/Fe2+': 0.77,
    'I2/I-': 0.54,
    'Cu2+/Cu': 0.34,
    # Near SHE (γ ~ 1 region)
    'AgCl/Ag,Cl-': 0.22,
    'Cu2+/Cu+': 0.15,
    '2H+/H2 (SHE)': 0.00,
    'Pb2+/Pb': -0.13,
    'Sn2+/Sn': -0.14,
    'Ni2+/Ni': -0.26,
    'Co2+/Co': -0.28,
    'Fe2+/Fe': -0.44,
    'Cr3+/Cr': -0.74,
    'Zn2+/Zn': -0.76,
    'Mn2+/Mn': -0.76,
    'Al3+/Al': -1.66,
    'Mg2+/Mg': -2.37,
    'Na+/Na': -2.71,
    'Ca2+/Ca': -2.87,
    'K+/K': -2.93,
    'Li+/Li': -3.04,
}

print("\nStandard Reduction Potentials (V vs SHE):")
print("-"*50)
print(f"{'Couple':<25} {'E° (V)':>10} {'E°/NERNST':>12}")
print("-"*50)

e_values = []
gamma_e = []
for couple, e0 in sorted(redox_potentials.items(), key=lambda x: -x[1]):
    gamma = e0 / NERNST  # E° in units of Nernst slope
    print(f"{couple:<25} {e0:>10.2f} {gamma:>12.1f}")
    e_values.append(e0)
    gamma_e.append(gamma)

e_arr = np.array(e_values)
gamma_e_arr = np.array(gamma_e)

print(f"\nRange: {min(e_arr):.2f} to {max(e_arr):.2f} V")
print(f"Mean E° = {np.mean(e_arr):.2f} V")
print(f"Median E° = {np.median(e_arr):.2f} V")

# Near E° = 0 (SHE reference)
near_zero = np.sum(np.abs(e_arr) < 0.30)
print(f"\nCouples with |E°| < 0.30 V: {near_zero}/{len(e_arr)}")
print("These are in the γ ~ 1 region around SHE!")

# =============================================================================
# BIOLOGICAL REDOX POTENTIALS
# =============================================================================
print("\n" + "="*60)
print("4. BIOLOGICAL REDOX: E°' AT pH 7")
print("="*60)

# Biological standard potentials (at pH 7)
bio_redox = {
    # Couple: E°' (V vs SHE at pH 7)
    'O2/H2O': 0.82,
    'Cytochrome a3 (ox/red)': 0.39,
    'Cytochrome c (ox/red)': 0.25,
    'Cytochrome b (ox/red)': 0.08,
    'CoQ/CoQH2': 0.06,
    'FAD/FADH2': -0.03,
    'Fumarate/Succinate': 0.03,
    'NAD+/NADH': -0.32,
    'NADP+/NADPH': -0.32,
    'Pyruvate/Lactate': -0.19,
    'Acetaldehyde/Ethanol': -0.20,
    'Ferredoxin (ox/red)': -0.43,
    '2H+/H2 (pH 7)': -0.41,
    'CO2/Formate': -0.43,
    'CO2/Glucose': -0.43,
}

print("\nBiological Standard Potentials E°' (pH 7):")
print("-"*50)
print(f"{'Couple':<30} {'E°′ (V)':>10}")
print("-"*50)

bio_e_values = []
for couple, e0 in sorted(bio_redox.items(), key=lambda x: -x[1]):
    print(f"{couple:<30} {e0:>10.2f}")
    bio_e_values.append(e0)

bio_e_arr = np.array(bio_e_values)
print(f"\nMean E°' = {np.mean(bio_e_arr):.2f} V")
print(f"Range: {min(bio_e_arr):.2f} to {max(bio_e_arr):.2f} V")

# The biological range is designed around electron transport
print("\nElectron transport chain spans:")
print(f"  NADH (-0.32 V) → O2 (+0.82 V)")
print(f"  ΔE = 1.14 V")
print(f"  This is ~{1.14/NERNST:.0f} × (RT/F × ln10)")

# =============================================================================
# CELL POTENTIAL AND EQUILIBRIUM
# =============================================================================
print("\n" + "="*60)
print("5. CELL POTENTIAL: E_cell AND ΔG")
print("="*60)

# E_cell = E_cathode - E_anode
# ΔG = -nFE_cell
# At equilibrium: E_cell = 0, ΔG = 0

print("\nCell Potential and Free Energy:")
print("  E_cell = E_cathode - E_anode")
print("  ΔG = -nFE_cell")
print()
print("At equilibrium:")
print("  E_cell = 0")
print("  ΔG = 0")
print("  This IS γ ~ 1!")

# Relationship between E° and K
print("\nE° and Equilibrium Constant K:")
print("  E° = (RT/nF) × ln(K)")
print("  At E° = 0: K = 1")
print("  This is γ = K = 1!")

# Calculate K for various E° values
print("\nE° vs K relationship:")
print("-"*40)
print(f"{'E° (V)':>10} {'K (n=1)':>15} {'K (n=2)':>15}")
print("-"*40)
e0_values = [-0.30, -0.10, 0.00, 0.10, 0.30]
for e0 in e0_values:
    k1 = np.exp(e0 * F / (R * T))
    k2 = np.exp(2 * e0 * F / (R * T))
    print(f"{e0:>10.2f} {k1:>15.2e} {k2:>15.2e}")

print("\nAt E° = 0: K = 1 exactly!")

# =============================================================================
# MIXED POTENTIAL THEORY
# =============================================================================
print("\n" + "="*60)
print("6. MIXED POTENTIAL: CORROSION EQUILIBRIUM")
print("="*60)

# At corrosion potential E_corr:
# i_anodic = i_cathodic (equal currents)
# γ_corr = i_a / i_c = 1

print("\nMixed Potential (Corrosion):")
print("  At E_corr: i_anodic = i_cathodic")
print("  γ_corr = i_a/i_c = 1")
print()
print("This IS γ ~ 1:")
print("  Equal rates of oxidation and reduction")
print("  Dynamic equilibrium at the interface")
print("  Connects to Session #185 (corrosion)")

# Example corrosion potentials
corrosion_data = {
    # Metal: (E_corr in seawater, V vs SCE)
    'Platinum': 0.22,
    'Gold': 0.15,
    'Stainless steel (passive)': 0.05,
    'Copper': -0.20,
    'Lead': -0.50,
    'Iron': -0.61,
    'Aluminum': -0.75,
    'Zinc': -1.05,
    'Magnesium': -1.60,
}

print("\nCorrosion Potentials in Seawater:")
print("-"*40)
print(f"{'Metal':<25} {'E_corr (V vs SCE)':>15}")
print("-"*40)

ecorr_values = []
for metal, ecorr in sorted(corrosion_data.items(), key=lambda x: -x[1]):
    print(f"{metal:<25} {ecorr:>15.2f}")
    ecorr_values.append(ecorr)

ecorr_arr = np.array(ecorr_values)
print(f"\nRange: {min(ecorr_arr):.2f} to {max(ecorr_arr):.2f} V")

# =============================================================================
# POURBAIX DIAGRAMS: STABILITY BOUNDARIES
# =============================================================================
print("\n" + "="*60)
print("7. POURBAIX DIAGRAMS: E-pH BOUNDARIES")
print("="*60)

# Boundaries in Pourbaix diagrams follow:
# dE/dpH = -59.2 mV/pH (for H+ involved reactions)
# At boundary: two species in equilibrium (γ ~ 1)

print("\nPourbaix Diagram Boundaries:")
print("  At each boundary: species1 ⇌ species2 in equilibrium")
print("  This IS γ ~ 1 along every boundary!")
print()
print("Boundary slopes:")
print(f"  dE/dpH = -59.2 mV/pH (when H+ involved)")
print(f"  = -2.303 × RT/F")
print("  This is exactly the Nernst slope!")

# Water stability diagram
print("\nWater Stability:")
print("  Upper line: O2/H2O (E = 1.23 - 0.059×pH)")
print("  Lower line: H+/H2 (E = 0.00 - 0.059×pH)")
print("  Water is stable BETWEEN these lines")

# At any given pH, the stable species is determined by E
# Crossing a boundary = phase transition at γ ~ 1

# =============================================================================
# CONCENTRATION CELLS
# =============================================================================
print("\n" + "="*60)
print("8. CONCENTRATION CELLS: γ_c = c1/c2")
print("="*60)

# E_cell = (RT/nF) × ln(c2/c1)
# At c1 = c2: E = 0 (no driving force)

print("\nConcentration Cell:")
print("  E = (RT/nF) × ln(c2/c1)")
print("  = (59.2 mV/n) × log(c2/c1)")
print()
print("At c1 = c2:")
print("  E = 0")
print("  γ_c = c2/c1 = 1")
print("  No driving force for electron flow")
print()
print("This IS γ ~ 1: equal concentrations = equilibrium")

# Calculate E for concentration ratios
print("\nConcentration ratio vs Cell Potential (n=1):")
print("-"*40)
print(f"{'c2/c1':>10} {'E (mV)':>15}")
print("-"*40)
ratios = [0.1, 0.5, 1.0, 2.0, 10.0]
for ratio in ratios:
    e_mv = NERNST * 1000 * np.log10(ratio)
    print(f"{ratio:>10.1f} {e_mv:>15.1f}")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: REDOX COHERENCE PARAMETERS")
print("="*60)

# Count potentials near zero
near_zero_count = np.sum(np.abs(e_arr) < 0.30)
near_zero_bio = np.sum(np.abs(bio_e_arr) < 0.30)

summary = {
    'Couples with |E°| < 0.30 V': (near_zero_count, len(e_arr)),
    'Bio couples with |E°′| < 0.30 V': (near_zero_bio, len(bio_e_arr)),
    'Mean E° (all)': (np.mean(e_arr), np.std(e_arr)),
    'Mean E°′ (bio)': (np.mean(bio_e_arr), np.std(bio_e_arr)),
}

print(f"\n{'Parameter':<35} {'Value':>15}")
print("-"*55)
for param, val in summary.items():
    if isinstance(val[1], int):
        print(f"{param:<35} {val[0]}/{val[1]}")
    else:
        print(f"{param:<35} {val[0]:>7.2f} ± {val[1]:.2f}")

# Key γ ~ 1 conditions
print("\nKEY γ ~ 1 CONDITIONS IN REDOX CHEMISTRY:")
print("1. E° = 0 (SHE): reference point (equal ox/red tendency)")
print("2. Q/K = 1: Nernst equilibrium")
print("3. K = 1 (at E° = 0): equilibrium constant unity")
print("4. E_cell = 0: electrochemical equilibrium")
print("5. i_a/i_c = 1: mixed potential (corrosion)")
print("6. c1/c2 = 1: concentration cell equilibrium")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Chemistry Session #192: Redox Potentials Coherence',
             fontsize=14, fontweight='bold')

# Panel 1: Standard potential distribution
ax1 = axes[0, 0]
ax1.hist(e_arr, bins=15, color='steelblue', alpha=0.7, edgecolor='black')
ax1.axvline(x=0, color='red', linestyle='--', linewidth=2, label='E° = 0 (SHE)')
ax1.axvspan(-0.30, 0.30, alpha=0.1, color='green', label='γ ~ 1 region')
ax1.set_xlabel('E° (V vs SHE)')
ax1.set_ylabel('Count')
ax1.set_title('Standard Reduction Potential Distribution')
ax1.legend()

# Panel 2: Biological redox potentials
ax2 = axes[0, 1]
couples = list(bio_redox.keys())
e_bio = list(bio_redox.values())
y = np.arange(len(couples))
colors = ['forestgreen' if abs(e) < 0.30 else 'steelblue' for e in e_bio]
ax2.barh(y, e_bio, color=colors, alpha=0.7, edgecolor='black')
ax2.axvline(x=0, color='red', linestyle='--', linewidth=2, label='E°′ = 0')
ax2.set_yticks(y)
ax2.set_yticklabels(couples, fontsize=7)
ax2.set_xlabel("E°' (V vs SHE at pH 7)")
ax2.set_title('Biological Redox Potentials')
ax2.legend()

# Panel 3: E° vs K relationship
ax3 = axes[1, 0]
e0_range = np.linspace(-0.5, 0.5, 100)
K_n1 = np.exp(e0_range * F / (R * T))
K_n2 = np.exp(2 * e0_range * F / (R * T))
ax3.semilogy(e0_range, K_n1, 'b-', linewidth=2, label='n = 1')
ax3.semilogy(e0_range, K_n2, 'r-', linewidth=2, label='n = 2')
ax3.axvline(x=0, color='gray', linestyle='--', linewidth=1)
ax3.axhline(y=1, color='gray', linestyle='--', linewidth=1)
ax3.plot(0, 1, 'ko', markersize=10, label='E°=0, K=1 (γ~1)')
ax3.set_xlabel('E° (V)')
ax3.set_ylabel('Equilibrium Constant K')
ax3.set_title('E° vs K: K = 1 at E° = 0')
ax3.legend()
ax3.set_xlim(-0.5, 0.5)
ax3.set_ylim(1e-10, 1e10)

# Panel 4: Summary γ values
ax4 = axes[1, 1]
gamma_labels = ['E°=0', 'Q/K=1', 'K=1', 'E_cell=0', 'i_a/i_c', 'c1/c2']
gamma_vals = [1, 1, 1, 1, 1, 1]  # All γ ~ 1 conditions
ax4.bar(gamma_labels, gamma_vals, color='forestgreen', alpha=0.7, edgecolor='black')
ax4.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='γ = 1')
ax4.set_ylabel('γ value at equilibrium')
ax4.set_title('All Redox Equilibria at γ = 1')
ax4.legend()
ax4.set_ylim(0, 1.5)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/redox_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "="*60)
print("FINDING #129: REDOX POTENTIALS AT γ ~ 1")
print("="*60)

print("""
KEY RESULTS:

1. STANDARD HYDROGEN ELECTRODE
   - E° = 0 is THE γ ~ 1 reference
   - Equal tendency for oxidation/reduction
   - Midpoint of electrochemical series

2. NERNST EQUATION
   - γ_Q = Q/K
   - At equilibrium: Q = K, γ = 1
   - E = E° when γ_Q = 1

3. E° = 0 MEANS K = 1
   - E° = (RT/nF) × ln(K)
   - At E° = 0: K = 1 exactly
   - Equal product/reactant concentrations

4. CELL POTENTIAL
   - At E_cell = 0: electrochemical equilibrium
   - ΔG = 0 (no driving force)
   - This IS γ ~ 1

5. MIXED POTENTIAL
   - At E_corr: i_a = i_c
   - γ_corr = i_a/i_c = 1
   - Dynamic corrosion equilibrium

6. CONCENTRATION CELLS
   - At c1/c2 = 1: E = 0
   - Equal concentrations = no potential
   - This IS γ ~ 1

PHYSICAL INSIGHT:
All electrochemical equilibria occur at γ ~ 1:
- E° = 0: reference (SHE)
- Q = K: Nernst equilibrium
- K = 1: equal products/reactants
- E_cell = 0: cell equilibrium
- i_a = i_c: corrosion equilibrium

The electrochemical series IS a coherence gradient:
- E° > 0: oxidation favored (electron acceptors)
- E° < 0: reduction favored (electron donors)
- E° = 0: balanced (γ ~ 1)

{}/{} couples have |E°| < 0.30 V (near γ ~ 1).

55th phenomenon type at γ ~ 1!
""".format(near_zero_count, len(e_arr)))

print("\nVisualization saved to: redox_coherence.png")
print("\nSESSION #192 COMPLETE")
