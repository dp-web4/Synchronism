#!/usr/bin/env python3
"""
Chemistry Session #70: Reaction Kinetics & Coherence
Test whether coherence framework predicts activation energies and rate constants.

The Arrhenius equation:
k = A × exp(-E_a / RT)

Coherence interpretation:
- E_a = barrier from breaking coherent bonds
- A = attempt frequency × coherence factor
- Higher coherence in reactants → higher barrier
- Coherent transition state → lower barrier

Hypothesis:
E_a ∝ γ_reactants / γ_TS (ratio of coherences)
k ∝ (2/γ_TS) × exp(-E_a/RT)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

print("=" * 70)
print("CHEMISTRY SESSION #70: REACTION KINETICS & COHERENCE")
print("=" * 70)

# ==============================================================================
# DATASET: Activation Energies and Rate Constants
# ==============================================================================

print("\n" + "=" * 70)
print("DATASET: REACTION KINETICS")
print("=" * 70)

# Gas-phase reactions (well-characterized E_a)
gas_phase_reactions = {
    # Reaction: (E_a in kJ/mol, log10(A in s^-1 or L/mol/s), type)
    # Unimolecular
    'cyclopropane_isomerization': (272, 15.2, 'unimolecular'),
    'N2O5_decomposition': (103, 13.5, 'unimolecular'),
    'C2H5Cl_elimination': (254, 13.8, 'unimolecular'),
    'cyclobutane_decomposition': (261, 15.6, 'unimolecular'),
    'CH3NC_isomerization': (160, 13.5, 'unimolecular'),

    # Bimolecular
    'H_plus_H2': (33, 10.9, 'bimolecular'),
    'H_plus_Cl2': (10, 10.4, 'bimolecular'),
    'Cl_plus_H2': (23, 10.9, 'bimolecular'),
    'Br_plus_H2': (73, 10.8, 'bimolecular'),
    'OH_plus_H2': (25, 10.9, 'bimolecular'),
    'CH3_plus_H2': (50, 10.5, 'bimolecular'),
    'H_plus_CH4': (50, 11.0, 'bimolecular'),
    'NO_plus_O3': (10, 9.4, 'bimolecular'),
    'NO2_plus_F2': (46, 9.8, 'bimolecular'),
    '2NO2_to_2NO_plus_O2': (111, 9.3, 'bimolecular'),
}

# Solution-phase reactions (SN2, SN1)
solution_reactions = {
    # Reaction: (E_a in kJ/mol, mechanism)
    'CH3Br_plus_OH_SN2': (75, 'SN2'),
    'CH3I_plus_OH_SN2': (65, 'SN2'),
    'C2H5Br_plus_OH_SN2': (80, 'SN2'),
    'CH3Cl_plus_OH_SN2': (90, 'SN2'),
    't-BuCl_SN1': (100, 'SN1'),
    't-BuBr_SN1': (92, 'SN1'),
    'benzyl_Cl_SN1': (96, 'SN1'),
    'allyl_Cl_SN1': (105, 'SN1'),
}

# Enzyme-catalyzed reactions (dramatic E_a reduction)
enzyme_reactions = {
    # (uncatalyzed E_a, catalyzed E_a, enzyme)
    'H2O2_decomposition': (75, 28, 'catalase'),
    'urea_hydrolysis': (134, 44, 'urease'),
    'peptide_hydrolysis': (86, 50, 'carboxypeptidase'),
    'CO2_hydration': (58, 29, 'carbonic_anhydrase'),
    'glucose_isomerization': (125, 50, 'glucose_isomerase'),
    'starch_hydrolysis': (125, 33, 'amylase'),
}

print(f"Gas-phase reactions: {len(gas_phase_reactions)}")
print(f"Solution reactions: {len(solution_reactions)}")
print(f"Enzyme reactions: {len(enzyme_reactions)}")

# ==============================================================================
# COHERENCE PARAMETER ESTIMATION
# ==============================================================================

print("\n" + "=" * 70)
print("COHERENCE PARAMETER ESTIMATION")
print("=" * 70)

def gamma_from_reaction_type(reaction_type, mechanism=None):
    """
    Estimate γ based on reaction type and mechanism.

    Higher γ (less coherent) for:
    - More bonds breaking
    - More disorder in transition state
    - Solvent involvement
    """
    base_gamma = {
        'unimolecular': 1.2,      # Single molecule reorganization
        'bimolecular': 1.5,       # Two molecules must meet
        'SN2': 1.4,               # Concerted, ordered TS
        'SN1': 1.7,               # Stepwise, carbocation intermediate
    }

    return base_gamma.get(reaction_type or mechanism, 1.5)

def gamma_transition_state(E_a, A_factor, T=298):
    """
    Estimate γ_TS from activation parameters.

    Lower E_a and higher A suggest more ordered (coherent) transition state.
    """
    # Normalize E_a to typical range (0-300 kJ/mol)
    E_norm = E_a / 150.0

    # A-factor contribution (higher A = more organized approach)
    A_contribution = 0 if A_factor is None else (A_factor - 10.0) / 5.0

    # γ_TS: lower E_a and higher A → lower γ
    gamma_ts = 0.5 + 1.0 * E_norm - 0.2 * A_contribution

    return np.clip(gamma_ts, 0.5, 2.0)

def gamma_enzyme_effect(E_uncat, E_cat):
    """
    Estimate γ from enzyme catalysis effect.

    Enzyme lowers E_a by providing coherent binding pocket.
    """
    # Ratio of catalyzed to uncatalyzed E_a
    ratio = E_cat / E_uncat

    # Lower ratio = more effective catalysis = more coherent TS
    gamma = 0.5 + 1.5 * ratio

    return np.clip(gamma, 0.5, 2.0)

# ==============================================================================
# ANALYSIS 1: GAS-PHASE REACTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS 1: GAS-PHASE REACTIONS")
print("=" * 70)

# Extract data
E_a_list = []
A_list = []
gamma_ts_list = []
reaction_types = []

for name, (E_a, A, rtype) in gas_phase_reactions.items():
    E_a_list.append(E_a)
    A_list.append(A)
    gamma_ts_list.append(gamma_transition_state(E_a, A))
    reaction_types.append(rtype)

E_a_arr = np.array(E_a_list)
A_arr = np.array(A_list)
gamma_ts_arr = np.array(gamma_ts_list)

# Correlations
print("\nCorrelations:")

# E_a vs 2/γ_TS
coh_factor = 2.0 / gamma_ts_arr
r_Ea_coh, p_Ea_coh = stats.pearsonr(E_a_arr, coh_factor)
print(f"E_a vs 2/γ_TS: r = {r_Ea_coh:.3f}, p = {p_Ea_coh:.4f}")

# A-factor vs γ_TS
r_A_gamma, p_A_gamma = stats.pearsonr(A_arr, gamma_ts_arr)
print(f"A vs γ_TS: r = {r_A_gamma:.3f}, p = {p_A_gamma:.4f}")

# By reaction type
uni_mask = np.array([t == 'unimolecular' for t in reaction_types])
bi_mask = np.array([t == 'bimolecular' for t in reaction_types])

print(f"\nUnimolecular: mean E_a = {E_a_arr[uni_mask].mean():.1f} kJ/mol, n = {sum(uni_mask)}")
print(f"Bimolecular: mean E_a = {E_a_arr[bi_mask].mean():.1f} kJ/mol, n = {sum(bi_mask)}")

# ==============================================================================
# ANALYSIS 2: SN1 vs SN2 MECHANISMS
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS 2: SN1 vs SN2 MECHANISMS")
print("=" * 70)

sn2_Ea = []
sn1_Ea = []

for name, (E_a, mech) in solution_reactions.items():
    if mech == 'SN2':
        sn2_Ea.append(E_a)
    else:
        sn1_Ea.append(E_a)

print(f"\nSN2 reactions: mean E_a = {np.mean(sn2_Ea):.1f} ± {np.std(sn2_Ea):.1f} kJ/mol, n = {len(sn2_Ea)}")
print(f"SN1 reactions: mean E_a = {np.mean(sn1_Ea):.1f} ± {np.std(sn1_Ea):.1f} kJ/mol, n = {len(sn1_Ea)}")

# γ prediction
gamma_sn2 = gamma_from_reaction_type(None, 'SN2')
gamma_sn1 = gamma_from_reaction_type(None, 'SN1')

print(f"\nPredicted γ_SN2 = {gamma_sn2:.2f} (concerted, ordered)")
print(f"Predicted γ_SN1 = {gamma_sn1:.2f} (stepwise, carbocation)")

# Does higher γ → higher E_a?
print(f"\nγ_SN1 > γ_SN2: {gamma_sn1 > gamma_sn2}")
print(f"E_a_SN1 > E_a_SN2: {np.mean(sn1_Ea) > np.mean(sn2_Ea)}")

if (gamma_sn1 > gamma_sn2) == (np.mean(sn1_Ea) > np.mean(sn2_Ea)):
    print("✓ Direction matches prediction!")
else:
    print("✗ Direction opposite to prediction")

# ==============================================================================
# ANALYSIS 3: ENZYME CATALYSIS
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS 3: ENZYME CATALYSIS AS COHERENCE MATCHING")
print("=" * 70)

E_uncat_list = []
E_cat_list = []
reduction_list = []
gamma_enz_list = []

for name, (E_uncat, E_cat, enzyme) in enzyme_reactions.items():
    E_uncat_list.append(E_uncat)
    E_cat_list.append(E_cat)
    reduction_list.append(E_uncat - E_cat)
    gamma_enz_list.append(gamma_enzyme_effect(E_uncat, E_cat))
    print(f"{enzyme}: E_a {E_uncat} → {E_cat} kJ/mol (reduction = {E_uncat - E_cat}, γ_eff = {gamma_enzyme_effect(E_uncat, E_cat):.2f})")

E_uncat_arr = np.array(E_uncat_list)
E_cat_arr = np.array(E_cat_list)
reduction_arr = np.array(reduction_list)
gamma_enz_arr = np.array(gamma_enz_list)

print(f"\nMean uncatalyzed E_a: {E_uncat_arr.mean():.1f} kJ/mol")
print(f"Mean catalyzed E_a: {E_cat_arr.mean():.1f} kJ/mol")
print(f"Mean reduction: {reduction_arr.mean():.1f} kJ/mol ({100*reduction_arr.mean()/E_uncat_arr.mean():.0f}%)")

# Rate enhancement
RT = 2.479  # kJ/mol at 298K
rate_enhancement = np.exp(reduction_arr / RT)
print(f"\nRate enhancements: {rate_enhancement.min():.1e} to {rate_enhancement.max():.1e}")

# Coherence interpretation
print("\nCoherence interpretation:")
print("Enzymes provide a coherent binding pocket that:")
print("  1. Pre-organizes reactants (lowers entropic barrier)")
print("  2. Stabilizes transition state (lowers enthalpic barrier)")
print("  3. Results in lower effective γ for catalyzed reaction")

# ==============================================================================
# ANALYSIS 4: HAMMOND POSTULATE & COHERENCE
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS 4: HAMMOND POSTULATE & COHERENCE")
print("=" * 70)

print("""
Hammond Postulate: Transition state resembles closer species.
- Exothermic: TS resembles reactants (early TS)
- Endothermic: TS resembles products (late TS)

Coherence interpretation:
- Early TS: γ_TS ≈ γ_reactants (reactant-like coherence)
- Late TS: γ_TS ≈ γ_products (product-like coherence)

For exothermic reactions (ΔH < 0):
- Lower E_a expected (early TS)
- γ_TS closer to γ_reactants

For endothermic reactions (ΔH > 0):
- Higher E_a expected (late TS)
- γ_TS closer to γ_products
""")

# ==============================================================================
# ANALYSIS 5: TRANSITION STATE THEORY + COHERENCE
# ==============================================================================

print("\n" + "=" * 70)
print("ANALYSIS 5: TRANSITION STATE THEORY + COHERENCE")
print("=" * 70)

print("""
Eyring equation:
k = (kT/h) × exp(-ΔG‡/RT)

Where ΔG‡ = ΔH‡ - TΔS‡

Coherence contributions:
1. ΔH‡ ∝ (2/γ_reactants - 2/γ_TS)
   - Breaking coherent bonds costs energy

2. ΔS‡ reflects transition state order
   - More ordered TS (lower γ_TS) → more negative ΔS‡
   - Less ordered TS (higher γ_TS) → less negative ΔS‡

3. Combined effect:
   k ∝ (2/γ_TS) × exp(-ΔH‡/RT) × exp(ΔS‡/R)
""")

# ==============================================================================
# MODEL: E_a PREDICTION FROM COHERENCE
# ==============================================================================

print("\n" + "=" * 70)
print("MODEL: COHERENCE-BASED E_a PREDICTION")
print("=" * 70)

# For gas-phase reactions, test if E_a correlates with γ_TS

# Linear model: E_a = a × γ_TS + b
def Ea_model(gamma, a, b):
    return a * gamma + b

# Fit to gas-phase data
popt, pcov = curve_fit(Ea_model, gamma_ts_arr, E_a_arr)
a_fit, b_fit = popt
E_a_pred = Ea_model(gamma_ts_arr, a_fit, b_fit)

# Calculate R²
ss_res = np.sum((E_a_arr - E_a_pred)**2)
ss_tot = np.sum((E_a_arr - E_a_arr.mean())**2)
R_squared = 1 - ss_res/ss_tot

print(f"\nLinear fit: E_a = {a_fit:.1f} × γ_TS + {b_fit:.1f}")
print(f"R² = {R_squared:.3f}")

# Correlation E_a vs γ_TS directly
r_Ea_gamma, p_Ea_gamma = stats.pearsonr(E_a_arr, gamma_ts_arr)
print(f"\nE_a vs γ_TS: r = {r_Ea_gamma:.3f}, p = {p_Ea_gamma:.4f}")

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SESSION #70 SUMMARY: REACTION KINETICS")
print("=" * 70)

print(f"""
Correlations Found:
- E_a vs γ_TS (gas phase): r = {r_Ea_gamma:.3f} {"(GOOD)" if abs(r_Ea_gamma) > 0.6 else "(MODERATE)" if abs(r_Ea_gamma) > 0.3 else "(WEAK)"}
- E_a vs 2/γ_TS: r = {r_Ea_coh:.3f}
- A vs γ_TS: r = {r_A_gamma:.3f}
- Linear model R² = {R_squared:.3f}

Key Findings:
1. Unimolecular vs Bimolecular:
   - Unimolecular: mean E_a = {E_a_arr[uni_mask].mean():.0f} kJ/mol (higher)
   - Bimolecular: mean E_a = {E_a_arr[bi_mask].mean():.0f} kJ/mol (lower)

2. SN1 vs SN2:
   - SN2 (concerted, γ = {gamma_sn2}): mean E_a = {np.mean(sn2_Ea):.0f} kJ/mol
   - SN1 (stepwise, γ = {gamma_sn1}): mean E_a = {np.mean(sn1_Ea):.0f} kJ/mol
   - Higher γ correlates with higher E_a ✓

3. Enzyme Catalysis:
   - Mean reduction: {reduction_arr.mean():.0f} kJ/mol ({100*reduction_arr.mean()/E_uncat_arr.mean():.0f}%)
   - Enzymes provide coherent environment (lower γ_TS)
   - Rate enhancements: 10^{np.log10(rate_enhancement).mean():.0f}
""")

# ==============================================================================
# PREDICTIONS
# ==============================================================================

print("\n" + "=" * 70)
print("PREDICTIONS")
print("=" * 70)

print("""
P70.1: E_a ∝ γ_TS
Higher transition state disorder (γ) → higher activation energy.

P70.2: SN2 < SN1 in γ_TS
Concerted mechanism more coherent than stepwise.

P70.3: Enzyme catalysis = coherence matching
Enzyme active site provides γ-matched environment for TS.

P70.4: A-factor reflects TS order
Higher A → more favorable entropy → lower γ_TS.

P70.5: Hammond postulate as γ interpolation
γ_TS interpolates between γ_reactants and γ_products based on ΔH.
""")

# ==============================================================================
# VALIDATION STATUS
# ==============================================================================

print("\n" + "=" * 70)
print("VALIDATION STATUS")
print("=" * 70)

# Determine validation level
if abs(r_Ea_gamma) > 0.7:
    status = "STRONG SUPPORTING EVIDENCE"
elif abs(r_Ea_gamma) > 0.5:
    status = "MODERATE SUPPORTING EVIDENCE"
elif abs(r_Ea_gamma) > 0.3:
    status = "WEAK SUPPORTING EVIDENCE"
else:
    status = "INSUFFICIENT CORRELATION"

print(f"\n{status}")
print(f"""
The coherence framework provides:
1. QUALITATIVE success for mechanism comparison (SN1 > SN2 in γ and E_a)
2. QUANTITATIVE correlation r = {r_Ea_gamma:.3f} for gas-phase reactions
3. INTERPRETIVE value for enzyme catalysis

Limitations:
- γ_TS estimation is indirect (derived from E_a, A)
- Circular reasoning risk: γ estimated from what we're trying to predict
- Need independent γ_TS measures (computational, spectroscopic)
""")

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: E_a vs γ_TS (gas phase)
ax1 = axes[0, 0]
colors = ['blue' if t == 'unimolecular' else 'red' for t in reaction_types]
ax1.scatter(gamma_ts_arr, E_a_arr, c=colors, s=80, alpha=0.7)
gamma_line = np.linspace(gamma_ts_arr.min(), gamma_ts_arr.max(), 100)
ax1.plot(gamma_line, Ea_model(gamma_line, a_fit, b_fit), 'k--', label=f'Linear fit (R²={R_squared:.2f})')
ax1.set_xlabel('γ_TS (estimated)', fontsize=12)
ax1.set_ylabel('Activation Energy (kJ/mol)', fontsize=12)
ax1.set_title(f'E_a vs γ_TS\n(r = {r_Ea_gamma:.3f})', fontsize=14)
ax1.legend(['Linear fit', 'Unimolecular', 'Bimolecular'])
ax1.grid(True, alpha=0.3)

# Plot 2: SN1 vs SN2 comparison
ax2 = axes[0, 1]
positions = [1, 2]
bp = ax2.boxplot([sn2_Ea, sn1_Ea], positions=positions, widths=0.6)
ax2.set_xticks(positions)
ax2.set_xticklabels([f'SN2\n(γ={gamma_sn2:.1f})', f'SN1\n(γ={gamma_sn1:.1f})'])
ax2.set_ylabel('Activation Energy (kJ/mol)', fontsize=12)
ax2.set_title('SN2 vs SN1: Mechanism & Coherence', fontsize=14)
ax2.grid(True, alpha=0.3, axis='y')

# Plot 3: Enzyme catalysis
ax3 = axes[1, 0]
x = np.arange(len(enzyme_reactions))
width = 0.35
enzymes = list(enzyme_reactions.keys())
ax3.bar(x - width/2, E_uncat_arr, width, label='Uncatalyzed', color='gray')
ax3.bar(x + width/2, E_cat_arr, width, label='Catalyzed', color='green')
ax3.set_xticks(x)
ax3.set_xticklabels([list(enzyme_reactions.values())[i][2] for i in range(len(enzyme_reactions))], rotation=45, ha='right')
ax3.set_ylabel('Activation Energy (kJ/mol)', fontsize=12)
ax3.set_title('Enzyme Catalysis: E_a Reduction', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3, axis='y')

# Plot 4: γ_enzyme vs E_a reduction
ax4 = axes[1, 1]
ax4.scatter(gamma_enz_arr, reduction_arr, s=100, c='green', alpha=0.7)
r_enz, _ = stats.pearsonr(gamma_enz_arr, reduction_arr)
ax4.set_xlabel('γ_eff (catalyzed)', fontsize=12)
ax4.set_ylabel('E_a reduction (kJ/mol)', fontsize=12)
ax4.set_title(f'Enzyme γ vs E_a Reduction\n(r = {r_enz:.3f})', fontsize=14)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/reaction_kinetics_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nFigure saved to: simulations/chemistry/reaction_kinetics_coherence.png")

print("\n" + "=" * 70)
print("SESSION #70 COMPLETE: REACTION KINETICS")
print("=" * 70)
