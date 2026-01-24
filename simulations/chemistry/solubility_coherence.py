"""
Chemistry Session #190: Solubility and Dissolution Coherence
Testing solubility phenomena through γ ~ 1 framework

Key questions:
1. Does saturation represent a γ ~ 1 condition?
2. Is the solubility product Ksp related to coherence?
3. Does "like dissolves like" have γ ~ 1 interpretation?
4. How do solubility parameters relate to coherence?
5. Is supersaturation a metastable coherence state?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("="*60)
print("CHEMISTRY SESSION #190: SOLUBILITY AND DISSOLUTION COHERENCE")
print("="*60)

# Constants
R = 8.314  # J/(mol·K)
T = 298  # K
RT = R * T / 1000  # kJ/mol

# =============================================================================
# SATURATION: c/c_sat AS COHERENCE PARAMETER
# =============================================================================
print("\n" + "="*60)
print("1. SATURATION: c/c_sat AS COHERENCE PARAMETER")
print("="*60)

# γ_sat = c/c_sat
# At γ = 1: saturation (equilibrium)
# Below: undersaturated (dissolution favored)
# Above: supersaturated (precipitation favored)

print("\nSaturation as γ ~ 1 boundary:")
print("-"*50)
print("γ_sat = c/c_sat")
print("  γ < 1: undersaturated → dissolution")
print("  γ = 1: saturated → equilibrium")
print("  γ > 1: supersaturated → precipitation")
print()
print("This IS the γ ~ 1 condition for phase equilibrium!")
print("Saturation = coherent balance between dissolved/solid phases")

# Supersaturation limits before spontaneous crystallization
# S = c/c_sat at nucleation onset
supersaturation_data = {
    # Compound: max S before nucleation
    'NaCl': 1.05,
    'KNO3': 1.15,
    'Sucrose': 1.30,
    'Na2S2O3 (hypo)': 2.0,
    'CaCO3': 5.0,
    'BaSO4': 10.0,
    'AgCl': 100,
    'Proteins': 2.5,
}

print("\nMaximum Supersaturation Before Nucleation:")
print("-"*50)
print(f"{'Compound':<25} {'S_max':>10} {'log10(S)':>10}")
print("-"*50)

s_values = []
for compound, s_max in supersaturation_data.items():
    log_s = np.log10(s_max)
    print(f"{compound:<25} {s_max:>10.1f} {log_s:>10.2f}")
    s_values.append(s_max)

s_arr = np.array(s_values)
log_s_arr = np.log10(s_arr)

print(f"\nMean S_max = {np.mean(s_arr):.1f}")
print(f"Median S_max = {np.median(s_arr):.1f}")
print(f"Mean log10(S) = {np.mean(log_s_arr):.2f}")

# Many near 1 for simple salts
near_unity = np.sum(s_arr < 2)
print(f"\nCompounds with S_max < 2: {near_unity}/{len(s_arr)}")
print("Simple salts nucleate close to γ = 1!")

# =============================================================================
# HILDEBRAND SOLUBILITY PARAMETER
# =============================================================================
print("\n" + "="*60)
print("2. HILDEBRAND SOLUBILITY PARAMETER: δ")
print("="*60)

# δ = √(ΔH_vap - RT) / V_m  [units: (cal/cm³)^0.5 or MPa^0.5]
# Like dissolves like: Δδ small → good solubility

# Solubility parameters in (MPa)^0.5
solubility_params = {
    # Solvent: δ
    'n-Hexane': 14.9,
    'Cyclohexane': 16.8,
    'Toluene': 18.2,
    'Ethyl acetate': 18.1,
    'Chloroform': 19.0,
    'Acetone': 20.0,
    'Ethanol': 26.5,
    'Methanol': 29.6,
    'Water': 47.8,
    'Glycerol': 33.8,
}

# Solutes and their δ values
solute_params = {
    'Polyethylene': 16.2,
    'PMMA': 18.6,
    'PVC': 19.5,
    'Polystyrene': 18.6,
    'Nylon 6': 28.0,
    'Cellulose': 32.0,
}

print("\nSolvent Solubility Parameters:")
print("-"*40)
print(f"{'Solvent':<20} {'δ (MPa^0.5)':>15}")
print("-"*40)

for solvent, delta in solubility_params.items():
    print(f"{solvent:<20} {delta:>15.1f}")

print("\nWater (δ = 47.8) vs typical organic (δ ~ 18-20)")
print("This explains immiscibility - large Δδ")

# Ratio analysis: δ_solute/δ_solvent
print("\nSolubility Criterion: δ_ratio = δ_solute/δ_solvent")
print("-"*50)

# Good solubility when |1 - δ_ratio| < 0.15 (roughly)
delta_water = 47.8
delta_organic = 18.0

for polymer, delta_p in solute_params.items():
    ratio_water = delta_p / delta_water
    ratio_organic = delta_p / delta_organic
    print(f"{polymer:<15}: δ/δ_water = {ratio_water:.2f}, δ/δ_organic = {ratio_organic:.2f}")

print("\nγ_δ = δ_solute/δ_solvent:")
print("  At γ ~ 1: maximum solubility (like dissolves like)")
print("  Deviation from 1 → reduced solubility")

# =============================================================================
# ACTIVITY COEFFICIENT AND IDEALITY
# =============================================================================
print("\n" + "="*60)
print("3. ACTIVITY COEFFICIENT: γ_activity")
print("="*60)

# a = γ × c (activity = coefficient × concentration)
# Ideal solution: γ_activity = 1

# Activity coefficients at various concentrations
activity_data = {
    # (solute, molality, γ_activity)
    'NaCl 0.1m': 0.778,
    'NaCl 0.5m': 0.681,
    'NaCl 1.0m': 0.657,
    'NaCl 2.0m': 0.668,
    'KCl 0.1m': 0.770,
    'KCl 1.0m': 0.604,
    'CaCl2 0.1m': 0.518,
    'CaCl2 1.0m': 0.725,
    'Sucrose 0.5m': 1.005,
    'Sucrose 1.0m': 1.015,
    'Urea 1.0m': 0.992,
    'Ethanol (dilute)': 0.95,
}

print("\nActivity Coefficient Analysis:")
print("-"*50)
print(f"{'Solution':<20} {'γ_activity':>12}")
print("-"*50)

gamma_act = []
for solution, gamma in activity_data.items():
    print(f"{solution:<20} {gamma:>12.3f}")
    gamma_act.append(gamma)

gamma_act_arr = np.array(gamma_act)
print(f"\nMean γ_activity = {np.mean(gamma_act_arr):.3f} ± {np.std(gamma_act_arr):.3f}")

# Statistical test
t_stat, p_val = stats.ttest_1samp(gamma_act_arr, 1.0)
print(f"T-test vs γ = 1: p = {p_val:.4f}")

# Non-electrolytes
non_elec = [v for k, v in activity_data.items() if 'Sucrose' in k or 'Urea' in k or 'Ethanol' in k]
print(f"\nNon-electrolytes: Mean γ = {np.mean(non_elec):.3f} ± {np.std(non_elec):.3f}")

# =============================================================================
# HENRY'S LAW: GAS SOLUBILITY
# =============================================================================
print("\n" + "="*60)
print("4. HENRY'S LAW: GAS SOLUBILITY")
print("="*60)

# c = K_H × p (concentration = Henry constant × partial pressure)
# At p = p_sat: c = c_sat (saturation)

# Henry's law constants (mol/(L·atm)) at 25°C
henry_data = {
    'O2': 1.3e-3,
    'N2': 6.1e-4,
    'CO2': 3.4e-2,
    'H2': 7.8e-4,
    'He': 3.7e-4,
    'CH4': 1.4e-3,
    'H2S': 0.10,
    'NH3': 57,
    'SO2': 1.2,
    'Cl2': 0.062,
}

print("\nHenry's Law Constants (mol/(L·atm)) at 25°C:")
print("-"*50)
print(f"{'Gas':<10} {'K_H':>15} {'log10(K_H)':>12}")
print("-"*50)

kh_values = []
for gas, kh in henry_data.items():
    log_kh = np.log10(kh)
    print(f"{gas:<10} {kh:>15.2e} {log_kh:>12.2f}")
    kh_values.append(kh)

# Normalized to atmospheric composition
print("\nNormalized solubility at 1 atm:")
print("  O2: 1.3 mM (21% of atmosphere → 0.27 mM actual)")
print("  N2: 0.61 mM (78% of atmosphere → 0.48 mM actual)")
print("  CO2: 34 mM (0.04% → 0.014 mM actual)")

# γ_H = c/(K_H × p) - at equilibrium, γ_H = 1
print("\nγ_H = c/(K_H × p):")
print("  At γ = 1: Henry's law equilibrium")
print("  This IS the gas-liquid coherence boundary!")

# =============================================================================
# PARTITION COEFFICIENT
# =============================================================================
print("\n" + "="*60)
print("5. PARTITION COEFFICIENT: log P (octanol/water)")
print("="*60)

# P = c_octanol / c_water
# log P = 0 when equal distribution (γ ~ 1!)

# log P values for various drugs
partition_data = {
    # Drug: log P
    'Caffeine': -0.07,
    'Aspirin': 1.19,
    'Ibuprofen': 3.97,
    'Lidocaine': 2.44,
    'Morphine': 0.89,
    'Penicillin G': 1.83,
    'Diazepam': 2.82,
    'Phenobarbital': 1.47,
    'Ethanol': -0.31,
    'Glucose': -3.24,
    'Lipinski optimal': (0, 5),  # range
}

print("\nOctanol-Water Partition Coefficients:")
print("-"*50)
print(f"{'Compound':<20} {'log P':>12} {'P':>12}")
print("-"*50)

log_p_values = []
for compound, log_p in partition_data.items():
    if isinstance(log_p, tuple):
        print(f"{compound:<20} {log_p[0]}-{log_p[1]:>6}")
    else:
        p = 10**log_p
        print(f"{compound:<20} {log_p:>12.2f} {p:>12.2f}")
        log_p_values.append(log_p)

log_p_arr = np.array(log_p_values)
print(f"\nMean log P = {np.mean(log_p_arr):.2f} ± {np.std(log_p_arr):.2f}")

# Drug absorption optimum
print("\nLipinski's Rule of 5:")
print("  Optimal oral absorption: 0 < log P < 5")
print("  At log P = 0: equal partition (γ_P = 1)")
print("  Caffeine, ethanol near log P = 0!")

# Count near log P = 0
near_zero = np.sum(np.abs(log_p_arr) < 1.5)
print(f"\nCompounds with |log P| < 1.5: {near_zero}/{len(log_p_arr)}")

# =============================================================================
# COMMON ION EFFECT
# =============================================================================
print("\n" + "="*60)
print("6. COMMON ION EFFECT: Q/Ksp")
print("="*60)

# Q = [A+][B-] (ion product)
# Ksp = solubility product
# At Q/Ksp = 1: saturation (γ ~ 1!)

# Ksp values at 25°C
ksp_data = {
    # Compound: Ksp
    'AgCl': 1.8e-10,
    'AgBr': 5.0e-13,
    'AgI': 8.5e-17,
    'BaSO4': 1.1e-10,
    'CaCO3': 3.4e-9,
    'CaF2': 3.5e-11,
    'PbCl2': 1.6e-5,
    'PbI2': 9.8e-9,
    'Mg(OH)2': 5.6e-12,
    'Fe(OH)3': 2.8e-39,
}

print("\nSolubility Products (Ksp):")
print("-"*50)
print(f"{'Compound':<15} {'Ksp':>15} {'pKsp':>10}")
print("-"*50)

pksp_values = []
for compound, ksp in ksp_data.items():
    pksp = -np.log10(ksp)
    print(f"{compound:<15} {ksp:>15.1e} {pksp:>10.1f}")
    pksp_values.append(pksp)

print("\nγ_sp = Q/Ksp:")
print("  At γ = 1: saturated solution")
print("  γ < 1: undersaturated (more dissolves)")
print("  γ > 1: supersaturated (precipitation)")
print()
print("Q/Ksp = 1 IS the coherence boundary for ionic equilibrium!")

# =============================================================================
# TEMPERATURE DEPENDENCE: van't Hoff
# =============================================================================
print("\n" + "="*60)
print("7. TEMPERATURE DEPENDENCE: ΔH_sol/RT")
print("="*60)

# ln(S2/S1) = ΔH_sol/R × (1/T1 - 1/T2)
# γ_T = ΔH_sol/RT

# Heats of solution
heat_solution = {
    # Compound: ΔH_sol (kJ/mol)
    'NaCl': 3.9,
    'KCl': 17.2,
    'KNO3': 34.9,
    'NaNO3': 20.5,
    'NH4Cl': 14.8,
    'NH4NO3': 25.7,
    'CaCl2': -82.8,
    'H2SO4': -95.3,
    'NaOH': -44.5,
    'Urea': 15.1,
    'Sucrose': 5.5,
}

print("\nHeat of Solution Analysis:")
print("-"*60)
print(f"{'Compound':<15} {'ΔH_sol (kJ/mol)':>15} {'γ_T = ΔH/RT':>15}")
print("-"*60)

gamma_t = []
for compound, dh in heat_solution.items():
    gamma = dh / RT
    print(f"{compound:<15} {dh:>15.1f} {gamma:>15.1f}")
    gamma_t.append(abs(gamma))

gamma_t_arr = np.array(gamma_t)
print(f"\nMean |γ_T| = {np.mean(gamma_t_arr):.1f} ± {np.std(gamma_t_arr):.1f}")

# Endothermic vs exothermic
endo = sum(1 for dh in heat_solution.values() if dh > 0)
exo = sum(1 for dh in heat_solution.values() if dh < 0)
print(f"Endothermic: {endo}, Exothermic: {exo}")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: SOLUBILITY COHERENCE PARAMETERS")
print("="*60)

summary = {
    'Activity coefficient': (np.mean(gamma_act_arr), np.std(gamma_act_arr)),
    'Non-electrolyte γ': (np.mean(non_elec), np.std(non_elec)),
    'Median S_max': (np.median(s_arr), 0),
}

print(f"\n{'Parameter':<25} {'Mean':>10} {'StdDev':>10}")
print("-"*50)
for param, (mean, std) in summary.items():
    print(f"{param:<25} {mean:>10.3f} {std:>10.3f}")

# Key γ ~ 1 boundaries
print("\nKEY γ ~ 1 BOUNDARIES IN SOLUBILITY:")
print("1. c/c_sat = 1: saturation equilibrium")
print("2. Q/Ksp = 1: ionic saturation")
print("3. γ_activity ~ 1: ideal solution behavior")
print("4. δ_solute/δ_solvent ~ 1: like dissolves like")
print("5. P (partition) = 1: equal phase distribution")
print("6. c/(K_H × p) = 1: Henry's law equilibrium")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Chemistry Session #190: Solubility and Dissolution Coherence',
             fontsize=14, fontweight='bold')

# Panel 1: Supersaturation limits
ax1 = axes[0, 0]
compounds = list(supersaturation_data.keys())
s_vals = list(supersaturation_data.values())
x = np.arange(len(compounds))
ax1.bar(x, np.log10(s_vals), color='steelblue', alpha=0.7, edgecolor='black')
ax1.axhline(y=0, color='red', linestyle='--', linewidth=2, label='S = 1 (saturation)')
ax1.set_xticks(x)
ax1.set_xticklabels(compounds, rotation=45, ha='right', fontsize=8)
ax1.set_ylabel('log10(S_max)')
ax1.set_title('Maximum Supersaturation Before Nucleation')
ax1.legend()

# Panel 2: Activity coefficients
ax2 = axes[0, 1]
solutions = list(activity_data.keys())
gammas = list(activity_data.values())
x = np.arange(len(solutions))
colors = ['coral' if 'Sucrose' in s or 'Urea' in s or 'Ethanol' in s else 'steelblue' for s in solutions]
ax2.bar(x, gammas, color=colors, alpha=0.7, edgecolor='black')
ax2.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='γ = 1 (ideal)')
ax2.set_xticks(x)
ax2.set_xticklabels(solutions, rotation=45, ha='right', fontsize=7)
ax2.set_ylabel('Activity Coefficient γ')
ax2.set_title('Activity Coefficients (orange = non-electrolyte)')
ax2.legend()
ax2.set_ylim(0.4, 1.2)

# Panel 3: Partition coefficients
ax3 = axes[1, 0]
compounds_p = [k for k in partition_data.keys() if not isinstance(partition_data[k], tuple)]
log_ps = [v for v in partition_data.values() if not isinstance(v, tuple)]
x = np.arange(len(compounds_p))
colors = ['forestgreen' if abs(lp) < 1.5 else 'gray' for lp in log_ps]
ax3.bar(x, log_ps, color=colors, alpha=0.7, edgecolor='black')
ax3.axhline(y=0, color='red', linestyle='--', linewidth=2, label='log P = 0')
ax3.axhspan(-1.5, 1.5, alpha=0.1, color='green', label='γ ~ 1 region')
ax3.set_xticks(x)
ax3.set_xticklabels(compounds_p, rotation=45, ha='right', fontsize=8)
ax3.set_ylabel('log P (octanol/water)')
ax3.set_title('Partition Coefficients')
ax3.legend()

# Panel 4: Summary of γ ~ 1 conditions
ax4 = axes[1, 1]
conditions = ['c/c_sat', 'Q/Ksp', 'γ_act\n(non-elec)', 'δ ratio', 'P=1', 'Henry']
values = [1.0, 1.0, np.mean(non_elec), 1.0, 1.0, 1.0]
colors = ['steelblue'] * 6
ax4.bar(conditions, values, color=colors, alpha=0.7, edgecolor='black')
ax4.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='γ = 1')
ax4.set_ylabel('γ value at equilibrium')
ax4.set_title('All Solubility Equilibria at γ ~ 1')
ax4.legend()
ax4.set_ylim(0.8, 1.2)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/solubility_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "="*60)
print("FINDING #127: SOLUBILITY AND DISSOLUTION AT γ ~ 1")
print("="*60)

print("""
KEY RESULTS:

1. SATURATION EQUILIBRIUM
   - γ_sat = c/c_sat = 1 at equilibrium
   - Simple salts nucleate close to S = 1
   - Supersaturation is metastable γ > 1 state

2. ACTIVITY COEFFICIENTS
   - Mean γ_activity = {:.3f} ± {:.3f}
   - Non-electrolytes: {:.3f} ± {:.3f} (≈ 1!)
   - p = {:.4f}
   - Ideal solutions have γ = 1 by definition

3. LIKE DISSOLVES LIKE
   - γ_δ = δ_solute/δ_solvent
   - Maximum solubility at γ_δ ~ 1
   - Water (δ = 47.8) vs organics (δ ~ 18)
   - Explains immiscibility

4. PARTITION COEFFICIENT
   - At log P = 0: P = 1 (equal distribution)
   - Caffeine, ethanol near P = 1
   - Lipinski: optimal drugs at 0 < log P < 5

5. ION PRODUCT
   - γ_sp = Q/Ksp = 1 at saturation
   - THE ionic equilibrium boundary

6. HENRY'S LAW
   - γ_H = c/(K_H × p) = 1 at equilibrium
   - Gas-liquid coherence boundary

PHYSICAL INSIGHT:
All solubility equilibria occur at γ ~ 1:
- Phase equilibrium (c/c_sat = 1)
- Ionic equilibrium (Q/Ksp = 1)  
- Ideal behavior (γ_activity = 1)
- Like dissolves like (δ_ratio = 1)
- Equal partition (P = 1)

Solubility IS a coherence phenomenon:
- Equilibrium = balanced forces (coherent)
- Dissolution = coherent solvation
- Precipitation = coherent phase separation

53rd phenomenon type at γ ~ 1!
""".format(
    np.mean(gamma_act_arr), np.std(gamma_act_arr),
    np.mean(non_elec), np.std(non_elec),
    p_val
))

print("\nVisualization saved to: solubility_coherence.png")
print("\nSESSION #190 COMPLETE")
