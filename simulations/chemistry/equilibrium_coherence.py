"""
Chemistry Session #193: Chemical Equilibrium Constants Coherence
Testing equilibrium phenomena through γ ~ 1 framework

Key questions:
1. Is K = 1 a special γ ~ 1 condition?
2. Does ΔG° = 0 represent coherence balance?
3. How does Le Chatelier's principle relate to γ ~ 1?
4. Are reaction quotient Q/K dynamics coherence dynamics?
5. What makes equilibrium "balanced"?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("="*60)
print("CHEMISTRY SESSION #193: CHEMICAL EQUILIBRIUM COHERENCE")
print("="*60)

# Constants
R = 8.314  # J/(mol·K)
T = 298  # K
RT = R * T / 1000  # kJ/mol = 2.478 kJ/mol

# =============================================================================
# K = 1: THE BALANCED EQUILIBRIUM
# =============================================================================
print("\n" + "="*60)
print("1. K = 1: THE BALANCED EQUILIBRIUM (γ ~ 1)")
print("="*60)

# ΔG° = -RT ln(K)
# At K = 1: ΔG° = 0 (equally favored products/reactants)

print("\nThermodynamic relationship:")
print("  ΔG° = -RT × ln(K)")
print("  At K = 1: ΔG° = 0")
print()
print("This IS γ ~ 1:")
print("  Equal free energy of products and reactants")
print("  No thermodynamic preference")
print("  System at coherent balance")

# Calculate ΔG° for various K values
print("\nK vs ΔG° relationship:")
print("-"*40)
print(f"{'K':>12} {'ΔG° (kJ/mol)':>15}")
print("-"*40)
k_values = [1e-6, 1e-3, 0.1, 1.0, 10, 1000, 1e6]
for k in k_values:
    dg = -RT * np.log(k)
    print(f"{k:>12.0e} {dg:>15.2f}")

print("\nAt K = 1: ΔG° = 0.00 kJ/mol exactly!")

# =============================================================================
# EQUILIBRIUM CONSTANTS NEAR UNITY
# =============================================================================
print("\n" + "="*60)
print("2. EQUILIBRIUM CONSTANTS: DISTRIBUTION")
print("="*60)

# Real equilibrium constants from chemistry
equilibrium_data = {
    # Reaction: K at 25°C
    # Gas phase
    'N2 + 3H2 ⇌ 2NH3': 6.0e5,
    'H2 + I2 ⇌ 2HI': 54,
    'H2 + Cl2 ⇌ 2HCl': 4e31,
    'N2O4 ⇌ 2NO2': 4.6e-3,
    'PCl5 ⇌ PCl3 + Cl2': 1.8,
    'CO + H2O ⇌ CO2 + H2': 5.1,
    '2SO2 + O2 ⇌ 2SO3': 4.0e24,
    # Aqueous
    'CH3COOH ⇌ CH3COO- + H+': 1.8e-5,
    'NH3 + H2O ⇌ NH4+ + OH-': 1.8e-5,
    'H2CO3 ⇌ HCO3- + H+': 4.5e-7,
    'HF ⇌ H+ + F-': 6.8e-4,
    'CH3NH2 + H2O ⇌ CH3NH3+ + OH-': 4.4e-4,
    # Near K = 1 (interesting cases)
    'Fe3+ + SCN- ⇌ FeSCN2+': 1.1e2,
    'AgCl(s) + 2NH3 ⇌ Ag(NH3)2+ + Cl-': 2.0e3,
    'Cu2+ + 4NH3 ⇌ Cu(NH3)4(2+)': 1.1e13,
}

print("\nEquilibrium Constants at 25°C:")
print("-"*60)
print(f"{'Reaction':<40} {'K':>12} {'log10(K)':>10}")
print("-"*60)

log_k_values = []
for rxn, k in equilibrium_data.items():
    log_k = np.log10(k)
    print(f"{rxn:<40} {k:>12.1e} {log_k:>10.1f}")
    log_k_values.append(log_k)

log_k_arr = np.array(log_k_values)
print(f"\nRange: log10(K) from {min(log_k_arr):.1f} to {max(log_k_arr):.1f}")
print(f"Mean log10(K) = {np.mean(log_k_arr):.1f}")

# Count near K = 1 (|log K| < 2)
near_unity = np.sum(np.abs(log_k_arr) < 2)
print(f"\nReactions with |log10(K)| < 2: {near_unity}/{len(log_k_arr)}")
print("These are in the γ ~ 1 region!")

# =============================================================================
# REACTION QUOTIENT: Q/K DYNAMICS
# =============================================================================
print("\n" + "="*60)
print("3. REACTION QUOTIENT: γ_Q = Q/K")
print("="*60)

# Q = [products]/[reactants] (at any time)
# γ_Q = Q/K
# At γ = 1: equilibrium

print("\nReaction Quotient Dynamics:")
print("  γ_Q = Q/K")
print()
print("  γ < 1 (Q < K): forward reaction favored")
print("  γ = 1 (Q = K): EQUILIBRIUM")
print("  γ > 1 (Q > K): reverse reaction favored")
print()
print("This IS γ ~ 1 dynamics!")
print("System evolves toward γ = 1 (equilibrium)")

# Le Chatelier response
print("\nLe Chatelier's Principle:")
print("  Perturbation moves Q away from K (γ ≠ 1)")
print("  System responds to restore Q = K (γ = 1)")
print("  This IS coherence restoration!")

# =============================================================================
# TEMPERATURE DEPENDENCE: VAN'T HOFF
# =============================================================================
print("\n" + "="*60)
print("4. VAN'T HOFF: TEMPERATURE DEPENDENCE")
print("="*60)

# d(ln K)/dT = ΔH°/(RT²)
# Integrated: ln(K2/K1) = -ΔH°/R × (1/T2 - 1/T1)

# At what temperature does K = 1?
# 0 = -ΔH°/R × (1/T*) + ΔS°/R
# T* = ΔH°/ΔS° (temperature where K = 1)

print("\nVan't Hoff equation:")
print("  ln(K) = -ΔH°/(RT) + ΔS°/R")
print()
print("At K = 1 (ln K = 0):")
print("  T* = ΔH°/ΔS°")
print()
print("Every reaction has a temperature where K = 1!")
print("This is the γ ~ 1 temperature.")

# Example calculations
thermo_data = {
    # Reaction: (ΔH° kJ/mol, ΔS° J/mol·K)
    'N2O4 ⇌ 2NO2': (57.2, 176),
    'NH4Cl(s) ⇌ NH3 + HCl': (176, 285),
    'CaCO3 ⇌ Caite and CO2': (178, 161),
    'H2O(l) ⇌ H2O(g)': (44.0, 119),
    'PCl5 ⇌ PCl3 + Cl2': (92.5, 181),
}

print("\nTemperature where K = 1:")
print("-"*50)
print(f"{'Reaction':<30} {'T* (K)':>10} {'T* (°C)':>10}")
print("-"*50)

t_star_values = []
for rxn, (dh, ds) in thermo_data.items():
    if ds != 0:
        t_star = (dh * 1000) / ds  # K
        print(f"{rxn:<30} {t_star:>10.0f} {t_star-273:>10.0f}")
        t_star_values.append(t_star)

print("\nAt T = T*, the reaction is at γ ~ 1 (K = 1)!")

# =============================================================================
# GIBBS FREE ENERGY SURFACE
# =============================================================================
print("\n" + "="*60)
print("5. GIBBS FREE ENERGY: ΔG AND EQUILIBRIUM")
print("="*60)

# ΔG = ΔG° + RT ln(Q)
# At equilibrium: ΔG = 0, Q = K

print("\nGibbs Free Energy Equation:")
print("  ΔG = ΔG° + RT × ln(Q)")
print()
print("At equilibrium (ΔG = 0):")
print("  0 = ΔG° + RT × ln(K)")
print("  Therefore: ΔG° = -RT × ln(K)")
print()
print("ΔG = 0 IS the γ ~ 1 condition!")
print("No driving force for change")

# =============================================================================
# COUPLED EQUILIBRIA
# =============================================================================
print("\n" + "="*60)
print("6. COUPLED EQUILIBRIA: K_overall = K1 × K2")
print("="*60)

# For coupled reactions: K_overall = product of K values
# γ_overall = γ1 × γ2

print("\nCoupled Equilibria:")
print("  A ⇌ B (K1)")
print("  B ⇌ C (K2)")
print("  A ⇌ C (K_overall = K1 × K2)")
print()
print("In log space:")
print("  log(K_overall) = log(K1) + log(K2)")
print()
print("If K1 × K2 = 1, overall reaction at γ ~ 1!")

# Example: Buffer system
print("\nExample: Carbonate Buffer")
print("  H2CO3 ⇌ HCO3- + H+ (K1 = 4.5e-7)")
print("  HCO3- ⇌ CO3(2-) + H+ (K2 = 4.7e-11)")
print("  H2CO3 ⇌ CO3(2-) + 2H+ (K = K1×K2 = 2.1e-17)")

# =============================================================================
# EXTENT OF REACTION: ξ
# =============================================================================
print("\n" + "="*60)
print("7. EXTENT OF REACTION: ξ/ξ_max")
print("="*60)

# ξ = extent of reaction (moles converted)
# At equilibrium: ξ = ξ_eq
# γ_ξ = ξ/ξ_max

print("\nExtent of Reaction:")
print("  ξ = moles of reaction progress")
print("  ξ_eq = extent at equilibrium")
print()
print("For K = 1:")
print("  ξ_eq/ξ_max = 0.5 (50% conversion)")
print("  Halfway between reactants and products")
print("  This IS γ ~ 1!")

# Calculate equilibrium conversion for various K
print("\nEquilibrium Conversion (A ⇌ B, simple case):")
print("-"*40)
print(f"{'K':>10} {'α_eq':>12} {'% Conversion':>15}")
print("-"*40)

k_test = [0.01, 0.1, 1.0, 10, 100]
for k in k_test:
    # For A ⇌ B: K = α/(1-α), so α = K/(1+K)
    alpha = k / (1 + k)
    print(f"{k:>10.2f} {alpha:>12.3f} {alpha*100:>15.1f}%")

print("\nAt K = 1: α = 0.500 (50% conversion exactly!)")

# =============================================================================
# EQUILIBRIUM ISOTOPE EFFECTS
# =============================================================================
print("\n" + "="*60)
print("8. ISOTOPE EFFECTS: K_H/K_D")
print("="*60)

# Equilibrium isotope effect: K_H ≠ K_D
# γ_iso = K_H/K_D

isotope_data = {
    # Reaction: K_H/K_D
    'H2O + HD ⇌ HDO + H2': 1.85,
    'CH4 + D2O ⇌ CH3D + HDO': 1.40,
    'HCl(g) + D2O ⇌ DCl + HDO': 2.2,
    'NH3 + D2O ⇌ NH2D + HDO': 1.33,
    'C-H bond cleavage (typical)': 2.0,
    'O-H bond cleavage (typical)': 3.0,
}

print("\nEquilibrium Isotope Effects:")
print("-"*50)
print(f"{'Reaction':<40} {'K_H/K_D':>10}")
print("-"*50)

kie_values = []
for rxn, kie in isotope_data.items():
    print(f"{rxn:<40} {kie:>10.2f}")
    kie_values.append(kie)

kie_arr = np.array(kie_values)
print(f"\nMean K_H/K_D = {np.mean(kie_arr):.2f} ± {np.std(kie_arr):.2f}")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: EQUILIBRIUM COHERENCE PARAMETERS")
print("="*60)

summary = {
    'Reactions with |log K| < 2': (near_unity, len(log_k_arr)),
    'Mean log10(K)': (np.mean(log_k_arr), np.std(log_k_arr)),
    'Mean K_H/K_D': (np.mean(kie_arr), np.std(kie_arr)),
}

print(f"\n{'Parameter':<35} {'Value':>15}")
print("-"*55)
for param, val in summary.items():
    if isinstance(val[1], int):
        print(f"{param:<35} {val[0]}/{val[1]}")
    else:
        print(f"{param:<35} {val[0]:>7.2f} ± {val[1]:.2f}")

# Key γ ~ 1 conditions
print("\nKEY γ ~ 1 CONDITIONS IN CHEMICAL EQUILIBRIUM:")
print("1. K = 1: equal products and reactants (ΔG° = 0)")
print("2. Q/K = 1: system at equilibrium (ΔG = 0)")
print("3. ξ/ξ_max = 0.5: 50% conversion (for K = 1)")
print("4. T = T* = ΔH°/ΔS°: temperature where K = 1")
print("5. Le Chatelier: system restores γ ~ 1 after perturbation")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Chemistry Session #193: Chemical Equilibrium Coherence',
             fontsize=14, fontweight='bold')

# Panel 1: K vs ΔG° relationship
ax1 = axes[0, 0]
log_k_range = np.linspace(-10, 10, 100)
dg_values = -RT * log_k_range * np.log(10)
ax1.plot(log_k_range, dg_values, 'b-', linewidth=2)
ax1.axhline(y=0, color='red', linestyle='--', linewidth=2, label='ΔG° = 0')
ax1.axvline(x=0, color='red', linestyle='--', linewidth=2, label='K = 1')
ax1.plot(0, 0, 'ko', markersize=10, label='K=1, ΔG°=0 (γ~1)')
ax1.set_xlabel('log₁₀(K)')
ax1.set_ylabel('ΔG° (kJ/mol)')
ax1.set_title('ΔG° = -RT ln(K): K = 1 at ΔG° = 0')
ax1.legend()
ax1.set_xlim(-10, 10)
ax1.set_ylim(-60, 60)
ax1.grid(alpha=0.3)

# Panel 2: Distribution of log K
ax2 = axes[0, 1]
ax2.hist(log_k_arr, bins=10, color='steelblue', alpha=0.7, edgecolor='black')
ax2.axvline(x=0, color='red', linestyle='--', linewidth=2, label='K = 1 (γ ~ 1)')
ax2.axvspan(-2, 2, alpha=0.1, color='green', label='γ ~ 1 region')
ax2.set_xlabel('log₁₀(K)')
ax2.set_ylabel('Count')
ax2.set_title('Distribution of Equilibrium Constants')
ax2.legend()

# Panel 3: Q/K dynamics
ax3 = axes[1, 0]
q_k_range = np.logspace(-2, 2, 100)
dg_reaction = RT * np.log(q_k_range)  # ΔG when ΔG° = 0
ax3.semilogx(q_k_range, dg_reaction, 'b-', linewidth=2)
ax3.axhline(y=0, color='red', linestyle='--', linewidth=2, label='ΔG = 0 (equilibrium)')
ax3.axvline(x=1, color='red', linestyle='--', linewidth=2, label='Q/K = 1')
ax3.fill_between([0.01, 1], -10, 0, alpha=0.2, color='green', label='Forward favored')
ax3.fill_between([1, 100], 0, 10, alpha=0.2, color='orange', label='Reverse favored')
ax3.set_xlabel('Q/K (γ)')
ax3.set_ylabel('ΔG (kJ/mol)')
ax3.set_title('Q/K Dynamics: System Evolves to γ = 1')
ax3.legend(fontsize=8)
ax3.set_xlim(0.01, 100)
ax3.set_ylim(-15, 15)

# Panel 4: Conversion vs K
ax4 = axes[1, 1]
k_range = np.logspace(-2, 2, 100)
alpha_eq = k_range / (1 + k_range)
ax4.semilogx(k_range, alpha_eq, 'b-', linewidth=2)
ax4.axhline(y=0.5, color='red', linestyle='--', linewidth=2, label='50% conversion')
ax4.axvline(x=1, color='red', linestyle='--', linewidth=2, label='K = 1')
ax4.plot(1, 0.5, 'ko', markersize=10, label='K=1, α=0.5 (γ~1)')
ax4.set_xlabel('Equilibrium Constant K')
ax4.set_ylabel('Equilibrium Conversion α')
ax4.set_title('Conversion: α = 0.5 at K = 1')
ax4.legend()
ax4.set_ylim(0, 1)
ax4.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/equilibrium_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("\n" + "="*60)
print("FINDING #130: CHEMICAL EQUILIBRIUM AT γ ~ 1")
print("="*60)

print("""
KEY RESULTS:

1. K = 1 IS γ ~ 1
   - At K = 1: ΔG° = 0 (no thermodynamic preference)
   - Equal free energy of products and reactants
   - Perfect coherent balance

2. Q/K DYNAMICS
   - γ_Q = Q/K
   - System evolves toward γ = 1 (equilibrium)
   - Le Chatelier restores γ ~ 1 after perturbation

3. TEMPERATURE T* WHERE K = 1
   - T* = ΔH°/ΔS°
   - Every reaction has a γ ~ 1 temperature
   - At T*, products and reactants equally favored

4. 50% CONVERSION AT K = 1
   - α_eq = K/(1+K)
   - At K = 1: α = 0.5 exactly
   - Halfway point between reactants and products

5. GIBBS FREE ENERGY
   - ΔG = 0 at equilibrium
   - No driving force
   - This IS γ ~ 1

6. COUPLED REACTIONS
   - K_overall = K1 × K2
   - log(K) additive
   - Can engineer γ ~ 1 through coupling

PHYSICAL INSIGHT:
Chemical equilibrium IS a coherence phenomenon:
- K = 1 represents balanced forces (γ ~ 1)
- Q/K dynamics drive toward γ = 1
- Le Chatelier's principle IS coherence restoration
- Every reaction has a T* where γ ~ 1

The equilibrium constant K IS a coherence parameter:
- K >> 1: products strongly favored (γ >> 1)
- K << 1: reactants strongly favored (γ << 1)
- K = 1: coherent balance (γ = 1)

{}/{} reactions have |log10(K)| < 2 (near γ ~ 1).

56th phenomenon type at γ ~ 1!
""".format(near_unity, len(log_k_arr)))

print("\nVisualization saved to: equilibrium_coherence.png")
print("\nSESSION #193 COMPLETE")
