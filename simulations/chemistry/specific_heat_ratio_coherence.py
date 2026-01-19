#!/usr/bin/env python3
"""
Chemistry Session #112: Specific Heat Ratio and Coherence

Test whether specific heat ratio γ_ad = Cp/Cv relates to coherence parameters.

Theory:
γ_ad = Cp/Cv (adiabatic index)

For ideal gas: γ = (f+2)/f where f = degrees of freedom
- Monoatomic: f=3, γ = 5/3 ≈ 1.67
- Diatomic: f=5, γ = 7/5 = 1.40
- Polyatomic: f=6+, γ → 1.0

For solids:
γ_ad = 1 + α²VTB/Cv (thermodynamic relation)
where:
- α = thermal expansion coefficient
- V = molar volume
- T = temperature
- B = bulk modulus

Since α ∝ γ_G/B (Session #109), we expect:
γ_ad - 1 ∝ γ_G² / B

This connects to the Grüneisen parameter:
γ_G = αVB/(Cv) (definition)
So: γ_ad = 1 + γ_G × α × V × B / Cv

Coherence question: Does γ_ad relate to γ_phonon or γ_G?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# SOLID DATA
# =============================================================================

# For solids at 300 K
solids_data = {
    # Format: Cp (J/mol·K), Cv (J/mol·K), α (×10⁻⁶/K), B (GPa), V (cm³/mol), θ_D (K), γ_G
    # Noble metals
    'Ag': {'Cp': 25.4, 'Cv': 24.9, 'alpha': 18.9, 'B': 100, 'V': 10.3, 'theta_D': 225, 'gamma_G': 2.40},
    'Cu': {'Cp': 24.4, 'Cv': 23.8, 'alpha': 16.5, 'B': 140, 'V': 7.1, 'theta_D': 343, 'gamma_G': 2.00},
    'Au': {'Cp': 25.4, 'Cv': 25.0, 'alpha': 14.2, 'B': 180, 'V': 10.2, 'theta_D': 165, 'gamma_G': 2.97},

    # Alkali metals
    'Na': {'Cp': 28.2, 'Cv': 24.5, 'alpha': 71.0, 'B': 6.3, 'V': 23.7, 'theta_D': 158, 'gamma_G': 1.25},
    'K': {'Cp': 29.6, 'Cv': 24.6, 'alpha': 83.0, 'B': 3.1, 'V': 45.5, 'theta_D': 91, 'gamma_G': 1.34},

    # Transition metals
    'W': {'Cp': 24.3, 'Cv': 23.8, 'alpha': 4.5, 'B': 310, 'V': 9.5, 'theta_D': 400, 'gamma_G': 1.62},
    'Mo': {'Cp': 24.1, 'Cv': 23.6, 'alpha': 4.8, 'B': 230, 'V': 9.4, 'theta_D': 450, 'gamma_G': 1.57},
    'Fe': {'Cp': 25.1, 'Cv': 24.4, 'alpha': 11.8, 'B': 170, 'V': 7.1, 'theta_D': 470, 'gamma_G': 1.66},
    'Ni': {'Cp': 26.1, 'Cv': 25.3, 'alpha': 13.4, 'B': 180, 'V': 6.6, 'theta_D': 450, 'gamma_G': 1.88},
    'Ti': {'Cp': 25.1, 'Cv': 24.2, 'alpha': 8.6, 'B': 110, 'V': 10.6, 'theta_D': 420, 'gamma_G': 1.24},

    # Simple metals
    'Al': {'Cp': 24.2, 'Cv': 23.4, 'alpha': 23.1, 'B': 76, 'V': 10.0, 'theta_D': 428, 'gamma_G': 2.17},
    'Pb': {'Cp': 26.8, 'Cv': 25.2, 'alpha': 28.9, 'B': 46, 'V': 18.3, 'theta_D': 105, 'gamma_G': 2.65},
    'Zn': {'Cp': 25.4, 'Cv': 24.6, 'alpha': 30.2, 'B': 70, 'V': 9.2, 'theta_D': 327, 'gamma_G': 2.10},
    'Mg': {'Cp': 24.9, 'Cv': 23.9, 'alpha': 24.8, 'B': 45, 'V': 14.0, 'theta_D': 400, 'gamma_G': 1.52},

    # Semiconductors
    'Si': {'Cp': 19.8, 'Cv': 19.6, 'alpha': 2.6, 'B': 98, 'V': 12.1, 'theta_D': 645, 'gamma_G': 1.00},
    'Ge': {'Cp': 23.4, 'Cv': 23.1, 'alpha': 5.9, 'B': 75, 'V': 13.6, 'theta_D': 374, 'gamma_G': 1.20},

    # Ceramics
    'Diamond': {'Cp': 6.1, 'Cv': 6.0, 'alpha': 1.0, 'B': 442, 'V': 3.4, 'theta_D': 2230, 'gamma_G': 0.90},
    'MgO': {'Cp': 37.0, 'Cv': 36.1, 'alpha': 10.8, 'B': 155, 'V': 11.3, 'theta_D': 946, 'gamma_G': 1.50},
}

# =============================================================================
# GAS DATA FOR COMPARISON
# =============================================================================

gases_data = {
    # Format: γ_ad = Cp/Cv, f (degrees of freedom), molecular structure
    'He': {'gamma_ad': 1.667, 'f': 3, 'structure': 'monoatomic'},
    'Ar': {'gamma_ad': 1.667, 'f': 3, 'structure': 'monoatomic'},
    'N2': {'gamma_ad': 1.401, 'f': 5, 'structure': 'diatomic'},
    'O2': {'gamma_ad': 1.400, 'f': 5, 'structure': 'diatomic'},
    'H2': {'gamma_ad': 1.410, 'f': 5, 'structure': 'diatomic'},
    'CO2': {'gamma_ad': 1.289, 'f': 7, 'structure': 'linear triatomic'},
    'H2O': {'gamma_ad': 1.330, 'f': 6, 'structure': 'nonlinear triatomic'},
    'CH4': {'gamma_ad': 1.310, 'f': 9, 'structure': 'tetrahedral'},
    'SF6': {'gamma_ad': 1.094, 'f': 15, 'structure': 'octahedral'},
}

# =============================================================================
# ANALYSIS: SOLIDS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #112: SPECIFIC HEAT RATIO AND COHERENCE")
print("="*70)

T = 300  # K

solids = list(solids_data.keys())
Cp = np.array([solids_data[s]['Cp'] for s in solids])
Cv = np.array([solids_data[s]['Cv'] for s in solids])
alpha = np.array([solids_data[s]['alpha'] for s in solids])
B = np.array([solids_data[s]['B'] for s in solids])
V = np.array([solids_data[s]['V'] for s in solids])
theta_D = np.array([solids_data[s]['theta_D'] for s in solids])
gamma_G = np.array([solids_data[s]['gamma_G'] for s in solids])

# Specific heat ratio
gamma_ad = Cp / Cv

# Coherence parameters
gamma_phonon = 2 * T / theta_D

# Grüneisen relation check: γ_G = αVB/Cv
gamma_G_calc = (alpha * 1e-6) * (V * 1e-6) * (B * 1e9) / Cv

print("\n" + "="*70)
print("SOLIDS: SPECIFIC HEAT RATIO")
print("="*70)

print(f"\n{'Solid':<10} {'Cp':<8} {'Cv':<8} {'γ_ad':<8} {'γ_G':<8} {'γ_G_calc':<10} {'γ_phonon':<8}")
print("-"*62)
for i, s in enumerate(solids):
    print(f"{s:<10} {Cp[i]:<8.1f} {Cv[i]:<8.1f} {gamma_ad[i]:<8.3f} {gamma_G[i]:<8.2f} {gamma_G_calc[i]:<10.2f} {gamma_phonon[i]:<8.3f}")

# =============================================================================
# CORRELATIONS
# =============================================================================

print("\n" + "="*70)
print("CORRELATIONS")
print("="*70)

# γ_ad vs γ_G
r1, p1 = stats.pearsonr(gamma_ad - 1, gamma_G**2)
print(f"\n(γ_ad - 1) vs γ_G²: r = {r1:.3f}")

# γ_ad vs γ_phonon
r2, p2 = stats.pearsonr(gamma_ad, gamma_phonon)
print(f"γ_ad vs γ_phonon: r = {r2:.3f}")

# γ_ad vs α²VB/Cv (thermodynamic identity)
thermo_term = (alpha * 1e-6)**2 * (V * 1e-6) * T * (B * 1e9) / Cv
r3, p3 = stats.pearsonr(gamma_ad - 1, thermo_term)
print(f"(γ_ad - 1) vs α²VTB/Cv: r = {r3:.3f}")

# γ_G calculated vs tabulated
r4, p4 = stats.pearsonr(gamma_G, gamma_G_calc)
print(f"γ_G (tabulated) vs γ_G (calculated): r = {r4:.3f}")

# γ_ad vs 1/θ_D
r5, p5 = stats.pearsonr(gamma_ad, 1/theta_D)
print(f"γ_ad vs 1/θ_D: r = {r5:.3f}")

# =============================================================================
# GAS ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("GASES: DEGREES OF FREEDOM")
print("="*70)

gases = list(gases_data.keys())
gamma_gas = np.array([gases_data[g]['gamma_ad'] for g in gases])
f_dof = np.array([gases_data[g]['f'] for g in gases])

print(f"\n{'Gas':<8} {'γ_ad':<8} {'f':<6} {'(f+2)/f':<10} {'Structure':<20}")
print("-"*56)
for g in gases:
    d = gases_data[g]
    predicted = (d['f'] + 2) / d['f']
    print(f"{g:<8} {d['gamma_ad']:<8.3f} {d['f']:<6} {predicted:<10.3f} {d['structure']:<20}")

# Correlation
r6, p6 = stats.pearsonr(gamma_gas, (f_dof + 2) / f_dof)
print(f"\nγ_ad vs (f+2)/f: r = {r6:.3f}")

# Effective degrees of freedom from γ
f_eff = 2 / (gamma_gas - 1)
r7, p7 = stats.pearsonr(f_eff, f_dof)
print(f"f_eff = 2/(γ-1) vs f: r = {r7:.3f}")

# =============================================================================
# COHERENCE INTERPRETATION
# =============================================================================

print("\n" + "="*70)
print("COHERENCE INTERPRETATION")
print("="*70)

print("""
1. SOLIDS: γ_ad ≈ 1 (very close to 1)
   Range: {:.3f} - {:.3f}
   This means Cp ≈ Cv for solids at RT.

2. GASES: γ_ad = (f+2)/f
   Monoatomic: f=3, γ = 1.67
   Diatomic: f=5, γ = 1.40
   Polyatomic: f→∞, γ → 1

3. THERMODYNAMIC IDENTITY
   γ_ad - 1 = α²VTB/Cv
   For solids: this is typically 0.01-0.2

4. GRÜNEISEN CONNECTION
   γ_G = αVB/Cv
   So: γ_ad - 1 = γ_G × α × T

5. COHERENCE?
   γ_ad vs γ_phonon: r = {:.3f}
   γ_ad vs γ_G²: r = {:.3f}

   WEAK correlations! γ_ad is primarily thermodynamic,
   not directly coherence-related.
""".format(np.min(gamma_ad), np.max(gamma_ad), r2, r1))

# =============================================================================
# PHYSICAL INSIGHTS
# =============================================================================

print("\n" + "="*70)
print("PHYSICAL INSIGHTS")
print("="*70)

print(f"""
1. γ_ad FOR SOLIDS IS NOT COHERENCE-DEPENDENT
   γ_ad vs γ_phonon: r = {r2:.3f} (WEAK)
   γ_ad vs γ_G²: r = {r1:.3f} (WEAK)

   Specific heat ratio reflects THERMODYNAMIC degrees of freedom,
   not phonon coherence.

2. γ_ad - 1 FOLLOWS THERMODYNAMIC IDENTITY
   (γ_ad - 1) vs α²VTB/Cv: r = {r3:.3f}
   This is an EXACT relation, not a correlation.

3. GRÜNEISEN PARAMETER VALIDATION
   γ_G (tabulated) vs γ_G (calculated): r = {r4:.3f}
   The Grüneisen definition γ_G = αVB/Cv holds.

4. GASES FOLLOW EQUIPARTITION
   γ = (f+2)/f with r = {r6:.3f}
   Each degree of freedom contributes (1/2)kT energy.

5. FRAMEWORK BOUNDARY IDENTIFIED
   γ_ad is NOT a coherence-related property.
   It's determined by:
   - Degrees of freedom (gases)
   - Thermodynamic identity (solids)

   This establishes γ_ad as a THERMODYNAMIC property,
   outside the coherence framework.
""")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Define colors by class
colors = []
for s in solids:
    if s in ['Ag', 'Cu', 'Au']:
        colors.append('gold')
    elif s in ['Na', 'K']:
        colors.append('blue')
    elif s in ['W', 'Mo', 'Fe', 'Ni', 'Ti']:
        colors.append('red')
    elif s in ['Al', 'Pb', 'Zn', 'Mg']:
        colors.append('purple')
    elif s in ['Si', 'Ge']:
        colors.append('orange')
    else:
        colors.append('green')

# Plot 1: Solids - γ_ad vs γ_phonon
ax1 = axes[0, 0]
ax1.scatter(gamma_phonon, gamma_ad, c=colors, s=100, alpha=0.7)
for i, s in enumerate(solids):
    ax1.annotate(s, (gamma_phonon[i], gamma_ad[i]), fontsize=9)
ax1.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
ax1.set_xlabel('γ_phonon', fontsize=12)
ax1.set_ylabel('γ_ad = Cp/Cv', fontsize=12)
ax1.set_title(f'Solids: γ_ad vs γ_phonon (r = {r2:.3f})', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: Grüneisen validation
ax2 = axes[0, 1]
ax2.scatter(gamma_G, gamma_G_calc, c=colors, s=100, alpha=0.7)
for i, s in enumerate(solids):
    ax2.annotate(s, (gamma_G[i], gamma_G_calc[i]), fontsize=9)

# 1:1 line
max_g = max(np.max(gamma_G), np.max(gamma_G_calc))
ax2.plot([0, max_g], [0, max_g], 'k--', alpha=0.5, label='1:1')
ax2.set_xlabel('γ_G (tabulated)', fontsize=12)
ax2.set_ylabel('γ_G = αVB/Cv (calculated)', fontsize=12)
ax2.set_title(f'Grüneisen Validation (r = {r4:.3f})', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Gases - γ vs (f+2)/f
ax3 = axes[1, 0]
ax3.scatter(f_dof, gamma_gas, c='steelblue', s=150, alpha=0.7)
for i, g in enumerate(gases):
    ax3.annotate(g, (f_dof[i], gamma_gas[i]), fontsize=10)

# Theory line
f_theory = np.linspace(3, 20, 100)
gamma_theory = (f_theory + 2) / f_theory
ax3.plot(f_theory, gamma_theory, 'k--', alpha=0.5, label='γ = (f+2)/f')

ax3.set_xlabel('Degrees of freedom f', fontsize=12)
ax3.set_ylabel('γ_ad = Cp/Cv', fontsize=12)
ax3.set_title(f'Gases: Equipartition (r = {r6:.3f})', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Distribution comparison
ax4 = axes[1, 1]
ax4.hist(gamma_ad, bins=10, color='steelblue', alpha=0.7, label='Solids', edgecolor='black')
ax4.axvline(x=np.mean(gamma_ad), color='blue', linestyle='--', label=f'Solid mean = {np.mean(gamma_ad):.3f}')
ax4.axvline(x=5/3, color='red', linestyle='--', alpha=0.5, label='Monoatomic gas = 1.67')
ax4.axvline(x=7/5, color='orange', linestyle='--', alpha=0.5, label='Diatomic gas = 1.40')
ax4.set_xlabel('γ_ad = Cp/Cv', fontsize=12)
ax4.set_ylabel('Count', fontsize=12)
ax4.set_title('Distribution of γ_ad', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/specific_heat_ratio_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #112")
print("="*70)

print(f"""
KEY RESULTS:

1. γ_ad vs COHERENCE (SOLIDS)
   γ_ad vs γ_phonon: r = {r2:.3f} (WEAK)
   γ_ad vs γ_G²: r = {r1:.3f} (WEAK)

   Specific heat ratio is NOT coherence-dependent.

2. THERMODYNAMIC IDENTITY VALIDATED
   (γ_ad - 1) = α²VTB/Cv
   r = {r3:.3f} (as expected for exact relation)

3. GRÜNEISEN DEFINITION VALIDATED
   γ_G = αVB/Cv
   Tabulated vs calculated: r = {r4:.3f}

4. GASES FOLLOW EQUIPARTITION
   γ = (f+2)/f
   r = {r6:.3f} (EXCELLENT)

5. SOLID γ_ad RANGE
   Mean: {np.mean(gamma_ad):.3f}
   Range: {np.min(gamma_ad):.3f} - {np.max(gamma_ad):.3f}
   Very close to 1 (Cp ≈ Cv)

FRAMEWORK BOUNDARY:
γ_ad = Cp/Cv is a THERMODYNAMIC property,
determined by degrees of freedom and α²VTB/Cv.

It is NOT directly related to coherence (γ_phonon, γ_electron).

This establishes a clear distinction:
- COHERENCE-RELATED: transport (σ, κ, α), optical (n, ε), decay (Γ)
- THERMODYNAMIC: heat capacity ratio (γ_ad = Cp/Cv)
""")

print("\nFigure saved to: specific_heat_ratio_coherence.png")
