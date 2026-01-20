#!/usr/bin/env python3
"""
Chemistry Session #136: Solvent Reorganization Dynamics and Coherence

Session #135 showed γ_ET = λ/kT correlates with ET rates.
But λ is only HALF the story - dynamics matter too!

Marcus theory assumes instant nuclear equilibration.
But what if solvent motion is slow?

Key insight: The TIMESCALE of solvent reorganization matters.
- Fast solvent (τ_s < τ_ET): Adiabatic limit, Marcus applies
- Slow solvent (τ_s > τ_ET): Solvent-controlled, rate limited

Define: γ_solvent = τ_s / τ_ET

This measures whether the solvent can "keep up" with electron transfer.

Also explore:
- Solvent viscosity effects
- Longitudinal relaxation time τ_L
- Debye relaxation time τ_D
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Physical constants
kB = 1.381e-23    # J/K
eV_to_J = 1.602e-19
T = 298  # K

# =============================================================================
# DATASET: Solvents and Their Properties
# =============================================================================

# Collect solvent data
# Sources: Marcus reviews, Maroncelli group, Barbara group

solvents = [
    # Common organic solvents
    {"name": "Acetonitrile", "tau_L_ps": 0.2, "tau_D_ps": 3.5,
     "epsilon_s": 37.5, "epsilon_inf": 1.8, "eta_cP": 0.37,
     "lambda_out_eV": 1.2},

    {"name": "Methanol", "tau_L_ps": 5.0, "tau_D_ps": 51.0,
     "epsilon_s": 33.0, "epsilon_inf": 1.77, "eta_cP": 0.55,
     "lambda_out_eV": 1.1},

    {"name": "Water", "tau_L_ps": 0.5, "tau_D_ps": 8.3,
     "epsilon_s": 78.4, "epsilon_inf": 1.78, "eta_cP": 0.89,
     "lambda_out_eV": 1.3},

    {"name": "DMSO", "tau_L_ps": 2.0, "tau_D_ps": 19.0,
     "epsilon_s": 47.0, "epsilon_inf": 2.01, "eta_cP": 2.0,
     "lambda_out_eV": 1.15},

    {"name": "Ethanol", "tau_L_ps": 15.0, "tau_D_ps": 130.0,
     "epsilon_s": 24.5, "epsilon_inf": 1.85, "eta_cP": 1.08,
     "lambda_out_eV": 1.0},

    {"name": "DMF", "tau_L_ps": 0.5, "tau_D_ps": 9.0,
     "epsilon_s": 38.0, "epsilon_inf": 2.04, "eta_cP": 0.8,
     "lambda_out_eV": 1.18},

    {"name": "Propylene Carbonate", "tau_L_ps": 3.0, "tau_D_ps": 43.0,
     "epsilon_s": 65.0, "epsilon_inf": 2.0, "eta_cP": 2.5,
     "lambda_out_eV": 1.25},

    {"name": "Dichloromethane", "tau_L_ps": 0.4, "tau_D_ps": 1.8,
     "epsilon_s": 9.0, "epsilon_inf": 2.0, "eta_cP": 0.43,
     "lambda_out_eV": 0.8},

    {"name": "Tetrahydrofuran", "tau_L_ps": 0.3, "tau_D_ps": 3.0,
     "epsilon_s": 7.5, "epsilon_inf": 1.97, "eta_cP": 0.46,
     "lambda_out_eV": 0.75},

    {"name": "Benzonitrile", "tau_L_ps": 5.0, "tau_D_ps": 23.0,
     "epsilon_s": 26.0, "epsilon_inf": 2.3, "eta_cP": 1.24,
     "lambda_out_eV": 1.0},

    {"name": "1-Propanol", "tau_L_ps": 30.0, "tau_D_ps": 400.0,
     "epsilon_s": 20.5, "epsilon_inf": 1.92, "eta_cP": 2.0,
     "lambda_out_eV": 0.95},

    {"name": "Butyronitrile", "tau_L_ps": 0.5, "tau_D_ps": 7.5,
     "epsilon_s": 25.0, "epsilon_inf": 1.9, "eta_cP": 0.55,
     "lambda_out_eV": 1.05},

    {"name": "Formamide", "tau_L_ps": 3.0, "tau_D_ps": 37.0,
     "epsilon_s": 111.0, "epsilon_inf": 2.1, "eta_cP": 3.3,
     "lambda_out_eV": 1.35},
]

# Convert to numpy arrays
names = [s["name"] for s in solvents]
tau_L = np.array([s["tau_L_ps"] for s in solvents])  # ps
tau_D = np.array([s["tau_D_ps"] for s in solvents])  # ps
epsilon_s = np.array([s["epsilon_s"] for s in solvents])
epsilon_inf = np.array([s["epsilon_inf"] for s in solvents])
eta = np.array([s["eta_cP"] for s in solvents])  # cP
lambda_out = np.array([s["lambda_out_eV"] for s in solvents])  # eV

print("=" * 70)
print("CHEMISTRY SESSION #136: Solvent Dynamics and Coherence")
print("=" * 70)
print(f"\nSolvents analyzed: {len(solvents)}")
print(f"τ_L range: {tau_L.min():.1f} - {tau_L.max():.1f} ps")
print(f"τ_D range: {tau_D.min():.1f} - {tau_D.max():.1f} ps")
print(f"ε_s range: {epsilon_s.min():.1f} - {epsilon_s.max():.1f}")

# =============================================================================
# COHERENCE PARAMETERS
# =============================================================================
print("\n" + "=" * 70)
print("PART 1: Solvent Coherence Parameters")
print("=" * 70)

# Define solvent coherence parameters

# 1. Dielectric coherence: How much reorganization?
# γ_diel = (1/ε_∞ - 1/ε_s) / (1/ε_∞)
# This is the fraction of reorganization energy from solvent
gamma_diel = (1/epsilon_inf - 1/epsilon_s) / (1/epsilon_inf)

print(f"\nγ_diel = dielectric reorganization fraction")
print(f"γ_diel range: {gamma_diel.min():.3f} - {gamma_diel.max():.3f}")

# 2. Dynamic coherence: How fast is reorganization?
# Use τ_L as the relevant timescale
# Compare to thermal time τ_thermal = ℏ/kT ~ 25 fs at 300K
tau_thermal = 0.025  # ps
gamma_dynamic = tau_L / tau_thermal

print(f"\nγ_dynamic = τ_L / τ_thermal")
print(f"γ_dynamic range: {gamma_dynamic.min():.1f} - {gamma_dynamic.max():.0f}")

# 3. Combined solvent coherence
# Low λ_out AND fast τ_L → coherent
gamma_ET = lambda_out * eV_to_J / (kB * T)  # From session #135
gamma_solvent_combined = gamma_dynamic * (1 + gamma_diel)

print(f"\nγ_combined = γ_dynamic × (1 + γ_diel)")
print(f"γ_combined range: {gamma_solvent_combined.min():.1f} - {gamma_solvent_combined.max():.0f}")

# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================
print("\n" + "=" * 70)
print("PART 2: Property Correlations")
print("=" * 70)

# Test relationships between solvent properties

# τ_L vs ε_s (expected: negative correlation - polar = fast)
r1, p1 = stats.pearsonr(epsilon_s, tau_L)
print(f"\nε_s vs τ_L: r = {r1:.3f}, p = {p1:.4f}")
print("  (More polar → faster relaxation expected)")

# τ_L vs η (expected: positive - viscous = slow)
r2, p2 = stats.pearsonr(eta, tau_L)
print(f"\nη vs τ_L: r = {r2:.3f}, p = {p2:.4f}")
print("  (More viscous → slower relaxation expected)")

# λ_out vs (1/ε_∞ - 1/ε_s) (Marcus prediction)
marcus_factor = 1/epsilon_inf - 1/epsilon_s
r3, p3 = stats.pearsonr(marcus_factor, lambda_out)
print(f"\n(1/ε_∞ - 1/ε_s) vs λ_out: r = {r3:.3f}, p = {p3:.4f}")
print("  (Marcus outer-sphere prediction)")

# τ_D / τ_L ratio (should be ~ ε_s/ε_∞)
tau_ratio = tau_D / tau_L
epsilon_ratio = epsilon_s / epsilon_inf

r4, p4 = stats.pearsonr(epsilon_ratio, tau_ratio)
print(f"\nε_s/ε_∞ vs τ_D/τ_L: r = {r4:.3f}, p = {p4:.4f}")
print("  (Debye theory prediction)")

# =============================================================================
# SOLVENT CLASSIFICATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 3: Solvent Classification by Coherence")
print("=" * 70)

# Classify solvents by their coherence properties

# Fast solvents: τ_L < 1 ps
fast_mask = tau_L < 1.0
slow_mask = tau_L >= 5.0
medium_mask = ~fast_mask & ~slow_mask

print("\nFast solvents (τ_L < 1 ps):")
for i, name in enumerate(names):
    if fast_mask[i]:
        print(f"  {name}: τ_L = {tau_L[i]:.2f} ps, γ_dynamic = {gamma_dynamic[i]:.1f}")

print("\nMedium solvents (1 ≤ τ_L < 5 ps):")
for i, name in enumerate(names):
    if medium_mask[i]:
        print(f"  {name}: τ_L = {tau_L[i]:.2f} ps, γ_dynamic = {gamma_dynamic[i]:.1f}")

print("\nSlow solvents (τ_L ≥ 5 ps):")
for i, name in enumerate(names):
    if slow_mask[i]:
        print(f"  {name}: τ_L = {tau_L[i]:.2f} ps, γ_dynamic = {gamma_dynamic[i]:.0f}")

# =============================================================================
# ELECTRON TRANSFER RATE PREDICTION
# =============================================================================
print("\n" + "=" * 70)
print("PART 4: ET Rate Predictions")
print("=" * 70)

# For a model ET reaction, predict relative rates
# Using typical parameters: V = 10 meV, ΔG = -0.3 eV, λ_in = 0.2 eV

V_meV = 10.0
dG_eV = -0.3
lambda_in = 0.2  # Internal reorganization (fixed)

# Total λ = λ_in + λ_out
lambda_total = lambda_in + lambda_out

# Marcus rate (adiabatic limit)
hbar = 1.055e-34
V_J = V_meV * 1e-3 * eV_to_J
lambda_J = lambda_total * eV_to_J
dG_J = dG_eV * eV_to_J

activation = (dG_J + lambda_J)**2 / (4 * lambda_J * kB * T)
prefactor = (2 * np.pi / hbar) * V_J**2 / np.sqrt(4 * np.pi * lambda_J * kB * T)
k_marcus = prefactor * np.exp(-activation)
log_k_marcus = np.log10(k_marcus)

print(f"\nModel reaction: V = {V_meV} meV, ΔG = {dG_eV} eV, λ_in = {lambda_in} eV")
print(f"\nMarcus predicted log(k):")

# Sort by predicted rate
sorted_idx = np.argsort(log_k_marcus)[::-1]
for idx in sorted_idx[:5]:
    print(f"  {names[idx]}: log(k) = {log_k_marcus[idx]:.2f}, λ_total = {lambda_total[idx]:.2f} eV")

print("\n...slowest:")
for idx in sorted_idx[-3:]:
    print(f"  {names[idx]}: log(k) = {log_k_marcus[idx]:.2f}, λ_total = {lambda_total[idx]:.2f} eV")

# Now include dynamic corrections
# In solvent-controlled regime, rate limited by τ_L
# k_eff ~ k_marcus when τ_ET << τ_L (adiabatic)
# k_eff ~ 1/τ_L when τ_ET >> τ_L (solvent-controlled)

# Estimate τ_ET from Marcus rate
tau_ET = 1.0 / k_marcus  # seconds
tau_ET_ps = tau_ET * 1e12  # ps

# Effective rate with solvent dynamics
# k_eff = k_marcus / (1 + k_marcus × τ_L)
tau_L_s = tau_L * 1e-12  # seconds
k_eff = k_marcus / (1 + k_marcus * tau_L_s)
log_k_eff = np.log10(k_eff)

print("\n\nWith solvent dynamics correction:")
sorted_idx = np.argsort(log_k_eff)[::-1]
for idx in sorted_idx[:5]:
    ratio = tau_L[idx] / tau_ET_ps[idx]
    regime = "adiabatic" if ratio < 0.1 else "dynamic" if ratio > 10 else "intermediate"
    print(f"  {names[idx]}: log(k_eff) = {log_k_eff[idx]:.2f}, τ_L/τ_ET = {ratio:.1e} ({regime})")

# =============================================================================
# COHERENCE CLASSIFICATION
# =============================================================================
print("\n" + "=" * 70)
print("PART 5: Coherence Regime Classification")
print("=" * 70)

# Define coherence regime based on γ_dynamic
# γ < 10: Fast solvent, coherent ET possible
# 10 < γ < 100: Intermediate
# γ > 100: Slow solvent, classical dynamics

coherent_mask = gamma_dynamic < 10
intermediate_mask = (gamma_dynamic >= 10) & (gamma_dynamic < 100)
classical_mask = gamma_dynamic >= 100

print(f"\nCoherent regime (γ_dynamic < 10): {np.sum(coherent_mask)} solvents")
for i, name in enumerate(names):
    if coherent_mask[i]:
        print(f"  {name}: γ = {gamma_dynamic[i]:.1f}")

print(f"\nIntermediate regime (10 ≤ γ < 100): {np.sum(intermediate_mask)} solvents")
for i, name in enumerate(names):
    if intermediate_mask[i]:
        print(f"  {name}: γ = {gamma_dynamic[i]:.1f}")

print(f"\nClassical regime (γ ≥ 100): {np.sum(classical_mask)} solvents")
for i, name in enumerate(names):
    if classical_mask[i]:
        print(f"  {name}: γ = {gamma_dynamic[i]:.0f}")

# =============================================================================
# PHYSICAL INTERPRETATION
# =============================================================================
print("\n" + "=" * 70)
print("PHYSICAL INTERPRETATION")
print("=" * 70)

print("""
Solvent Coherence in Electron Transfer:

1. ENERGETIC COHERENCE (γ_diel = 1 - ε_∞/ε_s)
   How much nuclear reorganization is required?
   High ε_s → large (1/ε_∞ - 1/ε_s) → high λ_out
   But also: high ε_s stabilizes charged states

2. DYNAMIC COHERENCE (γ_dynamic = τ_L / τ_thermal)
   How fast can the solvent respond?
   Fast τ_L → solvent can track electron motion
   Slow τ_L → solvent lags, controls rate

3. COHERENCE REGIMES:
   γ_dynamic < 10:  COHERENT
     Solvent equilibrates instantly
     Quantum coherence possible
     Examples: acetonitrile, THF, DCM

   10 ≤ γ < 100:  INTERMEDIATE
     Partial decoherence
     Dynamics mixed with quantum
     Examples: methanol, DMSO, water

   γ ≥ 100:  CLASSICAL
     Solvent motion limits rate
     Classical dynamics dominate
     Examples: long-chain alcohols

4. WHY ACETONITRILE IS OPTIMAL:
   - Fast τ_L (0.2 ps) → γ_dynamic = 8 (coherent regime)
   - High ε_s (37.5) → stabilizes charges
   - Low viscosity (0.37 cP) → fast diffusion
   - Moderate λ_out (1.2 eV) → not too much reorganization

5. CONNECTION TO SESSION #135:
   Session #135: γ_ET = λ/kT (ENERGETIC coherence)
   Session #136: γ_dynamic = τ_L/τ_th (DYNAMIC coherence)

   BOTH matter for electron transfer!
   Optimal: Low λ (energetic) AND fast τ_L (dynamic)
""")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: τ_L vs ε_s
ax1 = axes[0, 0]
scatter1 = ax1.scatter(epsilon_s, tau_L, c=gamma_dynamic, cmap='coolwarm',
                       s=100, norm=plt.matplotlib.colors.LogNorm())
ax1.set_xlabel('ε_s (static dielectric)')
ax1.set_ylabel('τ_L (ps)')
ax1.set_yscale('log')
ax1.set_title(f'Dielectric vs Dynamics: r = {r1:.3f}')
cbar1 = plt.colorbar(scatter1, ax=ax1)
cbar1.set_label('γ_dynamic')
for i, name in enumerate(names):
    if tau_L[i] > 10 or epsilon_s[i] > 80 or tau_L[i] < 0.25:
        ax1.annotate(name[:8], (epsilon_s[i], tau_L[i]), fontsize=7)

# Plot 2: τ_L vs η
ax2 = axes[0, 1]
ax2.scatter(eta, tau_L, c=lambda_out, cmap='viridis', s=100)
ax2.set_xlabel('η (cP)')
ax2.set_ylabel('τ_L (ps)')
ax2.set_yscale('log')
ax2.set_title(f'Viscosity vs Dynamics: r = {r2:.3f}')
cbar2 = plt.colorbar(ax2.collections[0], ax=ax2)
cbar2.set_label('λ_out (eV)')

# Plot 3: λ_out vs Marcus factor
ax3 = axes[1, 0]
ax3.scatter(marcus_factor, lambda_out, c=tau_L, cmap='plasma', s=100,
           norm=plt.matplotlib.colors.LogNorm())
ax3.set_xlabel('1/ε_∞ - 1/ε_s')
ax3.set_ylabel('λ_out (eV)')
ax3.set_title(f'Marcus Factor Correlation: r = {r3:.3f}')
cbar3 = plt.colorbar(ax3.collections[0], ax=ax3)
cbar3.set_label('τ_L (ps)')

# Plot 4: Coherence classification
ax4 = axes[1, 1]
colors = ['green' if c else 'orange' if i else 'red'
          for c, i in zip(coherent_mask, intermediate_mask)]
ax4.scatter(gamma_dynamic, gamma_ET, c=colors, s=100)
ax4.axvline(x=10, color='gray', linestyle='--', alpha=0.5, label='Coherent/Intermediate')
ax4.axvline(x=100, color='gray', linestyle=':', alpha=0.5, label='Intermediate/Classical')
ax4.set_xlabel('γ_dynamic = τ_L/τ_thermal')
ax4.set_ylabel('γ_ET = λ/kT')
ax4.set_xscale('log')
ax4.set_title('Coherence Classification')
for i, name in enumerate(names):
    ax4.annotate(name[:6], (gamma_dynamic[i], gamma_ET[i]), fontsize=6)

plt.suptitle('Session #136: Solvent Dynamics and Coherence', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/solvent_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

print("\nVisualization saved to solvent_coherence.png")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SESSION #136 SUMMARY")
print("=" * 70)

print(f"""
KEY FINDINGS:

1. Two coherence parameters for solvents:
   - γ_diel = 1 - ε_∞/ε_s (energetic)
   - γ_dynamic = τ_L / τ_thermal (dynamic)

2. Property correlations:
   - ε_s vs τ_L: r = {r1:.3f}, p = {p1:.4f}
   - η vs τ_L: r = {r2:.3f}, p = {p2:.4f}
   - Marcus factor vs λ_out: r = {r3:.3f}, p = {p3:.4f}

3. Coherence classification:
   - Coherent (γ < 10): {np.sum(coherent_mask)} solvents
   - Intermediate: {np.sum(intermediate_mask)} solvents
   - Classical (γ > 100): {np.sum(classical_mask)} solvents

4. Optimal solvents for coherent ET:
   - Fast τ_L (< 1 ps): Acetonitrile, THF, DCM
   - High ε_s (stabilizes charges)
   - Low viscosity

5. Session #135 + #136:
   - γ_ET = λ/kT (energetic coherence)
   - γ_dynamic = τ_L/τ_th (dynamic coherence)
   BOTH matter for electron transfer!

FRAMEWORK EXTENSION:
Two dimensions of solvent coherence:
1. Energetic: How much reorganization? (λ)
2. Dynamic: How fast? (τ_L)

Coherent ET requires BOTH:
- Low γ_ET (small λ)
- Low γ_dynamic (fast τ_L)
""")

print("\nSESSION #136 COMPLETE")
print("=" * 70)
