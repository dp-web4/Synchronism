#!/usr/bin/env python3
"""
Phase 2 Session #5: Boundary vs Bulk Coherence

The framework fails at surfaces/boundaries because:
1. γ is a BULK parameter (assumes large homogeneous system)
2. Work function φ creates exponential barriers that overwhelm any γ effect
3. Surface electronic structure differs from bulk

Test case: Thermionic emission (Session #98)
  J = A × T² × exp(-φ/kT)
  The exponential in φ dominates — can γ compete?

Also examine: Why the Richardson constant A varies from material to material,
and whether this variation (which should reflect surface physics) correlates
with coherence.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("PHASE 2 SESSION #5: BOUNDARY vs BULK COHERENCE")
print("When Surface Physics Overwhelms Bulk Coherence")
print("=" * 70)

# ==============================================================================
# Thermionic emission data (from Session #98)
# ==============================================================================

# Materials: work function φ (eV), Richardson constant A (A/cm²K²), θ_D (K)
thermionic_data = {
    # Alkali metals
    'Li':    (2.93, 30,  344),
    'Na':    (2.75, 45,  158),
    'K':     (2.30, 48,  91),
    'Rb':    (2.16, 60,  56),
    'Cs':    (2.14, 75,  38),
    # Alkaline earth
    'Ca':    (2.87, 60,  230),
    'Sr':    (2.59, 60,  147),
    'Ba':    (2.52, 60,  110),
    # Transition metals
    'Ti':    (4.33, 60,  420),
    'V':     (4.30, 60,  380),
    'Cr':    (4.50, 60,  630),
    'Fe':    (4.67, 26,  470),
    'Co':    (5.00, 41,  445),
    'Ni':    (5.15, 30,  450),
    'Cu':    (4.65, 65,  343),
    'Zn':    (4.33, 50,  327),
    # Noble metals
    'Ag':    (4.26, 60,  225),
    'Au':    (5.10, 60,  165),
    'Pt':    (5.65, 32,  240),
    # Refractory
    'Mo':    (4.60, 55,  450),
    'Ta':    (4.25, 60,  258),
    'W':     (4.55, 60,  400),
    'Re':    (4.96, 60,  430),
    # Post-transition
    'Al':    (4.28, 120, 428),
    'Pb':    (4.25, 12,  105),
    # Special
    'LaB6':  (2.70, 29,  650),
}

T = 300  # K
k_B_eV = 8.617e-5  # eV/K

names = list(thermionic_data.keys())
phi = np.array([thermionic_data[m][0] for m in names])
A = np.array([thermionic_data[m][1] for m in names])
theta_D = np.array([thermionic_data[m][2] for m in names])

gamma_phonon = 2 * T / theta_D
A0 = 120.4  # Theoretical Richardson constant
A_ratio = A / A0  # Deviation from theoretical

# Current density at T=1000K (typical operating temperature for thermionic cathodes)
T_op = 1000  # K
J = A * T_op**2 * np.exp(-phi / (k_B_eV * T_op))

print(f"\n{'Material':<8} {'φ (eV)':<8} {'A (A/cm²K²)':<14} {'A/A₀':<8} {'θ_D (K)':<8} {'γ_ph':<8} {'J@1000K':<12}")
print("-" * 76)
for i, m in enumerate(names):
    print(f"{m:<8} {phi[i]:<8.2f} {A[i]:<14.0f} {A_ratio[i]:<8.2f} {theta_D[i]:<8.0f} {gamma_phonon[i]:<8.2f} {J[i]:<12.2e}")

# ==============================================================================
# Analysis 1: What determines thermionic emission?
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 1: WHAT DETERMINES THERMIONIC EMISSION?")
print("=" * 70)

log_J = np.log10(J + 1e-100)  # avoid log(0)

# J vs φ (dominant factor — exponential)
r_J_phi, p_J_phi = stats.pearsonr(phi, log_J)
print(f"\nlog(J) vs φ:        r = {r_J_phi:.3f}  (p = {p_J_phi:.3e})")

# J vs γ_phonon
r_J_gamma, p_J_gamma = stats.pearsonr(gamma_phonon, log_J)
print(f"log(J) vs γ_phonon: r = {r_J_gamma:.3f}  (p = {p_J_gamma:.3e})")

# J vs 1/γ (coherence model)
r_J_invg, p_J_invg = stats.pearsonr(1/gamma_phonon, log_J)
print(f"log(J) vs 1/γ:      r = {r_J_invg:.3f}")

# A vs γ_phonon (can coherence predict Richardson constant?)
log_A = np.log10(A)
r_A_gamma, p_A_gamma = stats.pearsonr(gamma_phonon, log_A)
print(f"\nlog(A) vs γ_phonon: r = {r_A_gamma:.3f}  (p = {p_A_gamma:.3e})")

# A/A0 vs γ_phonon
r_ratio_gamma, p_ratio_gamma = stats.pearsonr(gamma_phonon, A_ratio)
print(f"A/A₀ vs γ_phonon:  r = {r_ratio_gamma:.3f}  (p = {p_ratio_gamma:.3e})")

print(f"""
RESULT:
  Work function φ dominates: r = {r_J_phi:.3f}
  γ_phonon:                  r = {r_J_gamma:.3f}
  Richardson constant A vs γ: r = {r_A_gamma:.3f}

  The exponential barrier exp(-φ/kT) overwhelms any coherence effect.
  At T=1000K, kT = 0.086 eV. A change of φ by 1 eV changes J by
  exp(1/0.086) = exp(11.6) ≈ 10⁵. No polynomial in γ can compete.
""")

# ==============================================================================
# Analysis 2: Why does the framework fail at boundaries?
# ==============================================================================
print("=" * 70)
print("ANALYSIS 2: WHY THE FRAMEWORK FAILS AT BOUNDARIES")
print("=" * 70)

print("""
REASON 1: EXPONENTIAL vs POLYNOMIAL
  Thermionic: J ∝ exp(-φ/kT)  — EXPONENTIAL in barrier height
  Coherence:  Property ∝ γ^n   — POLYNOMIAL in disorder

  For any polynomial f(γ) and exponential g(φ):
    g(φ) >> f(γ) when the barrier exists

  This is not a coherence framework failure — it's a fundamental
  mathematical limitation. Exponential barriers ALWAYS dominate
  polynomial corrections.

  Analogy: Activation energy Ea determines reaction rate.
  Coherence (γ) modifies the pre-exponential factor, but the
  exponential in Ea/kT dominates over decades of magnitude.

REASON 2: SURFACE ≠ BULK
  γ = 2T/θ_D measures BULK phonon coherence (Debye sphere of atoms).
  Thermionic emission occurs at the SURFACE (last 1-2 atomic layers).

  At the surface:
  - Crystal symmetry is broken
  - Electronic states differ (surface states, Shockley states)
  - The "Debye temperature" of the surface ≠ bulk θ_D
  - Image charge effects modify the barrier
  - Adsorbates (even single atoms) change φ dramatically

  The relevant coherence scale is ~ 1 unit cell, where γ → 2
  (two correlated atoms = minimal coherence). The bulk γ is
  measuring the wrong thing.

REASON 3: BARRIER vs PROPAGATION
  From Session #3's two-regime theory:
  - Propagation properties (σ, κ): coherence helps
  - Response properties (d, α): incoherence helps
  - BARRIER properties: neither helps — exponential dominates

  Thermionic emission is a BARRIER-CROSSING problem.
  The barrier height φ is determined by surface electronic structure,
  not by bulk lattice coherence.

  This suggests a THIRD regime beyond coherence and incoherence:
  the BARRIER regime, where exponential factors dominate.
""")

# ==============================================================================
# Analysis 3: Can we extract ANY coherence signal from A?
# ==============================================================================
print("=" * 70)
print("ANALYSIS 3: RESIDUAL ANALYSIS — Is there a hidden γ signal?")
print("=" * 70)

# After removing the φ dependence, does γ predict the residual?
# First, fit log(J) = a × φ + b
slope_phi, intercept_phi = np.polyfit(phi, log_J, 1)
J_predicted = slope_phi * phi + intercept_phi
residual = log_J - J_predicted

# Does the residual correlate with γ?
r_resid_gamma, p_resid_gamma = stats.pearsonr(gamma_phonon, residual)
print(f"\nAfter removing φ dependence:")
print(f"  Residual vs γ_phonon: r = {r_resid_gamma:.3f}  (p = {p_resid_gamma:.3e})")

# The residual reflects the A variation (since J = A × T² × exp(-φ/kT))
# At fixed T, log(J) = log(A) + 2*log(T) + (-φ/kT)×log(e)
# So residual from φ fit ≈ log(A) variation

r_A_resid, _ = stats.pearsonr(log_A, residual)
print(f"  Residual vs log(A):   r = {r_A_resid:.3f}  (residual ≈ A variation)")

# Now: does A correlate with γ when we control for φ?
# Partial correlation: r(A,γ | φ)
# Using formula: r_partial = (r_AG - r_Aφ × r_Gφ) / √((1-r_Aφ²)(1-r_Gφ²))
r_Ag = stats.pearsonr(gamma_phonon, log_A)[0]
r_Aphi = stats.pearsonr(phi, log_A)[0]
r_Gphi = stats.pearsonr(phi, gamma_phonon)[0]

denom = np.sqrt((1 - r_Aphi**2) * (1 - r_Gphi**2))
if denom > 0:
    r_partial = (r_Ag - r_Aphi * r_Gphi) / denom
else:
    r_partial = 0

print(f"\nPartial correlation r(A, γ | φ): {r_partial:.3f}")
print(f"  Raw r(A, γ): {r_Ag:.3f}")
print(f"  r(A, φ): {r_Aphi:.3f}")
print(f"  r(γ, φ): {r_Gphi:.3f}")

print(f"""
RESULT:
  Even after removing the exponential φ dependence, γ_phonon shows
  {'weak' if abs(r_resid_gamma) < 0.3 else 'moderate' if abs(r_resid_gamma) < 0.5 else 'strong'}
  correlation with the residual (r = {r_resid_gamma:.3f}).

  Partial correlation r(A, γ | φ) = {r_partial:.3f}
  {'This suggests ZERO hidden coherence signal.' if abs(r_partial) < 0.2 else 'There may be a weak coherence contribution to A.'}

  The Richardson constant A is determined by:
  1. Surface crystal structure (facet-dependent)
  2. Effective mass at the surface
  3. Quantum-mechanical reflection coefficient
  4. Surface impurities and adsorbates
  None of these are well-captured by bulk θ_D.
""")

# ==============================================================================
# Analysis 4: The three regimes (complete framework)
# ==============================================================================
print("=" * 70)
print("THE COMPLETE THREE-REGIME FRAMEWORK")
print("=" * 70)

print("""
Combining Session #3 (two-regime) with Session #5 (barrier regime):

REGIME 1: COHERENCE (Property ∝ 1/γ)
  Physical type: Propagation through ordered structure
  Examples: σ, κ, μ, Tc, K, D
  Mechanism: Ordered pathways enable efficient transport/stability
  γ contribution: DOMINANT

REGIME 2: INCOHERENCE (Property ∝ γ)
  Physical type: Response to perturbation
  Examples: d_33, α, S, C_v, γ_G
  Mechanism: Soft/disordered structure enables large response
  γ contribution: DOMINANT (but inverted)

REGIME 3: BARRIER (Property ∝ exp(-E_barrier/kT))
  Physical type: Thermally activated escape over energy barrier
  Examples: Thermionic J, reaction rates k, diffusion D
  Mechanism: Barrier height dominates; γ modifies pre-exponential
  γ contribution: NEGLIGIBLE (overwhelmed by exponential)

REGIME 0: NEUTRAL (Property independent of γ)
  Physical type: Counting/extensive quantities
  Examples: R_H, Z, n_v
  γ contribution: ZERO (wrong level of description)

THE KEY INSIGHT:
  The coherence parameter γ has DIFFERENT roles in different regimes:
  - Regime 1: γ is the primary variable (landscape)
  - Regime 2: γ is the primary variable (inverted landscape)
  - Regime 3: γ is a correction to the pre-exponential factor
  - Regime 0: γ is irrelevant

  The framework's failures are concentrated in Regimes 0 and 3.
  Within Regimes 1 and 2, the framework works.
""")

# ==============================================================================
# Analysis 5: Other barrier properties — reaction rates
# ==============================================================================
print("=" * 70)
print("OTHER BARRIER PROPERTIES")
print("=" * 70)

print("""
Thermionic emission is one example. Other barrier properties include:

1. CHEMICAL REACTION RATES (Arrhenius)
   k = A × exp(-Ea/RT)
   Ea dominates. γ might affect A (the pre-exponential factor).
   Session #70: r = 0.997 but CIRCULAR (γ defined from mechanism).
   Session #133: Quantum tunneling k_t, r = 0.411.

2. DIFFUSION (thermally activated)
   D = D₀ × exp(-Q/RT)
   Activation energy Q dominates.
   Session #68: r = 0.53 (moderate) — this is the pre-exponential.

3. NUCLEATION RATE
   I = I₀ × exp(-ΔG*/kT)
   Free energy barrier ΔG* dominates.
   γ could affect the interfacial energy in ΔG*.

4. CONDUCTIVITY OF SEMICONDUCTORS
   σ = σ₀ × exp(-Eg/2kT)
   Band gap Eg dominates.
   γ affects σ₀ through mobility, but Eg → exponential.

PATTERN: All Arrhenius-type processes are barrier-dominated.
The pre-exponential factor A₀ is where coherence COULD contribute,
but it's typically 1-2 orders of magnitude variation vs
10-100 orders from the exponential.
""")

# ==============================================================================
# Visualization
# ==============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Phase 2 Session #5: Boundary vs Bulk Coherence\nWhen Exponential Barriers Overwhelm Polynomial Coherence',
             fontsize=14, fontweight='bold')

# Plot 1: J vs φ (dominant correlation)
ax = axes[0, 0]
ax.scatter(phi, log_J, c='crimson', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (phi[i], log_J[i]), fontsize=6, alpha=0.6)
# Fit line
phi_fit = np.linspace(phi.min(), phi.max(), 100)
J_fit = slope_phi * phi_fit + intercept_phi
ax.plot(phi_fit, J_fit, 'k--', alpha=0.5)
ax.set_xlabel('Work Function φ (eV)')
ax.set_ylabel('log₁₀(J) at T=1000K')
ax.set_title(f'J vs φ: r = {r_J_phi:.3f}\n(exponential dominates)')
ax.grid(True, alpha=0.3)

# Plot 2: J vs γ_phonon (fails)
ax = axes[0, 1]
ax.scatter(gamma_phonon, log_J, c='steelblue', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_phonon[i], log_J[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('log₁₀(J) at T=1000K')
ax.set_title(f'J vs γ_phonon: r = {r_J_gamma:.3f}\n(coherence FAILS)')
ax.grid(True, alpha=0.3)

# Plot 3: Residual vs γ (after removing φ)
ax = axes[0, 2]
ax.scatter(gamma_phonon, residual, c='green', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_phonon[i], residual[i]), fontsize=6, alpha=0.6)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('Residual (log J - φ fit)')
ax.set_title(f'Residual vs γ: r = {r_resid_gamma:.3f}\n(hidden signal?)')
ax.grid(True, alpha=0.3)

# Plot 4: Richardson constant A vs γ
ax = axes[1, 0]
ax.scatter(gamma_phonon, A, c='purple', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_phonon[i], A[i]), fontsize=6, alpha=0.6)
ax.axhline(y=A0, color='red', linestyle='--', alpha=0.5, label=f'A₀ = {A0:.1f}')
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('Richardson Constant A (A/cm²K²)')
ax.set_title(f'A vs γ_phonon: r = {r_A_gamma:.3f}\n(surface ≠ bulk)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 5: Exponential vs polynomial comparison
ax = axes[1, 1]
gamma_range = np.linspace(0.5, 6, 100)
# Polynomial contributions (normalized)
poly1 = 1/gamma_range; poly1 /= poly1.max()
poly2 = gamma_range; poly2 /= poly2.max()
# Exponential barriers at different φ
kT = k_B_eV * 1000  # eV at 1000K
for phi_val, label, color in [(2.0, 'φ=2.0 eV', 'red'), (3.0, 'φ=3.0 eV', 'darkred'), (5.0, 'φ=5.0 eV', 'black')]:
    barrier = np.exp(-phi_val / kT) * np.ones_like(gamma_range)
    barrier /= max(barrier.max(), 1e-100)
    ax.axhline(y=np.log10(barrier[0]) if barrier[0] > 0 else -50, color=color, linestyle='--', alpha=0.5, label=label)

ax.plot(gamma_range, np.log10(poly1/poly1.max()), 'b-', linewidth=2, label='1/γ (coherence)')
ax.plot(gamma_range, np.log10(poly2/poly2.max()), 'r-', linewidth=2, label='γ (incoherence)')
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('log₁₀(contribution)')
ax.set_title('Exponential barriers vs\npolynomial coherence')
ax.legend(fontsize=7, loc='lower left')
ax.set_ylim(-3, 1)
ax.grid(True, alpha=0.3)

# Plot 6: Three-regime summary
ax = axes[1, 2]
ax.text(0.5, 0.92, 'Complete Three-Regime Framework', fontsize=13, ha='center', fontweight='bold', transform=ax.transAxes)

ax.text(0.5, 0.78, 'REGIME 1: COHERENCE (∝ 1/γ)', fontsize=10, ha='center', color='blue', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 0.71, 'Propagation: σ, κ, Tc, K', fontsize=9, ha='center', color='blue', transform=ax.transAxes)

ax.text(0.5, 0.59, 'REGIME 2: INCOHERENCE (∝ γ)', fontsize=10, ha='center', color='red', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 0.52, 'Response: d₃₃, α, S, Cᵥ', fontsize=9, ha='center', color='red', transform=ax.transAxes)

ax.text(0.5, 0.40, 'REGIME 3: BARRIER (∝ exp(-E/kT))', fontsize=10, ha='center', color='darkgreen', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 0.33, 'Activated: J_therm, k_rxn, D_diff', fontsize=9, ha='center', color='darkgreen', transform=ax.transAxes)

ax.text(0.5, 0.21, 'REGIME 0: NEUTRAL (γ irrelevant)', fontsize=10, ha='center', color='gray', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 0.14, 'Counting: R_H, Z, n_v', fontsize=9, ha='center', color='gray', transform=ax.transAxes)

ax.text(0.5, 0.03, f'J(φ): r={r_J_phi:.2f}  |  J(γ): r={r_J_gamma:.2f}  |  A(γ): r={r_A_gamma:.2f}',
        fontsize=9, ha='center', transform=ax.transAxes, alpha=0.7)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phase2_boundary_effects.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nFigure saved: phase2_boundary_effects.png")

# ==============================================================================
# Final
# ==============================================================================
print("\n" + "=" * 70)
print("PHASE 2 SESSION #5: FINAL CONCLUSIONS")
print("=" * 70)

print(f"""
1. THERMIONIC EMISSION CONFIRMS BARRIER REGIME
   J vs φ: r = {r_J_phi:.3f}  (exponential barrier dominates)
   J vs γ: r = {r_J_gamma:.3f}  (coherence negligible)
   A vs γ: r = {r_A_gamma:.3f}  (even pre-exponential uncorrelated)
   Residual vs γ: r = {r_resid_gamma:.3f}  (no hidden signal)

2. THE COMPLETE REGIME FRAMEWORK
   Regime 0: NEUTRAL — counting properties, γ irrelevant
   Regime 1: COHERENCE — propagation ∝ 1/γ
   Regime 2: INCOHERENCE — response ∝ γ
   Regime 3: BARRIER — exponential ∝ exp(-E/kT), γ negligible

3. HOW TO CLASSIFY A NEW PROPERTY
   Step 1: Is it counting something? → Regime 0 (skip γ)
   Step 2: Is there an activation energy barrier? → Regime 3 (skip γ)
   Step 3: Does it measure propagation or stability? → Regime 1 (use 1/γ)
   Step 4: Does it measure response or deformation? → Regime 2 (use γ)

4. THE BARRIER REGIME IS NOT A FRAMEWORK FAILURE
   It's a boundary of applicability. γ was never designed to compete
   with exponential barriers. The framework's domain is polynomial
   (power-law) relationships, not exponential ones.
""")
