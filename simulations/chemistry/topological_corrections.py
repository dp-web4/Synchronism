#!/usr/bin/env python3
"""
Chemistry Session #43: Topological Material Corrections

Session #42 found systematic under-prediction of γ for topological materials.
This session investigates WHY and develops corrections.

Key observation:
- Bi2Se3: pred=0.20, obs=0.60 (error=0.40)
- Cd3As2: pred=0.03, obs=0.40 (error=0.37)

Both are large fractional errors. What's different about topological materials?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

print("=" * 70)
print("Chemistry Session #43: Topological Material Corrections")
print("=" * 70)
print()

# =============================================================================
# PART 1: THE DISCREPANCY
# =============================================================================

print("-" * 70)
print("PART 1: THE DISCREPANCY")
print("-" * 70)
print()

# From Session #42
topological_data = {
    "Bi2Se3 (TI)": {
        "d": 3, "d_lower": 0, "z": 1.5,
        "xi_over_a": 10, "gamma_pred": 0.20, "gamma_obs": 0.60
    },
    "Cd3As2 (Weyl)": {
        "d": 3, "d_lower": 0, "z": 1.0,
        "xi_over_a": 15, "gamma_pred": 0.03, "gamma_obs": 0.40
    },
}

# Non-topological for comparison
conventional_data = {
    "BiFeO3": {"gamma_pred": 1.40, "gamma_obs": 1.50},
    "TbMnO3": {"gamma_pred": 1.52, "gamma_obs": 1.70},
    "YBCO": {"gamma_pred": 1.41, "gamma_obs": 1.10},
    "Fe(Se,Te)": {"gamma_pred": 1.56, "gamma_obs": 1.40},
}

print("Topological materials:")
for name, data in topological_data.items():
    ratio = data["gamma_obs"] / data["gamma_pred"]
    print(f"  {name}: pred={data['gamma_pred']:.2f}, obs={data['gamma_obs']:.2f}, ratio={ratio:.1f}x")

print()
print("Conventional materials:")
for name, data in conventional_data.items():
    ratio = data["gamma_obs"] / data["gamma_pred"]
    print(f"  {name}: pred={data['gamma_pred']:.2f}, obs={data['gamma_obs']:.2f}, ratio={ratio:.1f}x")

print()

# =============================================================================
# PART 2: HYPOTHESIS - SURFACE STATE CONTRIBUTION
# =============================================================================

print("-" * 70)
print("PART 2: HYPOTHESIS - SURFACE STATE CONTRIBUTION")
print("-" * 70)
print()

print("Topological insulators have PROTECTED SURFACE STATES.")
print()
print("These add ADDITIONAL degrees of freedom that don't participate")
print("in bulk coherence, effectively increasing the entropy/fluctuations.")
print()

# Surface contribution model
# γ_total² = γ_bulk² + γ_surface²

print("Model: γ_total² = γ_bulk² + γ_surface²")
print()
print("This is analogous to adding incoherent fluctuations in quadrature.")
print()

# Extract surface contribution
print("Extracted surface contributions:")
for name, data in topological_data.items():
    gamma_bulk = data["gamma_pred"]
    gamma_total = data["gamma_obs"]
    gamma_surface = np.sqrt(gamma_total**2 - gamma_bulk**2)
    print(f"  {name}: γ_surface = √({gamma_total:.2f}² - {gamma_bulk:.2f}²) = {gamma_surface:.2f}")

print()

# =============================================================================
# PART 3: SURFACE STATE PHYSICS
# =============================================================================

print("-" * 70)
print("PART 3: SURFACE STATE PHYSICS")
print("-" * 70)
print()

print("In a topological insulator:")
print("  - Bulk is gapped (insulator)")
print("  - Surface has gapless Dirac cone")
print("  - Surface states span momentum space")
print()

# Surface d_eff calculation
# For 2D surface with linear dispersion (z = 1):
# d_eff_surface = (d_surface - 0) / z = (2 - 0) / 1 = 2

print("Surface state dimensionality:")
print("  d_surface = 2 (2D surface)")
print("  z_surface = 1 (linear Dirac dispersion)")
print("  d_lower_surface = 0 (no lower critical dimension for gapless)")
print()
print("  d_eff_surface = (2 - 0) / 1 = 2")
print()

# Calculate γ_surface from d_eff_surface
def calc_gamma(xi_over_a, d_eff):
    if d_eff == 0:
        return 2.0
    N_corr = (xi_over_a) ** d_eff
    return 2 / np.sqrt(N_corr)

print("Surface coherence depends on surface correlation length ξ_s.")
print("Surface is typically less correlated than bulk (disorder, scattering).")
print()

# =============================================================================
# PART 4: RATIO OF SURFACE TO BULK
# =============================================================================

print("-" * 70)
print("PART 4: RATIO OF SURFACE TO BULK")
print("-" * 70)
print()

# The surface-to-bulk ratio determines the contribution
# R = (surface DOFs) / (bulk DOFs) ∝ L^(d-1) / L^d = 1/L

# For a sample of size L lattice constants:
# N_bulk = L^d
# N_surface = L^(d-1)
# R = N_surface / N_bulk = 1/L

print("Surface-to-bulk ratio:")
print("  R = N_surface / N_bulk = L^(d-1) / L^d = 1/L")
print()
print("For typical samples:")
print("  L ~ 1000 lattice constants → R ~ 0.001")
print()
print("BUT: topological surface states are PROTECTED and contribute")
print("disproportionately to measurable fluctuations!")
print()

# Enhanced surface contribution
# Topological protection means surface fluctuations don't average away
# Effective R_eff >> R_geometric

print("Effective surface ratio (from data):")
for name, data in topological_data.items():
    gamma_bulk = data["gamma_pred"]
    gamma_total = data["gamma_obs"]
    # γ_total² = γ_bulk² + R² × γ_surface²
    # Assume γ_surface ~ 2 (uncorrelated surface)
    gamma_surface_assumed = 2.0
    R_eff = np.sqrt((gamma_total**2 - gamma_bulk**2) / gamma_surface_assumed**2)
    print(f"  {name}: R_eff = {R_eff:.2f}")

print()

# =============================================================================
# PART 5: MODIFIED d_eff FORMULA FOR TOPOLOGICAL MATERIALS
# =============================================================================

print("-" * 70)
print("PART 5: MODIFIED d_eff FORMULA")
print("-" * 70)
print()

print("For topological materials:")
print()
print("  γ_topo = √(γ_bulk² + f_s × γ_surface²)")
print()
print("Where:")
print("  f_s = surface fraction (effective, accounting for protection)")
print("  γ_bulk = 2 / √N_corr_bulk (from standard d_eff)")
print("  γ_surface = 2 / √N_corr_surface")
print()

# Fit f_s from data
def gamma_topo_model(gamma_bulk, f_s, gamma_surface=2.0):
    """Model for topological γ."""
    return np.sqrt(gamma_bulk**2 + f_s * gamma_surface**2)

# Data for fitting
gamma_bulk_arr = np.array([data["gamma_pred"] for data in topological_data.values()])
gamma_obs_arr = np.array([data["gamma_obs"] for data in topological_data.values()])

# Find optimal f_s
def residual(f_s):
    pred = gamma_topo_model(gamma_bulk_arr, f_s)
    return np.sum((pred - gamma_obs_arr)**2)

from scipy.optimize import minimize_scalar
result = minimize_scalar(residual, bounds=(0, 1), method='bounded')
f_s_opt = result.x

print(f"Optimal surface fraction: f_s = {f_s_opt:.3f}")
print()

# Verify fit
print("Verification with f_s correction:")
for name, data in topological_data.items():
    gamma_bulk = data["gamma_pred"]
    gamma_corr = gamma_topo_model(gamma_bulk, f_s_opt)
    gamma_obs = data["gamma_obs"]
    error = abs(gamma_corr - gamma_obs)
    print(f"  {name}: pred_corr={gamma_corr:.2f}, obs={gamma_obs:.2f}, error={error:.2f}")

print()

# =============================================================================
# PART 6: PHYSICAL INTERPRETATION OF f_s
# =============================================================================

print("-" * 70)
print("PART 6: PHYSICAL INTERPRETATION OF f_s")
print("-" * 70)
print()

print(f"f_s = {f_s_opt:.3f} means:")
print()
print("1. Surface contributes √f_s ~ {:.1f}% to fluctuation amplitude".format(np.sqrt(f_s_opt) * 100))
print()
print("2. This is MUCH larger than geometric ratio (~ 0.1%)")
print()
print("3. Topological protection 'amplifies' surface contribution")
print()
print("4. Surface states are ROBUST against disorder")
print("   → don't average to zero")
print("   → contribute coherently to observables")
print()

# =============================================================================
# PART 7: PREDICTIONS FOR OTHER TOPOLOGICAL MATERIALS
# =============================================================================

print("-" * 70)
print("PART 7: PREDICTIONS FOR OTHER TOPOLOGICAL MATERIALS")
print("-" * 70)
print()

# New topological materials to predict
new_topo = {
    "Bi2Te3 (TI)": {"d": 3, "d_lower": 0, "z": 1.5, "xi_over_a": 8},
    "Sb2Te3 (TI)": {"d": 3, "d_lower": 0, "z": 1.5, "xi_over_a": 7},
    "TaAs (Weyl)": {"d": 3, "d_lower": 0, "z": 1.0, "xi_over_a": 12},
    "NbAs (Weyl)": {"d": 3, "d_lower": 0, "z": 1.0, "xi_over_a": 10},
    "MoTe2 (Type-II Weyl)": {"d": 3, "d_lower": 0, "z": 1.2, "xi_over_a": 6},
    "WTe2 (Type-II Weyl)": {"d": 3, "d_lower": 0, "z": 1.2, "xi_over_a": 5},
    "Na3Bi (Dirac)": {"d": 3, "d_lower": 0, "z": 1.0, "xi_over_a": 8},
}

print(f"{'System':<25} | {'d_eff':>5} | {'γ_bulk':>6} | {'γ_corr':>6}")
print("-" * 55)

predictions = []
for name, data in new_topo.items():
    d = data["d"]
    d_lower = data["d_lower"]
    z = data["z"]
    xi_over_a = data["xi_over_a"]

    d_eff = (d - d_lower) / z
    gamma_bulk = calc_gamma(xi_over_a, d_eff)
    gamma_corr = gamma_topo_model(gamma_bulk, f_s_opt)

    print(f"{name:<25} | {d_eff:>5.2f} | {gamma_bulk:>6.2f} | {gamma_corr:>6.2f}")
    predictions.append({"name": name, "gamma_bulk": gamma_bulk, "gamma_corr": gamma_corr})

print()

# =============================================================================
# PART 8: SCALING WITH SAMPLE SIZE
# =============================================================================

print("-" * 70)
print("PART 8: SCALING WITH SAMPLE SIZE")
print("-" * 70)
print()

print("Surface contribution should scale with sample geometry:")
print()
print("  f_s ∝ (surface area) / (volume) ∝ 1/L")
print()
print("For thin films (thickness t):")
print("  f_s ∝ 2/t (top + bottom surfaces)")
print()

# Predict γ vs film thickness
thicknesses = np.array([5, 10, 20, 50, 100, 200, 500])  # nm
lattice_a = 0.5  # nm typical

print(f"Predicted γ vs film thickness (Bi2Se3-like):")
print(f"{'Thickness (nm)':<15} | {'γ_corr':>6}")
print("-" * 25)

for t in thicknesses:
    # f_s scales as 2/L where L is thickness in lattice units
    L = t / lattice_a
    f_s_thick = f_s_opt * (100 / L)  # normalize to ~100 lattice constant reference
    f_s_thick = min(f_s_thick, 0.5)  # cap at reasonable value

    gamma_bulk = 0.20  # Bi2Se3 bulk value
    gamma_thick = gamma_topo_model(gamma_bulk, f_s_thick)
    print(f"{t:<15} | {gamma_thick:>6.2f}")

print()
print("Prediction: Thinner films have LARGER γ (more surface contribution)")
print()

# =============================================================================
# PART 9: COMPARISON TO CONVENTIONAL MATERIALS
# =============================================================================

print("-" * 70)
print("PART 9: WHY CONVENTIONAL MATERIALS DON'T NEED CORRECTION")
print("-" * 70)
print()

print("Conventional materials have:")
print("  1. No protected surface states")
print("  2. Surface disorder averages fluctuations to zero")
print("  3. Bulk dominates by factor of L^d / L^(d-1) ~ 1000")
print()
print("Topological materials have:")
print("  1. Protected surface states (no backscattering)")
print("  2. Surface fluctuations don't average away")
print("  3. Effective surface contribution ~ f_s = {:.3f}".format(f_s_opt))
print()

# =============================================================================
# PART 10: VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Corrected vs uncorrected predictions
ax1 = axes[0]
gamma_bulk_plot = np.linspace(0.01, 0.5, 100)
gamma_uncorr = gamma_bulk_plot
gamma_corr_plot = gamma_topo_model(gamma_bulk_plot, f_s_opt)

ax1.plot(gamma_bulk_plot, gamma_uncorr, 'b--', label='Uncorrected')
ax1.plot(gamma_bulk_plot, gamma_corr_plot, 'r-', label=f'Corrected (f_s={f_s_opt:.2f})', linewidth=2)

# Add data points
ax1.scatter(gamma_bulk_arr, gamma_obs_arr, s=100, c='green', label='Observed', zorder=5)

ax1.set_xlabel('γ_bulk (from d_eff)')
ax1.set_ylabel('γ_total')
ax1.set_title('Surface State Correction')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 0.5)
ax1.set_ylim(0, 1.0)

# Plot 2: γ vs sample thickness
ax2 = axes[1]
L_plot = np.linspace(10, 500, 100)
f_s_plot = f_s_opt * (100 / L_plot)
f_s_plot = np.minimum(f_s_plot, 0.5)
gamma_vs_L = gamma_topo_model(0.20, f_s_plot)

ax2.plot(L_plot * 0.5, gamma_vs_L, 'b-', linewidth=2)  # Convert to nm
ax2.axhline(y=0.20, color='gray', linestyle='--', label='Bulk limit')
ax2.axhline(y=0.60, color='green', linestyle=':', label='Observed (thick)')

ax2.set_xlabel('Film thickness (nm)')
ax2.set_ylabel('γ')
ax2.set_title('Size Dependence')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 250)
ax2.set_ylim(0, 1.0)

# Plot 3: Surface vs bulk contributions
ax3 = axes[2]
f_s_range = np.linspace(0, 0.3, 100)
gamma_bulk_ref = 0.20  # Reference

gamma_total_range = gamma_topo_model(gamma_bulk_ref, f_s_range)
gamma_surface_contrib = np.sqrt(f_s_range) * 2

ax3.plot(f_s_range, gamma_total_range, 'b-', label='γ_total', linewidth=2)
ax3.fill_between(f_s_range, gamma_bulk_ref, gamma_total_range, alpha=0.3, label='Surface contribution')
ax3.axhline(y=gamma_bulk_ref, color='gray', linestyle='--', label='Bulk only')

ax3.axvline(x=f_s_opt, color='red', linestyle=':', label=f'Fitted f_s={f_s_opt:.2f}')

ax3.set_xlabel('Surface fraction f_s')
ax3.set_ylabel('γ')
ax3.set_title('Bulk + Surface Decomposition')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 0.3)
ax3.set_ylim(0, 1.0)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/topological_corrections.png',
            dpi=150, bbox_inches='tight')
plt.close()

print()
print("Figure saved to topological_corrections.png")

# =============================================================================
# SUMMARY
# =============================================================================

print()
print("-" * 70)
print("SUMMARY")
print("-" * 70)
print()

print("Session #43 resolves topological material discrepancy:")
print()
print("1. PROBLEM: Standard d_eff under-predicts γ for TIs and Weyl semimetals")
print()
print("2. CAUSE: Protected surface states add incoherent fluctuations")
print("   Surface states don't average to zero (topological protection)")
print()
print("3. FORMULA: γ_topo = √(γ_bulk² + f_s × γ_surface²)")
print(f"   With f_s = {f_s_opt:.3f} from fitting")
print()
print("4. PREDICTIONS:")
print("   - Thin films have larger γ (more surface)")
print("   - Weyl semimetals: γ ~ 0.4-0.5")
print("   - TIs: γ ~ 0.5-0.6")
print()
print("5. TEST: Measure γ vs film thickness")
print("   Should see γ increase as film gets thinner")
print()

print("=" * 70)
print("SESSION #43 COMPLETE: TOPOLOGICAL CORRECTIONS")
print("=" * 70)
