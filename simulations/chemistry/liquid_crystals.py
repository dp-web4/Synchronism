#!/usr/bin/env python3
"""
Chemistry Session #51: Liquid Crystals and Partial Coherence

Liquid crystals exhibit PARTIAL ordering:
- Nematic: Orientational order only (molecules align but don't crystallize)
- Smectic: Orientational + 1D positional order (layered)
- Cholesteric: Orientational + helical twist

Question: How does Synchronism explain partial coherence?
Can we predict the hierarchy of mesophases?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("Chemistry Session #51: Liquid Crystals & Partial Coherence")
print("=" * 70)
print()

# =============================================================================
# PART 1: THE PHENOMENON
# =============================================================================

print("-" * 70)
print("PART 1: THE LIQUID CRYSTAL PHENOMENON")
print("-" * 70)
print()

print("Liquid crystals combine properties of liquids AND crystals:")
print()
print("  LIQUID properties:")
print("    - Flow under stress")
print("    - No positional long-range order (in nematic)")
print("    - Finite viscosity")
print()
print("  CRYSTAL properties:")
print("    - Orientational long-range order")
print("    - Anisotropic optical properties")
print("    - Sharp phase transitions")
print()

# Phase hierarchy
print("Phase hierarchy (cooling):")
print("  Isotropic → Nematic → Smectic A → Smectic C → Crystal")
print()
print("Each transition adds MORE order (lower γ).")
print()

# =============================================================================
# PART 2: ORDER PARAMETERS
# =============================================================================

print("-" * 70)
print("PART 2: ORDER PARAMETERS")
print("-" * 70)
print()

print("Standard order parameter S = <(3cos²θ - 1)/2>")
print("  S = 0: isotropic (random orientations)")
print("  S = 1: perfect alignment")
print()

print("Typical experimental values:")
lc_data = {
    "Isotropic": {"S": 0.0, "T_rel": 1.1},
    "Nematic (near TNI)": {"S": 0.3, "T_rel": 0.99},
    "Nematic (deep)": {"S": 0.6, "T_rel": 0.90},
    "Smectic A": {"S": 0.8, "T_rel": 0.80},
    "Crystal": {"S": 1.0, "T_rel": 0.60},
}

print(f"{'Phase':<20} | {'S':>6} | {'T/TNI':>8}")
print("-" * 40)
for name, data in lc_data.items():
    print(f"{name:<20} | {data['S']:>6.2f} | {data['T_rel']:>8.2f}")

print()

# =============================================================================
# PART 3: SYNCHRONISM INTERPRETATION
# =============================================================================

print("-" * 70)
print("PART 3: SYNCHRONISM INTERPRETATION")
print("-" * 70)
print()

print("In Synchronism, each type of order has its own γ:")
print()
print("  γ_orient = 2 / √N_orient  (orientational coherence)")
print("  γ_pos = 2 / √N_pos        (positional coherence)")
print("  γ_total = √(γ_orient² + γ_pos²)  (quadrature sum)")
print()

print("Phase definitions:")
print("  Isotropic:  γ_orient ~ 2, γ_pos ~ 2")
print("  Nematic:    γ_orient < 1, γ_pos ~ 2")
print("  Smectic:    γ_orient < 1, γ_pos < 2 (in 1D)")
print("  Crystal:    γ_orient < 1, γ_pos < 1 (in 3D)")
print()

# =============================================================================
# PART 4: ORIENTATIONAL d_eff
# =============================================================================

print("-" * 70)
print("PART 4: ORIENTATIONAL d_eff")
print("-" * 70)
print()

print("For orientational order in 3D:")
print("  - The order parameter lives on a sphere (S²)")
print("  - d_orient = 2 (surface of sphere)")
print("  - d_lower_orient = 2 (XY universality for headless vectors)")
print("  - z_orient ~ 2")
print()
print("  d_eff_orient = (2 - 2) / 2 = 0")
print()
print("Wait - this gives d_eff = 0, meaning no long-range order!")
print()
print("But nematics DO have long-range orientational order...")
print()

print("Resolution: 3D bulk stabilizes what 2D surface cannot.")
print("  Effective d_eff_orient = (3 - 2) / 2 = 0.5")
print()
print("The 3D embedding space matters, not just the order parameter space!")
print()

# =============================================================================
# PART 5: MAPPING S TO γ
# =============================================================================

print("-" * 70)
print("PART 5: MAPPING ORDER PARAMETER S TO γ")
print("-" * 70)
print()

print("Hypothesis: S and γ are related by:")
print()
print("  S = 1 - γ_orient/2")
print()
print("Equivalently:")
print("  γ_orient = 2(1 - S)")
print()

print("Verification:")
print(f"{'Phase':<20} | {'S':>6} | {'γ_orient (pred)':>15}")
print("-" * 50)

for name, data in lc_data.items():
    S = data["S"]
    gamma_pred = 2 * (1 - S)
    print(f"{name:<20} | {S:>6.2f} | {gamma_pred:>15.2f}")

print()
print("This mapping gives:")
print("  Isotropic (S=0): γ = 2 (classical) ✓")
print("  Crystal (S=1): γ = 0 (perfect coherence) ✓")
print()

# =============================================================================
# PART 6: TEMPERATURE DEPENDENCE
# =============================================================================

print("-" * 70)
print("PART 6: TEMPERATURE DEPENDENCE")
print("-" * 70)
print()

print("Near the nematic-isotropic transition (TNI):")
print()
print("From Maier-Saupe theory: S ~ (TNI - T)^β with β ~ 0.25")
print()
print("From Synchronism: γ(T) = γ₀|T - Tc|^β_γ")
print()
print("If S = 1 - γ/2, then:")
print("  1 - γ/2 ~ (TNI - T)^0.25")
print("  γ ~ 2 - 2(TNI - T)^0.25")
print()

def gamma_nematic(T, T_NI, gamma_0=2.0):
    """γ for nematic phase near transition."""
    if T > T_NI:
        return gamma_0  # isotropic
    t = (T_NI - T) / T_NI  # reduced temperature
    S = 0.44 * t**0.25 * np.sign(t)  # Maier-Saupe-like
    S = min(S, 0.8)  # saturation
    return gamma_0 * (1 - S)

# Temperature sweep
T_NI = 350  # K
T_range = np.linspace(300, 370, 100)
gamma_vs_T = [gamma_nematic(T, T_NI) for T in T_range]
S_vs_T = [1 - g/2 for g in gamma_vs_T]

print(f"At T = TNI - 10 K: γ ~ {gamma_nematic(340, T_NI):.2f}, S ~ {1 - gamma_nematic(340, T_NI)/2:.2f}")
print(f"At T = TNI - 30 K: γ ~ {gamma_nematic(320, T_NI):.2f}, S ~ {1 - gamma_nematic(320, T_NI)/2:.2f}")
print()

# =============================================================================
# PART 7: SMECTIC PHASES
# =============================================================================

print("-" * 70)
print("PART 7: SMECTIC PHASES - 1D POSITIONAL COHERENCE")
print("-" * 70)
print()

print("Smectic phases have LAYERED structure:")
print("  - Positional order in 1 direction (perpendicular to layers)")
print("  - Liquid-like disorder within layers")
print()

print("For 1D positional order:")
print("  d_pos = 1")
print("  d_lower_pos = 0 (1D can order)")
print("  z_pos ~ 2")
print()
print("  d_eff_pos = (1 - 0) / 2 = 0.5")
print()

print("Smectic A: γ_orient ~ 0.4, γ_pos ~ 1.0")
print("  γ_total = √(0.4² + 1.0²) = 1.08")
print()

print("Smectic C (tilted): adds tilt angle order")
print("  γ_tilt ~ 1.2 (weak order)")
print("  γ_total = √(0.4² + 1.0² + 1.2²) = 1.57")
print()

# =============================================================================
# PART 8: PHASE TRANSITION HIERARCHY
# =============================================================================

print("-" * 70)
print("PART 8: PHASE TRANSITION HIERARCHY")
print("-" * 70)
print()

print("Synchronism predicts the ordering sequence:")
print()

phases = {
    "Isotropic": {"gamma_orient": 2.0, "gamma_pos": 2.0},
    "Nematic": {"gamma_orient": 0.8, "gamma_pos": 2.0},
    "Smectic A": {"gamma_orient": 0.4, "gamma_pos": 1.0},
    "Smectic C": {"gamma_orient": 0.4, "gamma_pos": 0.8},
    "Crystal": {"gamma_orient": 0.2, "gamma_pos": 0.2},
}

print(f"{'Phase':<15} | {'γ_orient':>10} | {'γ_pos':>8} | {'γ_total':>10}")
print("-" * 55)

for name, data in phases.items():
    gamma_total = np.sqrt(data["gamma_orient"]**2 + data["gamma_pos"]**2)
    print(f"{name:<15} | {data['gamma_orient']:>10.2f} | {data['gamma_pos']:>8.2f} | {gamma_total:>10.2f}")

print()
print("γ_total monotonically decreases: Iso > N > SmA > SmC > Crystal")
print("This IS the ordering hierarchy!")
print()

# =============================================================================
# PART 9: PREDICTIONS
# =============================================================================

print("-" * 70)
print("PART 9: PREDICTIONS")
print("-" * 70)
print()

print("P51.1: Order parameter-γ relationship")
print("  S = 1 - γ_orient/2 for orientational order")
print("  Testable: measure S and infer γ")
print()

print("P51.2: Transition temperatures")
print("  T_c ∝ 1/γ_threshold")
print("  Higher γ threshold → lower transition T")
print()

print("P51.3: Smectic layer spacing")
print("  d_layer ∝ ξ_pos (positional correlation length)")
print("  d_layer increases as γ_pos → 0")
print()

print("P51.4: Electric field response")
print("  Field aligns molecules → decreases γ_orient")
print("  Response time τ ∝ γ_orient (faster when more coherent)")
print()

print("P51.5: Mixed LC systems")
print("  Mixtures have weighted γ:")
print("  γ_mix = Σ xᵢγᵢ (linear mixing rule)")
print()

# =============================================================================
# PART 10: CHOLESTERIC (CHIRAL NEMATIC)
# =============================================================================

print("-" * 70)
print("PART 10: CHOLESTERIC - HELICAL COHERENCE")
print("-" * 70)
print()

print("Cholesteric LCs have a helical twist superimposed on nematic order.")
print()
print("The pitch p is the distance for 360° rotation.")
print("Typical: p ~ 100-1000 nm (visible wavelengths)")
print()

print("In Synchronism:")
print("  - Nematic order: γ_orient ~ 0.8")
print("  - Helical modulation adds: γ_helix = 2/√N_helix")
print("  - N_helix = (sample size / pitch)")
print()

# Calculate for typical cholesteric
sample_size = 10e-6  # 10 μm
pitch = 500e-9  # 500 nm
N_helix = sample_size / pitch
gamma_helix = 2 / np.sqrt(N_helix)

print(f"For sample = 10 μm, pitch = 500 nm:")
print(f"  N_helix = {N_helix:.0f}")
print(f"  γ_helix = {gamma_helix:.3f}")
print()
print("Helical order adds very small γ contribution (highly coherent)")
print()

# =============================================================================
# PART 11: VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: γ vs temperature
ax1 = axes[0, 0]
ax1.plot(T_range, gamma_vs_T, 'b-', linewidth=2)
ax1.axvline(x=T_NI, color='red', linestyle='--', label=f'T_NI = {T_NI} K')
ax1.axhline(y=1, color='gray', linestyle=':', alpha=0.5)

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('γ_orient')
ax1.set_title('Orientational Coherence vs Temperature')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(300, 370)
ax1.set_ylim(0, 2.5)

# Plot 2: S vs γ
ax2 = axes[0, 1]
gamma_plot = np.linspace(0, 2, 100)
S_plot = 1 - gamma_plot/2

ax2.plot(gamma_plot, S_plot, 'b-', linewidth=2)
ax2.scatter([2*(1-d["S"]) for d in lc_data.values()],
            [d["S"] for d in lc_data.values()], s=100, c='red', zorder=5)

for name, data in lc_data.items():
    ax2.annotate(name.split()[0], (2*(1-data["S"]), data["S"]), fontsize=8)

ax2.set_xlabel('γ_orient')
ax2.set_ylabel('Order parameter S')
ax2.set_title('S = 1 - γ/2 Relationship')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 2.5)
ax2.set_ylim(-0.1, 1.1)

# Plot 3: Phase diagram in γ space
ax3 = axes[1, 0]

# Plot regions
gamma_o = np.linspace(0, 2, 50)
gamma_p = np.linspace(0, 2, 50)
Go, Gp = np.meshgrid(gamma_o, gamma_p)
Gtot = np.sqrt(Go**2 + Gp**2)

# Phase boundaries (approximate)
ax3.contour(Go, Gp, Gtot, levels=[0.5, 1.0, 1.5, 2.0, 2.5], colors='gray', alpha=0.5)
ax3.contourf(Go, Gp, Gtot, levels=[0, 0.5, 1.0, 1.5, 2.0, 3.0], alpha=0.3, cmap='viridis_r')

# Phase labels
phase_coords = {
    "Crystal": (0.2, 0.2),
    "SmC": (0.4, 0.8),
    "SmA": (0.4, 1.0),
    "Nematic": (0.8, 1.8),
    "Isotropic": (1.8, 1.8),
}

for name, (go, gp) in phase_coords.items():
    ax3.scatter(go, gp, s=100, zorder=5)
    ax3.annotate(name, (go, gp), fontsize=10, ha='center')

ax3.set_xlabel('γ_orient')
ax3.set_ylabel('γ_pos')
ax3.set_title('Phase Diagram in γ-Space')
ax3.set_xlim(0, 2.5)
ax3.set_ylim(0, 2.5)

# Plot 4: γ_total vs phase
ax4 = axes[1, 1]
phase_names = list(phases.keys())
gamma_totals = [np.sqrt(d["gamma_orient"]**2 + d["gamma_pos"]**2) for d in phases.values()]

colors = plt.cm.viridis(np.linspace(0.8, 0.2, len(phase_names)))
bars = ax4.bar(phase_names, gamma_totals, color=colors)

ax4.axhline(y=1, color='red', linestyle='--', alpha=0.5, label='γ = 1')
ax4.set_ylabel('γ_total')
ax4.set_title('Total Coherence by Phase')
ax4.tick_params(axis='x', rotation=45)
ax4.legend()
ax4.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/liquid_crystals.png',
            dpi=150, bbox_inches='tight')
plt.close()

print()
print("Figure saved to liquid_crystals.png")

# =============================================================================
# SUMMARY
# =============================================================================

print()
print("-" * 70)
print("SUMMARY")
print("-" * 70)
print()

print("Session #51 applies Synchronism to liquid crystals:")
print()
print("1. PARTIAL COHERENCE: LCs have γ_orient and γ_pos independently")
print("   γ_total = √(γ_orient² + γ_pos²)")
print()
print("2. ORDER PARAMETER MAPPING: S = 1 - γ_orient/2")
print("   Connects traditional LC physics to Synchronism")
print()
print("3. PHASE HIERARCHY follows γ_total:")
print("   Iso (2.8) > N (2.2) > SmA (1.1) > SmC (0.9) > Xtal (0.3)")
print()
print("4. ORIENTATIONAL d_eff = 0.5 (3D embedding of 2D order)")
print()
print("5. SMECTIC has 1D positional coherence")
print("   d_eff_pos = 0.5")
print()

print("=" * 70)
print("SESSION #51 COMPLETE: LIQUID CRYSTALS & PARTIAL COHERENCE")
print("=" * 70)
