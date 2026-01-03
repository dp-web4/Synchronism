#!/usr/bin/env python3
"""
Session #217: Extended Analysis - The Exponent Discrepancy
============================================================

The exponent to exactly match MOND's a₀ is ~1.47, not φ = 1.618.
This 9% difference could indicate:
1. Measurement uncertainty in Ω_m or a₀_MOND
2. A more fundamental relationship we haven't found
3. The golden ratio is approximate, not exact

Author: Autonomous Research Agent
Date: January 3, 2026
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
c = 2.998e8
H0 = 67.4
H0_SI = H0 * 1e3 / (3.086e22)
a_H = c * H0_SI

phi = (1 + np.sqrt(5)) / 2
a0_mond = 1.2e-10

print("=" * 70)
print("Session #217: Extended Analysis - The Exponent Discrepancy")
print("=" * 70)

# =============================================================================
# Part 1: What Exponent Exactly Matches MOND?
# =============================================================================

print("\n" + "=" * 70)
print("Part 1: Solving for the Exact Exponent")
print("=" * 70)

Omega_m_planck = 0.315

# Solve: a₀_MOND = c × H₀ × Ω_m^α
# α = ln(a₀_MOND / (c × H₀)) / ln(Ω_m)
alpha_exact = np.log(a0_mond / a_H) / np.log(Omega_m_planck)

print(f"\nTo match MOND's a₀ = {a0_mond:.2e} m/s²:")
print(f"  Required exponent α = {alpha_exact:.6f}")
print(f"  Golden ratio φ = {phi:.6f}")
print(f"  Difference: {(alpha_exact - phi)/phi * 100:.2f}%")

# What Ω_m would make φ exact?
# a₀_MOND = c × H₀ × Ω_m^φ
# Ω_m = (a₀_MOND / (c × H₀))^(1/φ)
Omega_m_for_phi = (a0_mond / a_H) ** (1/phi)
print(f"\nTo make φ exact, would need:")
print(f"  Ω_m = {Omega_m_for_phi:.4f}")
print(f"  Current Ω_m = {Omega_m_planck:.4f}")
print(f"  Difference: {(Omega_m_for_phi - Omega_m_planck)/Omega_m_planck * 100:.2f}%")

# =============================================================================
# Part 2: Uncertainty Analysis
# =============================================================================

print("\n" + "=" * 70)
print("Part 2: Uncertainty Analysis")
print("=" * 70)

# MOND's a₀ has uncertainty
a0_mond_range = [1.0e-10, 1.2e-10, 1.4e-10]

# Ω_m has uncertainty
Omega_m_range = [0.307, 0.315, 0.323]  # Planck 2018 ± 0.008

# H₀ has uncertainty (and tension!)
H0_range = [67.4, 70.0, 73.0]

print("\nExponent α for different parameter combinations:")
print("-" * 60)

for a0 in a0_mond_range:
    for Om in Omega_m_range:
        alpha = np.log(a0 / a_H) / np.log(Om)
        print(f"  a₀={a0:.1e}, Ω_m={Om:.3f}: α = {alpha:.4f} (vs φ: {(alpha-phi)/phi*100:+.1f}%)")

# =============================================================================
# Part 3: Alternative Golden Ratio Formulas
# =============================================================================

print("\n" + "=" * 70)
print("Part 3: Alternative Formulas with Golden Ratio")
print("=" * 70)

print("\nTrying different formulas involving φ:")

formulas = [
    ("c × H₀ × Ω_m^φ", lambda Om: c * H0_SI * Om**phi),
    ("c × H₀ × Ω_m^(1/φ)", lambda Om: c * H0_SI * Om**(1/phi)),
    ("c × H₀ × Ω_m^(φ-1)", lambda Om: c * H0_SI * Om**(phi-1)),
    ("c × H₀ × Ω_m × φ^(-2)", lambda Om: c * H0_SI * Om * phi**(-2)),
    ("c × H₀ / (2π)", lambda Om: c * H0_SI / (2*np.pi)),
    ("c × H₀ × Ω_m / φ^2", lambda Om: c * H0_SI * Om / phi**2),
    ("c × H₀ × (1-Ω_m)^φ", lambda Om: c * H0_SI * (1-Om)**phi),
    ("c × H₀ × Ω_m^φ × Ω_b/Ω_m", lambda Om: c * H0_SI * Om**phi * 0.049/0.315),
]

print(f"\n{'Formula':<35} | {'Value (m/s²)':<15} | {'Ratio to MOND':<15}")
print("-" * 70)
for name, func in formulas:
    val = func(Omega_m_planck)
    print(f"{name:<35} | {val:.3e} | {val/a0_mond:.4f}")

# =============================================================================
# Part 4: The α = 3/2 Hypothesis
# =============================================================================

print("\n" + "=" * 70)
print("Part 4: Is α = 3/2 the True Exponent?")
print("=" * 70)

alpha_32 = 1.5
a0_32 = c * H0_SI * Omega_m_planck**alpha_32

print(f"""
The exact exponent to match MOND is α = {alpha_exact:.4f}.

What if the true exponent is a simple fraction?

  α = 3/2 = 1.5 gives:
    a₀ = {a0_32:.3e} m/s²
    Ratio to MOND: {a0_32/a0_mond:.4f}

  α = φ = 1.618 gives:
    a₀ = {c * H0_SI * Omega_m_planck**phi:.3e} m/s²
    Ratio to MOND: {c * H0_SI * Omega_m_planck**phi/a0_mond:.4f}

  α = √2 = 1.414 gives:
    a₀ = {c * H0_SI * Omega_m_planck**np.sqrt(2):.3e} m/s²
    Ratio to MOND: {c * H0_SI * Omega_m_planck**np.sqrt(2)/a0_mond:.4f}

OBSERVATION: 3/2 gives a better match than φ!

But why would 3/2 appear?
  - 3D space / 2D surface (holographic)
  - Virial theorem scaling
  - Half-integer spin (fermions)
""")

# =============================================================================
# Part 5: The Ω_m Coincidence
# =============================================================================

print("\n" + "=" * 70)
print("Part 5: What If Ω_m Is Special?")
print("=" * 70)

# Check if Ω_m ≈ 1/π
print(f"\nΩ_m = {Omega_m_planck:.4f}")
print(f"1/π = {1/np.pi:.4f}")
print(f"1/e = {1/np.e:.4f}")
print(f"Ω_m / (1/π) = {Omega_m_planck * np.pi:.4f}")

# If Ω_m = 1/π exactly, what a₀ would we get?
Om_pi = 1/np.pi
a0_if_pi = c * H0_SI * Om_pi**phi
print(f"\nIf Ω_m = 1/π = {Om_pi:.4f}:")
print(f"  a₀ = c × H₀ × (1/π)^φ = {a0_if_pi:.3e} m/s²")
print(f"  Ratio to MOND: {a0_if_pi/a0_mond:.4f}")

# Combined formula check: a₀ = c × H₀ / (2π)
a0_2pi = c * H0_SI / (2 * np.pi)
print(f"\nSimplest formula: a₀ = c × H₀ / (2π)")
print(f"  a₀ = {a0_2pi:.3e} m/s²")
print(f"  Ratio to MOND: {a0_2pi/a0_mond:.4f}")
print(f"  This is {(a0_2pi-a0_mond)/a0_mond*100:.1f}% off from MOND")

# =============================================================================
# Part 6: Visualization
# =============================================================================

print("\n" + "=" * 70)
print("Part 6: Creating Visualization")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Session #217: Exponent Analysis", fontsize=14)

# Panel 1: a₀ vs exponent with key values marked
ax1 = axes[0]
exponents = np.linspace(1.0, 2.0, 100)
a_values = c * H0_SI * Omega_m_planck**exponents

ax1.semilogy(exponents, a_values, 'b-', linewidth=2)
ax1.axhline(y=a0_mond, color='red', linestyle='--', linewidth=2, label=f'MOND a₀')
ax1.axvline(x=phi, color='green', linestyle=':', linewidth=2, label=f'φ = {phi:.3f}')
ax1.axvline(x=alpha_exact, color='orange', linestyle=':', linewidth=2, label=f'α_exact = {alpha_exact:.3f}')
ax1.axvline(x=1.5, color='purple', linestyle=':', linewidth=2, label='3/2')

ax1.set_xlabel('Exponent α in Ω_m^α', fontsize=12)
ax1.set_ylabel('a = c × H₀ × Ω_m^α (m/s²)', fontsize=12)
ax1.set_title('Finding the Right Exponent')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Panel 2: Uncertainty region
ax2 = axes[1]
# Plot a₀ vs Ω_m for different exponents
Om_range = np.linspace(0.28, 0.35, 100)

for alpha, color, label in [(1.5, 'purple', 'α = 3/2'),
                              (phi, 'green', 'α = φ'),
                              (alpha_exact, 'orange', f'α = {alpha_exact:.3f}')]:
    a_vals = c * H0_SI * Om_range**alpha
    ax2.plot(Om_range, a_vals*1e10, color=color, linewidth=2, label=label)

ax2.axhline(y=a0_mond*1e10, color='red', linestyle='--', linewidth=2, label='MOND a₀')
ax2.axvline(x=Omega_m_planck, color='gray', linestyle=':', alpha=0.7)
ax2.axvspan(0.307, 0.323, color='blue', alpha=0.1, label='Ω_m uncertainty')

ax2.set_xlabel('Ω_m', fontsize=12)
ax2.set_ylabel('a₀ (× 10⁻¹⁰ m/s²)', fontsize=12)
ax2.set_title('a₀ vs Ω_m for Different Exponents')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session217_exponent_analysis.png', dpi=150)
plt.close()

print("Saved: session217_exponent_analysis.png")

# =============================================================================
# Part 7: Conclusions
# =============================================================================

print("\n" + "=" * 70)
print("Session #217 Extended: CONCLUSIONS")
print("=" * 70)

print(f"""
KEY FINDINGS ON THE EXPONENT:

1. EXACT MATCH EXPONENT:
   To match MOND's a₀ = 1.2×10⁻¹⁰ m/s² exactly:
   α = {alpha_exact:.4f}

2. COMPARISON OF CANDIDATES:
   α = 3/2 = 1.500: Ratio to MOND = {c * H0_SI * Omega_m_planck**1.5/a0_mond:.3f}
   α = φ = 1.618:   Ratio to MOND = {c * H0_SI * Omega_m_planck**phi/a0_mond:.3f}
   α = √2 = 1.414:  Ratio to MOND = {c * H0_SI * Omega_m_planck**np.sqrt(2)/a0_mond:.3f}

3. THE 3/2 HYPOTHESIS:
   α = 3/2 gives a BETTER match than φ!
   This might indicate:
   - Virial theorem connection (KE/PE = 1/2 → 3/2)
   - Dimensional analysis (3D/2D holographic)
   - Simpler fundamental relationship

4. THE AMBIGUITY:
   Given uncertainties in a₀_MOND and Ω_m, we cannot
   definitively distinguish between:
   - α = φ (golden ratio hypothesis)
   - α = 3/2 (virial theorem hypothesis)
   - α ~ 1.47 (pure empirical fit)

5. RECOMMENDATION:
   Keep both formulas as viable:
   - Synchronism A: a₀ = c × H₀ × Ω_m^φ (golden ratio scaling)
   - Synchronism B: a₀ = c × H₀ × Ω_m^(3/2) (virial scaling)

   Both give a₀ ~ 10⁻¹⁰ m/s², which is the key prediction.
   Future precision on Ω_m and a₀ may distinguish them.

6. THE DEEP MYSTERY:
   Why is a₀ ≈ c × H₀ / (2π) at all?
   This is the robust prediction regardless of exact exponent.
""")

print("=" * 70)
print("Session #217 Extended: COMPLETE")
print("=" * 70)
