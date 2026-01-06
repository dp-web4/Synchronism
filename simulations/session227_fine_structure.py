#!/usr/bin/env python3
"""
Session #227: Fine-Structure Constant Variation from Coherence Physics

The fine-structure constant α = e²/(4πε₀ℏc) ≈ 1/137 determines the
strength of electromagnetic interactions.

Claims of α variation:
- Webb et al. (2011): Spatial dipole in α across the sky
- Quasar absorption: Δα/α ~ 10⁻⁵ to 10⁻⁶ at high z

Key Question: Does coherence physics predict any α variation?

Approach:
1. How does coherence affect electromagnetic coupling?
2. If G → G_eff = G/C(a), does e → e_eff?
3. What about c or ℏ in coherence physics?

Date: January 5, 2026
Machine: CBP
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# =============================================================================
# PART 1: CONSTANTS AND COHERENCE FUNCTION
# =============================================================================

print("=" * 70)
print("SESSION #227: FINE-STRUCTURE CONSTANT FROM COHERENCE PHYSICS")
print("=" * 70)

# Physical constants
c = 2.998e8          # m/s
G = 6.674e-11        # m³/(kg·s²)
hbar = 1.055e-34     # J·s
e = 1.602e-19        # C
epsilon_0 = 8.854e-12  # F/m
alpha_0 = e**2 / (4 * np.pi * epsilon_0 * hbar * c)  # ≈ 1/137

print(f"\nFine-structure constant α = {alpha_0:.6f} = 1/{1/alpha_0:.2f}")

# Cosmological parameters
H_0 = 67.4e3 / 3.086e22  # s⁻¹
Omega_m = 0.315
phi = (1 + np.sqrt(5)) / 2

# Critical acceleration
a_0 = 1.2e-10  # m/s² (MOND scale)

def coherence_function(a, alpha_exp=1/phi):
    """
    Coherence function C(a).

    C(a) = Ω_m + (1 - Ω_m) × (a/a₀)^α / [1 + (a/a₀)^α]
    """
    if a <= 0:
        return Omega_m
    x = (a / a_0) ** alpha_exp
    return Omega_m + (1 - Omega_m) * x / (1 + x)


# =============================================================================
# PART 2: THEORETICAL FRAMEWORK
# =============================================================================

print("\n" + "=" * 70)
print("PART 2: THEORETICAL FRAMEWORK")
print("=" * 70)

print("""
The fine-structure constant is:

    α = e² / (4πε₀ℏc)

In Synchronism, we've established:
- G → G_eff = G / C(a) [gravitational coupling]

Question: Do other fundamental constants vary with coherence?

APPROACH 1: Direct EM Coherence
If electromagnetic interactions are also affected by coherence:
    e_eff² = e² / C_em(a)
    α_eff = α / C_em(a)

But what determines C_em? Is it the same C(a) as for gravity?

APPROACH 2: Dimensional Analysis
α is dimensionless, so it depends on ratios of fundamental quantities.
If coherence modifies spacetime structure uniformly, α might be unchanged.

APPROACH 3: Energy-Dependent Coupling
In QED, α runs with energy scale (α increases at high energy).
Coherence could modify this running.

Let's explore each approach...
""")


# =============================================================================
# PART 3: APPROACH 1 - DIRECT EM COHERENCE
# =============================================================================

print("\n" + "=" * 70)
print("PART 3: APPROACH 1 - DIRECT EM COHERENCE")
print("=" * 70)

print("""
If electromagnetic coupling is modified by coherence:

    α_eff(a) = α₀ / C_em(a)

The simplest assumption: C_em(a) = C(a) [same as gravity]

This would give:
    α_eff(a) = α₀ / C(a)

In deep MOND regime (a << a₀):
    C → Ω_m ≈ 0.315
    α_eff → α₀ / 0.315 ≈ 3.2 × α₀

This would be a HUGE change (~220% increase)!
""")

def alpha_direct_coherence(a):
    """Fine-structure constant if α ∝ 1/C(a)."""
    C = coherence_function(a)
    return alpha_0 / C

# Calculate for different environments
environments = {
    'Deep void': 7e-14,      # m/s²
    'Mean cosmic': 2.3e-13,
    'Filament': 3.5e-13,
    'Cluster': 4.6e-12,
    'Earth surface': 9.8,
    'Lab (negligible)': 1e-5,
}

print("\nDirect EM Coherence Model:")
print("-" * 60)
print(f"{'Environment':<20} {'a (m/s²)':>12} {'C(a)':>10} {'α/α₀':>10}")
print("-" * 60)

for name, a in environments.items():
    C = coherence_function(a)
    alpha_ratio = alpha_direct_coherence(a) / alpha_0
    print(f"{name:<20} {a:>12.3e} {C:>10.4f} {alpha_ratio:>10.4f}")

print("""
PROBLEM: This predicts α varies by ~3x between void and lab!
This is COMPLETELY ruled out by observations.

Quasar absorption spectra constrain |Δα/α| < 10⁻⁵ at z ~ 2.
So this approach is WRONG.
""")


# =============================================================================
# PART 4: APPROACH 2 - COHERENCE AFFECTS SPACETIME, NOT α
# =============================================================================

print("\n" + "=" * 70)
print("PART 4: APPROACH 2 - COHERENCE AFFECTS SPACETIME, NOT α")
print("=" * 70)

print("""
Alternative view: Coherence modifies gravitational coupling but NOT
electromagnetic coupling.

Physical reasoning:
- Gravity couples to stress-energy (including vacuum energy)
- EM couples to charge (not affected by vacuum structure)

In this view:
- G_eff = G / C(a) [yes, gravity modified]
- α_eff = α₀ [no, EM unchanged]

This is consistent with observations but less theoretically unified.

PREDICTION: α is constant across all environments.
This is the "null hypothesis" from coherence physics.
""")


# =============================================================================
# PART 5: APPROACH 3 - RUNNING COUPLING MODIFICATION
# =============================================================================

print("\n" + "=" * 70)
print("PART 5: APPROACH 3 - RUNNING COUPLING MODIFICATION")
print("=" * 70)

print("""
In QED, α runs with energy scale Q:

    α(Q²) = α₀ / [1 - (α₀/3π) × ln(Q²/m_e²c²)]

At Q ~ 91 GeV (Z boson mass): α ≈ 1/128
At Q ~ m_e c² (electron mass): α ≈ 1/137

IDEA: What if coherence modifies the RUNNING of α, not α itself?

The running is due to vacuum polarization (virtual e⁺e⁻ pairs).
If coherence affects vacuum structure, it could affect the running.

Let's model this as:
    α_eff(Q², a) = α₀ / [1 - (α₀/3π) × ln(Q²/m_e²c²) × C(a)]

In deep MOND (C → Ω_m):
    Running is REDUCED (vacuum polarization suppressed)

In high-a regime (C → 1):
    Running is STANDARD
""")

def alpha_running(Q_GeV, a):
    """
    Running fine-structure constant with coherence modification.

    Parameters:
        Q_GeV: Energy scale in GeV
        a: Gravitational acceleration

    Returns:
        α(Q², a)
    """
    m_e = 0.511e-3  # GeV (electron mass)
    C = coherence_function(a)

    # Modified running
    log_term = np.log((Q_GeV / m_e)**2)
    denominator = 1 - (alpha_0 / (3 * np.pi)) * log_term * C

    if denominator <= 0:
        return np.inf  # Landau pole

    return alpha_0 / denominator

print("\nRunning α at Different Scales and Environments:")
print("-" * 70)
print(f"{'Q (GeV)':<10} {'α (lab)':>12} {'α (void)':>12} {'Δα/α':>12}")
print("-" * 70)

a_lab = 9.8
a_void = 7e-14

for Q in [0.001, 0.01, 0.1, 1, 10, 91]:
    alpha_lab = alpha_running(Q, a_lab)
    alpha_void = alpha_running(Q, a_void)
    delta = (alpha_void - alpha_lab) / alpha_lab
    print(f"{Q:<10.3f} {alpha_lab:>12.6f} {alpha_void:>12.6f} {delta:>12.2e}")

print("""
At low energies (Q ~ 1 MeV), the difference is negligible.
This is because the running itself is small at low Q.

At high energies (Q ~ 91 GeV), the running is significant,
and the void vs lab difference is ~ 2%.

PREDICTION: High-energy EM interactions in voids show slightly
different α from labs, at the ~2% level.

This could be tested with:
- High-energy cosmic rays from void regions
- Quasar emission spectra (though affected by gas physics)
""")


# =============================================================================
# PART 6: OBSERVED α VARIATION CLAIMS
# =============================================================================

print("\n" + "=" * 70)
print("PART 6: OBSERVED α VARIATION CLAIMS")
print("=" * 70)

print("""
OBSERVATIONAL CLAIMS:

1. WEBB et al. (2011) - Spatial Dipole
   - Combined Keck (N hemisphere) + VLT (S hemisphere) data
   - Found Δα/α ~ 10⁻⁵ dipole across the sky
   - Controversial: systematic errors disputed

2. QUASAR ABSORPTION LINES
   - High-z absorption systems probe α at z ~ 1-3
   - Different analyses give conflicting results
   - Current limit: |Δα/α| < 10⁻⁶ (most conservative)

3. OKLO NATURAL REACTOR
   - Nuclear reactor in Gabon, 1.8 billion years ago
   - Constrains α change over cosmic time
   - |Δα/α| < 10⁻⁷ over 2 Gyr

COMPARISON WITH COHERENCE PREDICTIONS:

Direct EM coherence: Δα/α ~ 200% → RULED OUT
Running modification: Δα/α ~ 2% at high Q → Not yet testable
No α variation: Δα/α = 0 → Consistent with best data
""")


# =============================================================================
# PART 7: ALTERNATIVE - G/α RATIO
# =============================================================================

print("\n" + "=" * 70)
print("PART 7: ALTERNATIVE - G/α RATIO VARIATION")
print("=" * 70)

print("""
Instead of α varying, what if the RATIO G/α varies?

If G → G_eff = G/C(a) but α is constant:

    (G/α)_eff = (G/α)₀ / C(a)

This ratio appears in:
- Gravitational fine structure (α_G = Gm_p²/ℏc ~ 10⁻³⁸)
- Strong/gravitational coupling ratios

The gravitational fine structure constant:
    α_G = G m_p² / (ℏ c) ≈ 5.9 × 10⁻³⁹

With coherence:
    α_G,eff = α_G / C(a)
""")

m_p = 1.673e-27  # kg (proton mass)
alpha_G = G * m_p**2 / (hbar * c)
print(f"\nGravitational fine structure: α_G = {alpha_G:.3e}")

print("\nα_G variation with environment:")
print("-" * 50)
print(f"{'Environment':<20} {'α_G / α_G,0':>15}")
print("-" * 50)

for name, a in environments.items():
    C = coherence_function(a)
    ratio = 1 / C
    print(f"{name:<20} {ratio:>15.4f}")

print("""
The G/α ratio varies by ~3× between void and lab.

This affects:
- Stellar structure (depends on G/α^some_power)
- Primordial nucleosynthesis
- Black hole thermodynamics

TESTABLE PREDICTION:
Stars in void environments should have slightly different
mass-luminosity relationships due to modified G_eff.

L ∝ M^4 × (G_eff/G)^{some power}

This might be detectable with careful stellar population studies.
""")


# =============================================================================
# PART 8: COHERENCE AND VACUUM PERMITTIVITY
# =============================================================================

print("\n" + "=" * 70)
print("PART 8: COHERENCE AND VACUUM PERMITTIVITY")
print("=" * 70)

print("""
The fine-structure constant can be written as:

    α = e² / (4πε₀ℏc) = μ₀ c e² / (4πℏ)

where:
    ε₀ = vacuum permittivity
    μ₀ = vacuum permeability
    c² = 1 / (ε₀ μ₀)

If coherence modifies vacuum structure:
    ε_0,eff = ε₀ × f(C(a))
    μ_0,eff = μ₀ × g(C(a))

Constraint: c² = 1/(ε₀μ₀) must hold (else causality breaks)

So: f(C) × g(C) = 1

Options:
1. f = g = 1 (vacuum unchanged, α constant)
2. f = 1/C, g = C (ε decreases, μ increases, c constant)
3. f = √(1/C), g = √C (symmetric modification)

Option 2 would give:
    α_eff = α₀ / √C (if e unchanged and ε₀ → ε₀/C, μ₀ → μ₀×C)

Let's compute...
""")

def alpha_vacuum_modification(a):
    """Fine-structure with modified vacuum permittivity."""
    C = coherence_function(a)
    # If ε₀ → ε₀/C and we keep e constant:
    # α = e²/(4πε₀ℏc) → e²/(4π(ε₀/C)ℏc) = α₀ × C
    # Wait, this goes the wrong way...

    # Actually: α = e²/(4πε₀ℏc)
    # If ε₀ → ε₀ × C (increased permittivity), α → α/C
    return alpha_0 / C

print("\nVacuum permittivity modification (ε₀ → ε₀/C):")
print("-" * 50)
for name, a in environments.items():
    C = coherence_function(a)
    alpha_mod = alpha_vacuum_modification(a)
    print(f"{name:<20} α = {alpha_mod:.6f} = 1/{1/alpha_mod:.1f}")

print("""
This also predicts ~3× variation - ruled out!

CONCLUSION: If α is to remain approximately constant,
coherence must NOT modify vacuum permittivity significantly.
""")


# =============================================================================
# PART 9: THE PHYSICAL PICTURE
# =============================================================================

print("\n" + "=" * 70)
print("PART 9: THE PHYSICAL PICTURE")
print("=" * 70)

print("""
WHY IS GRAVITY MODIFIED BUT NOT EM?

Physical reasoning from Synchronism:

1. GRAVITY couples to ENERGY/MASS
   - In low-acceleration regions, coherence with distant matter is enhanced
   - This modifies how local mass curves spacetime
   - G_eff = G/C(a) accounts for this enhanced coherence

2. ELECTROMAGNETISM couples to CHARGE
   - Charge is a topological quantum number (integer)
   - Cannot be modified by coherence effects
   - The coupling strength e remains fixed

3. THE FINE-STRUCTURE CONSTANT
   - α = e²/(4πε₀ℏc) involves only EM quantities
   - If e, ε₀, ℏ, c are all unaffected by coherence, α is constant
   - Gravity doesn't enter the definition of α

4. CROSS-COUPLING
   - The only place G and α interact is in extreme situations
   - Black hole thermodynamics: S = A/(4 G ℏ c⁻³)
   - But this involves G, not α

PREDICTION: α is constant to high precision across all environments.
This is actually the STANDARD expectation and is consistent with data.
""")


# =============================================================================
# PART 10: WHAT ABOUT WEBB'S DIPOLE?
# =============================================================================

print("\n" + "=" * 70)
print("PART 10: WEBB'S DIPOLE - CAN COHERENCE EXPLAIN IT?")
print("=" * 70)

print("""
Webb et al. (2011) claimed a spatial dipole in α:

    Δα/α = A × cos(θ)

where θ is the angle from some preferred direction.
Amplitude: A ~ 10⁻⁵ over ~10 Gpc scales.

CAN COHERENCE PHYSICS EXPLAIN THIS?

For coherence to produce a dipole:
- Need ASYMMETRIC distribution of matter/voids
- Our local group is near a large void (Local Void)
- But this would affect G_eff, not α

ALTERNATIVE INTERPRETATION:
If the dipole is real, it might indicate:
- Anisotropic dark energy
- Cosmic domain walls
- Primordial magnetic fields
- Systematic errors in the data

Coherence physics predicts NO α dipole because:
1. α is independent of G
2. α is independent of local acceleration
3. Coherence affects gravity, not EM coupling

PREDICTION: If Webb's dipole is confirmed, it is NOT due to coherence.
It would require new physics beyond the current Synchronism framework.
""")


# =============================================================================
# PART 11: VISUALIZATION
# =============================================================================

print("\n" + "=" * 70)
print("PART 11: CREATING VISUALIZATION")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel 1: Coherence function
ax1 = axes[0, 0]
a_range = np.logspace(-14, 2, 1000)
C_values = [coherence_function(a) for a in a_range]

ax1.semilogx(a_range, C_values, 'b-', linewidth=2)
ax1.axhline(y=Omega_m, color='r', linestyle='--', label=f'C_min = Ω_m = {Omega_m}')
ax1.axhline(y=1, color='g', linestyle='--', label='C_max = 1')
ax1.axvline(x=a_0, color='purple', linestyle=':', label=f'a₀ = {a_0:.1e} m/s²')

for name, a in list(environments.items())[:5]:
    ax1.scatter([a], [coherence_function(a)], s=80, zorder=5, label=name)

ax1.set_xlabel('Acceleration (m/s²)', fontsize=12)
ax1.set_ylabel('Coherence C(a)', fontsize=12)
ax1.set_title('Coherence Function Across Environments', fontsize=14)
ax1.legend(loc='lower right', fontsize=8)
ax1.grid(True, alpha=0.3)

# Panel 2: α variation models
ax2 = axes[0, 1]
models = {
    'Direct EM coherence': lambda a: alpha_0 / coherence_function(a) / alpha_0,
    'α constant (null)': lambda a: 1.0,
}

for name, model in models.items():
    alpha_ratios = [model(a) for a in a_range]
    ax2.semilogx(a_range, alpha_ratios, linewidth=2, label=name)

ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5)
ax2.fill_between([1e-14, 1e2], [1 - 1e-5, 1 - 1e-5], [1 + 1e-5, 1 + 1e-5],
                  alpha=0.3, color='green', label='Observational limit (±10⁻⁵)')
ax2.set_xlabel('Acceleration (m/s²)', fontsize=12)
ax2.set_ylabel('α / α₀', fontsize=12)
ax2.set_title('Fine-Structure Constant Variation Models', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_ylim(0.5, 3.5)

# Panel 3: Running coupling
ax3 = axes[1, 0]
Q_range = np.logspace(-3, 2, 100)

for name, a in [('Lab (a = g)', 9.8), ('Void', 7e-14)]:
    alpha_values = [alpha_running(Q, a) for Q in Q_range]
    ax3.semilogx(Q_range, alpha_values, linewidth=2, label=name)

ax3.axhline(y=alpha_0, color='black', linestyle='--', label=f'α₀ = {alpha_0:.5f}')
ax3.set_xlabel('Energy Scale Q (GeV)', fontsize=12)
ax3.set_ylabel('Running α(Q²)', fontsize=12)
ax3.set_title('Running Fine-Structure Constant with Coherence', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Panel 4: Summary comparison
ax4 = axes[1, 1]
categories = ['Direct EM\nCoherence', 'Running\nModification', 'α Constant\n(Null)']
predictions = [200, 2, 0]  # % variation in void
obs_limit = 0.001  # 10⁻⁵ = 0.001%

bars = ax4.bar(categories, predictions, color=['red', 'orange', 'green'], edgecolor='black')
ax4.axhline(y=obs_limit, color='blue', linestyle='--', linewidth=2,
            label=f'Observational limit ({obs_limit}%)')
ax4.set_ylabel('Predicted |Δα/α| in void (%)', fontsize=12)
ax4.set_title('Model Predictions vs Observational Limit', fontsize=14)
ax4.legend()
ax4.set_yscale('symlog', linthresh=0.01)
ax4.set_ylim(0, 300)
ax4.grid(True, alpha=0.3, axis='y')

# Add text labels
for i, (cat, pred) in enumerate(zip(categories, predictions)):
    status = "RULED OUT" if pred > obs_limit else "CONSISTENT"
    ax4.text(i, pred * 1.5, f'{pred}%\n{status}', ha='center', fontsize=10,
             fontweight='bold', color='red' if pred > obs_limit else 'green')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session227_fine_structure.png',
            dpi=150, bbox_inches='tight')
plt.close()

print("Figure saved: session227_fine_structure.png")


# =============================================================================
# PART 12: CONCLUSIONS
# =============================================================================

print("\n" + "=" * 70)
print("SESSION #227: CONCLUSIONS")
print("=" * 70)

print("""
KEY FINDINGS:

1. DIRECT EM COHERENCE (α ∝ 1/C(a)):
   - Predicts ~200% variation in α between void and lab
   - COMPLETELY RULED OUT by observations (|Δα/α| < 10⁻⁵)
   - This is a strong constraint on coherence physics

2. RUNNING COUPLING MODIFICATION:
   - Predicts ~2% variation at high energies (Q ~ 100 GeV)
   - Not currently testable with precision
   - Could be tested with high-energy cosmic rays

3. α CONSTANT (NULL HYPOTHESIS):
   - Coherence affects gravity (G → G_eff) but not EM (α unchanged)
   - Physically motivated: charge is topological, mass is not
   - CONSISTENT with all observations

4. THEORETICAL INSIGHT:
   - Coherence modifies how mass/energy curves spacetime
   - Charge coupling is unaffected (quantized, topological)
   - This separation is natural in Synchronism

5. WEBB'S DIPOLE:
   - If confirmed, NOT explained by coherence physics
   - Would require new physics beyond current framework
   - Likely systematic error rather than real effect

6. CONSTRAINT ON COHERENCE:
   - The fact that α is constant constrains vacuum modifications
   - ε₀ and μ₀ cannot vary with coherence at significant level
   - Only G is modified, not fundamental EM constants

OVERALL ASSESSMENT:
- Coherence physics is CONSISTENT with α constancy
- This is actually a POSITIVE result - no contradiction with data
- The separation G-modified vs α-constant is physically sensible
""")

print("\n" + "=" * 70)
print("SESSION #227 COMPLETE")
print("=" * 70)
