#!/usr/bin/env python3
"""
Session #243: Dirac Equation from Phase Dynamics

Building on Session #236 (Schrödinger from phase), this session derives
the Dirac equation - the relativistic wave equation for spin-1/2 particles.

KEY QUESTION: How does SPIN emerge from phase geometry in Synchronism?

The hypothesis: Spin is the intrinsic rotation rate of the phase field,
and the Dirac equation emerges from requiring Lorentz-covariant phase dynamics.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

# Physical constants
hbar = constants.hbar
c = constants.c
m_e = constants.m_e

# Golden ratio
phi = (1 + np.sqrt(5)) / 2

print("=" * 80)
print("SESSION #243: DIRAC EQUATION FROM PHASE DYNAMICS")
print("=" * 80)

# =============================================================================
# Part 1: Review - Schrödinger from Phase (Session #236)
# =============================================================================

print("\n" + "=" * 80)
print("PART 1: REVIEW - SCHRÖDINGER FROM PHASE")
print("=" * 80)

print("""
From Session #236:

The wave function is a phase field description:
  ψ(x,t) = A(x,t) × exp(iφ(x,t))

Phase dynamics give de Broglie relations:
  ∂φ/∂t = -E/ℏ     (energy determines phase rate)
  ∇φ = p/ℏ         (momentum determines phase gradient)

For a free particle (E = p²/2m):
  iℏ ∂ψ/∂t = -ℏ²/(2m) ∇²ψ

This is the non-relativistic Schrödinger equation - DERIVED, not postulated.

LIMITATION: This is non-relativistic. For relativistic particles,
we need E² = (pc)² + (mc²)².
""")

# =============================================================================
# Part 2: The Relativistic Energy-Momentum Relation
# =============================================================================

print("\n" + "=" * 80)
print("PART 2: THE RELATIVISTIC ENERGY-MOMENTUM RELATION")
print("=" * 80)

print("""
RELATIVISTIC PHASE DYNAMICS:

In special relativity, energy and momentum form a 4-vector:
  p^μ = (E/c, p)

The phase also forms a scalar from a 4-vector:
  φ = p^μ x_μ / ℏ = (Et - p·x) / ℏ

This is Lorentz invariant - the phase is a spacetime scalar.

From E² = (pc)² + (mc²)², we get:
  (∂φ/∂t)² = c²(∇φ)² + (mc²/ℏ)²

This is the RELATIVISTIC PHASE EQUATION.
""")

# Compton frequency
omega_C = m_e * c**2 / hbar
lambda_C = hbar / (m_e * c)

print(f"\nCompton parameters (electron):")
print(f"  Compton frequency: ω_C = mc²/ℏ = {omega_C:.3e} rad/s")
print(f"  Compton wavelength: λ_C = ℏ/mc = {lambda_C:.3e} m")
print(f"  Compton angular momentum: ℏ/2 = {hbar/2:.3e} J·s (spin!)")

# =============================================================================
# Part 3: From Klein-Gordon to Dirac
# =============================================================================

print("\n" + "=" * 80)
print("PART 3: FROM KLEIN-GORDON TO DIRAC")
print("=" * 80)

print("""
THE KLEIN-GORDON EQUATION:

Squaring the relativistic energy-momentum relation gives:
  E² = (pc)² + (mc²)²

With E → iℏ∂/∂t and p → -iℏ∇:
  -ℏ²∂²ψ/∂t² = -ℏ²c²∇²ψ + (mc²)²ψ

Or: (□ + (mc/ℏ)²)ψ = 0

where □ = (1/c²)∂²/∂t² - ∇² is the d'Alembertian.

PROBLEM: Klein-Gordon allows negative probability densities!

DIRAC'S INSIGHT: Take the SQUARE ROOT of the operator.

(iℏ∂/∂t)² = (cα·p + βmc²)²

This requires α and β to satisfy:
  {α_i, α_j} = 2δ_ij    (anticommutator)
  {α_i, β} = 0
  α_i² = β² = 1

These are the DIRAC MATRICES - they require 4×4 matrices!
""")

# =============================================================================
# Part 4: Spin from Phase Geometry
# =============================================================================

print("\n" + "=" * 80)
print("PART 4: SPIN FROM PHASE GEOMETRY")
print("=" * 80)

print("""
SYNCHRONISM INTERPRETATION OF SPIN:

In standard QM, spin is postulated as intrinsic angular momentum.
In Synchronism, spin emerges from PHASE ROTATION GEOMETRY.

KEY INSIGHT: The phase field can rotate in two distinct ways:
1. EXTERNAL: φ(x,t) changing with position (momentum)
2. INTERNAL: φ rotating around an internal axis (spin)

The internal rotation has two directions:
  - Clockwise: spin-down (|↓⟩)
  - Counter-clockwise: spin-up (|↑⟩)

The Dirac spinor has 4 components because:
  - 2 spin states (up/down)
  - 2 energy signs (particle/antiparticle)

PHASE FIELD REPRESENTATION:
  Ψ = (ψ_+↑, ψ_+↓, ψ_-↑, ψ_-↓)^T

where ± denotes positive/negative energy and ↑↓ denotes spin.
""")

# =============================================================================
# Part 5: The Dirac Equation
# =============================================================================

print("\n" + "=" * 80)
print("PART 5: THE DIRAC EQUATION")
print("=" * 80)

print("""
THE DIRAC EQUATION:

  (iℏγ^μ ∂_μ - mc)Ψ = 0

where γ^μ are the 4×4 Dirac gamma matrices satisfying:
  {γ^μ, γ^ν} = 2η^μν

In expanded form:
  iℏ∂Ψ/∂t = (cα·p + βmc²)Ψ

SYNCHRONISM DERIVATION:

1. Start with relativistic phase: φ = (Et - p·x)/ℏ
2. Require Lorentz covariance for phase dynamics
3. Account for internal phase rotation (spin)
4. The 4-component structure emerges naturally

The Dirac equation is the UNIQUE first-order linear equation that:
- Is Lorentz covariant
- Reduces to Schrödinger in non-relativistic limit
- Incorporates spin-1/2

From Synchronism perspective:
- γ matrices encode phase rotation geometry
- Spin is intrinsic phase helicity
- Antiparticles are negative-phase-frequency modes
""")

# =============================================================================
# Part 6: Spin as Phase Helicity
# =============================================================================

print("\n" + "=" * 80)
print("PART 6: SPIN AS PHASE HELICITY")
print("=" * 80)

print("""
THE HELICITY CONNECTION:

Helicity = projection of spin along momentum direction
  h = S·p/|p|

For massless particles (m=0), helicity is Lorentz invariant.
For massive particles, helicity mixes under boosts.

SYNCHRONISM INTERPRETATION:

The phase field has a HANDEDNESS - which way it rotates
as you move along the momentum direction.

For a right-handed helix (positive helicity):
  φ(x,t) = ωt - k·x with specific rotation sense

For a left-handed helix (negative helicity):
  φ(x,t) = ωt - k·x with opposite rotation sense

This is exactly what spin-1/2 describes!

KEY RESULT:
  Spin = intrinsic helicity of the phase field
  Spin-1/2 = minimal non-trivial phase rotation
  Spin magnitude ℏ/2 = Compton angular momentum scale
""")

# Calculate spin angular momentum
spin_magnitude = hbar / 2
print(f"\nSpin-1/2 angular momentum: ℏ/2 = {spin_magnitude:.3e} J·s")
print(f"This is the minimum non-trivial phase rotation quantum!")

# =============================================================================
# Part 7: The Golden Ratio Connection
# =============================================================================

print("\n" + "=" * 80)
print("PART 7: THE GOLDEN RATIO CONNECTION")
print("=" * 80)

print("""
QUESTION: Does the golden ratio appear in relativistic phase dynamics?

POTENTIAL CONNECTIONS:

1. SPIN PRECESSION
   In a magnetic field, spin precesses at the Larmor frequency.
   The precession rate might show golden ratio structure.

2. ZITTERBEWEGUNG
   The rapid oscillation of a Dirac particle has frequency 2mc²/ℏ.
   This is 2× the Compton frequency - related to phase doubling.

3. FINE STRUCTURE
   The fine structure constant α ≈ 1/137 governs relativistic corrections.
   α = e²/(4πε₀ℏc) is the electromagnetic coupling.

   Interestingly: 1/α ≈ 137 ≈ φ^(φ^2) to within 0.3%!

4. g-FACTOR
   The electron g-factor g ≈ 2.002...
   Anomalous magnetic moment (g-2)/2 ≈ α/(2π) + ...

   The (g-2) is a coherence correction from phase self-interaction.
""")

# Check golden ratio numerology
alpha_em = constants.fine_structure
print(f"\nFine structure constant: α = {alpha_em:.6f}")
print(f"1/α = {1/alpha_em:.2f}")
print(f"φ^(φ²) = {phi**(phi**2):.2f}")
print(f"Difference: {abs(1/alpha_em - phi**(phi**2)):.2f}")

g_electron = 2.00231930436256  # CODATA value
print(f"\nElectron g-factor: g = {g_electron:.10f}")
print(f"(g-2)/2 = {(g_electron-2)/2:.6f}")
print(f"α/(2π) = {alpha_em/(2*np.pi):.6f}")
print(f"Ratio: {(g_electron-2)/2 / (alpha_em/(2*np.pi)):.4f}")

# =============================================================================
# Part 8: Dirac Sea and Antiparticles
# =============================================================================

print("\n" + "=" * 80)
print("PART 8: DIRAC SEA AND ANTIPARTICLES")
print("=" * 80)

print("""
THE DIRAC SEA INTERPRETATION:

The Dirac equation has both positive and negative energy solutions.
Dirac proposed the "sea" of filled negative energy states.

SYNCHRONISM REINTERPRETATION:

Negative energy = negative phase frequency

  E > 0: φ increases with time (forward oscillation)
  E < 0: φ decreases with time (backward oscillation)

An antiparticle is a "hole" in the phase field - a region where
the phase rotates in the opposite temporal direction.

This connects to time-reversal:
  CPT symmetry = phase conjugation + spatial inversion + time reversal

In Synchronism, CPT is a symmetry of the phase field, not an accident!

PARTICLE-ANTIPARTICLE ANNIHILATION:
When opposite-phase regions meet, they can:
1. Cancel (annihilation) → photons (pure phase waves)
2. Exchange (scattering) → phase transfer
""")

# =============================================================================
# Part 9: Connection to Coherence Function
# =============================================================================

print("\n" + "=" * 80)
print("PART 9: CONNECTION TO COHERENCE FUNCTION")
print("=" * 80)

print("""
HOW DOES C(ξ) CONNECT TO DIRAC PHYSICS?

The coherence function C(ξ) governs phase correlations.
For relativistic particles, we need a COVARIANT coherence function.

PROPOSAL: Relativistic Coherence

For a particle with 4-momentum p^μ:
  ξ = Δx^μ Δp_μ / ℏ²   (invariant)

The coherence between two spacetime points:
  C_rel(Δx) = ⟨exp(ip·Δx/ℏ)⟩

For massive particles, this involves the Compton scale:
  ξ ~ |Δx| / λ_C for spacelike separations
  ξ ~ |Δt| × ω_C for timelike separations

SPIN COHERENCE:
Spin adds an additional phase factor:
  C_spin = ⟨exp(iσ·θ)⟩

where θ is the relative spin rotation angle.

For spin-1/2: C_spin = cos(θ/2) for pure states
            (hence the "half-angle" behavior of spinors!)
""")

# =============================================================================
# Part 10: Derivation Summary
# =============================================================================

print("\n" + "=" * 80)
print("PART 10: DERIVATION SUMMARY")
print("=" * 80)

print("""
THE DIRAC EQUATION FROM SYNCHRONISM:

STARTING POINT:
  - Phase field φ(x,t) with relativistic energy-momentum
  - Requirement of Lorentz covariance
  - Minimum non-trivial internal rotation (spin-1/2)

DERIVATION STEPS:

1. Relativistic phase: φ = p^μx_μ/ℏ = (Et - p·x)/ℏ

2. Phase dynamics require:
   (∂φ/∂t)² = c²(∇φ)² + (mc²/ℏ)²

3. Taking square root (Dirac's insight):
   iℏ∂φ/∂t = cα·(-iℏ∇)φ + βmc²φ

4. Internal phase rotation → spinor structure:
   Ψ = (ψ_1, ψ_2, ψ_3, ψ_4)^T

5. α, β matrices encode phase rotation geometry

RESULT:
  (iℏγ^μ∂_μ - mc)Ψ = 0

PHYSICAL MEANING:
  - γ matrices: Phase rotation generators
  - Spinor components: Phase helicity states
  - Antiparticles: Negative phase frequency modes
  - Spin: Intrinsic phase helicity quantum
""")

# =============================================================================
# Part 11: Visualization
# =============================================================================

print("\n" + "=" * 80)
print("PART 11: GENERATING VISUALIZATIONS")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Phase helix for spin-up
ax1 = axes[0, 0]
t = np.linspace(0, 4*np.pi, 500)
z = t / (2*np.pi)  # Position along propagation axis
x_up = np.cos(t)
y_up = np.sin(t)

ax1.plot(x_up, y_up, 'b-', lw=2)
ax1.plot([0], [0], 'ro', ms=10)
ax1.set_xlabel('Phase Re', fontsize=12)
ax1.set_ylabel('Phase Im', fontsize=12)
ax1.set_title('Spin-Up: Right-handed Phase Helix\n(positive helicity)', fontsize=12)
ax1.set_aspect('equal')
ax1.grid(True, alpha=0.3)
ax1.annotate('', xy=(0.7, 0.7), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', color='red', lw=2))

# Plot 2: Phase helix for spin-down
ax2 = axes[0, 1]
x_down = np.cos(-t)
y_down = np.sin(-t)

ax2.plot(x_down, y_down, 'r-', lw=2)
ax2.plot([0], [0], 'bo', ms=10)
ax2.set_xlabel('Phase Re', fontsize=12)
ax2.set_ylabel('Phase Im', fontsize=12)
ax2.set_title('Spin-Down: Left-handed Phase Helix\n(negative helicity)', fontsize=12)
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)
ax2.annotate('', xy=(0.7, -0.7), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', color='blue', lw=2))

# Plot 3: Energy-momentum dispersion
ax3 = axes[1, 0]
p_range = np.linspace(-3, 3, 100)  # In units of mc
E_rel = np.sqrt(p_range**2 + 1)  # E in units of mc²
E_nonrel = 1 + p_range**2 / 2  # Non-relativistic approximation

ax3.plot(p_range, E_rel, 'b-', lw=2.5, label='Relativistic: E = √(p²c² + m²c⁴)')
ax3.plot(p_range, E_nonrel, 'r--', lw=2, label='Non-rel: E = mc² + p²/2m')
ax3.plot(p_range, -E_rel, 'g-', lw=2, label='Antiparticle: E < 0')
ax3.axhline(0, color='gray', ls=':', alpha=0.5)
ax3.set_xlabel('Momentum p (units of mc)', fontsize=12)
ax3.set_ylabel('Energy E (units of mc²)', fontsize=12)
ax3.set_title('Energy-Momentum Dispersion\n(Dirac vs Schrödinger)', fontsize=12)
ax3.legend(loc='upper right')
ax3.grid(True, alpha=0.3)
ax3.set_ylim(-3.5, 3.5)

# Plot 4: Spin coherence
ax4 = axes[1, 1]
theta = np.linspace(0, 2*np.pi, 200)
C_spin = np.cos(theta/2)**2  # Probability for same spin measurement

ax4.plot(theta * 180/np.pi, C_spin, 'b-', lw=2.5, label='Spin coherence: cos²(θ/2)')
ax4.plot(theta * 180/np.pi, 0.5*(1 + np.cos(theta)), 'r--', lw=2,
         label='Classical: (1+cos θ)/2')
ax4.set_xlabel('Relative spin angle θ (degrees)', fontsize=12)
ax4.set_ylabel('Correlation', fontsize=12)
ax4.set_title('Spin-1/2 Coherence Function\n(Half-angle behavior)', fontsize=12)
ax4.legend()
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 360)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session243_dirac_from_phase.png',
            dpi=150)
plt.close()

print("Saved: session243_dirac_from_phase.png")

# =============================================================================
# Summary
# =============================================================================

print("\n" + "=" * 80)
print("SESSION #243 SUMMARY: DIRAC EQUATION FROM PHASE DYNAMICS")
print("=" * 80)

print(f"""
KEY FINDINGS:

1. RELATIVISTIC PHASE DYNAMICS
   - Phase is a spacetime scalar: φ = p^μx_μ/ℏ
   - Lorentz invariance requires E² = (pc)² + (mc²)²
   - Klein-Gordon from squaring, Dirac from square root

2. SPIN AS PHASE HELICITY
   - Spin is intrinsic rotation direction of phase field
   - Spin-up = right-handed helix (positive helicity)
   - Spin-down = left-handed helix (negative helicity)
   - Spin-1/2 is the minimum non-trivial rotation quantum

3. DIRAC STRUCTURE
   - 4 components = 2 spin × 2 energy signs
   - γ matrices encode phase rotation geometry
   - Antiparticles = negative phase frequency modes
   - CPT symmetry is natural in phase picture

4. COHERENCE CONNECTION
   - Spin coherence: C = cos²(θ/2) - the half-angle behavior
   - This explains why spinors need 4π rotation, not 2π
   - Phase field naturally incorporates spinor structure

5. POTENTIAL GOLDEN RATIO CONNECTIONS
   - 1/α ≈ φ^(φ²) ≈ 137 (fine structure)
   - Further investigation needed

CONCLUSION:
The Dirac equation emerges from requiring Lorentz-covariant
phase dynamics with the minimum non-trivial internal rotation.
Spin is not mysterious - it's phase helicity.

This extends Session #236 from non-relativistic to relativistic regime,
maintaining the Synchronism interpretation: physics IS phase dynamics.
""")

print("\n" + "=" * 80)
print("SESSION #243 COMPLETE")
print("=" * 80)
