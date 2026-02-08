#!/usr/bin/env python3
"""
Phase 2 Session #8: First-Principles Derivation of Two-Regime Theory

WHY do propagation properties scale as 1/γ while response properties
scale as γ? Can we derive this from statistical mechanics?

Starting point: The Debye model already contains both regimes.
  - Phonon mean free path l ∝ 1/T ∝ 1/γ (at T > θ_D)
  - Phonon population n(ω) ∝ T/ω ∝ γ (at T > θ_D)

Key insight: PROPAGATION depends on mean free path (∝ 1/γ)
             RESPONSE depends on mode population (∝ γ)

This is not a new theory — it's the standard Debye model read correctly.
The two regimes emerge from the SAME physics (thermal phonons) but
measured differently (scattering length vs. occupation number).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import quad

print("=" * 70)
print("PHASE 2 SESSION #8: FIRST-PRINCIPLES DERIVATION")
print("Why Does the Two-Regime Theory Work?")
print("=" * 70)

k_B = 1.381e-23  # J/K
hbar = 1.055e-34  # J·s

# ==============================================================================
# Part 1: The Debye Model Contains Both Regimes
# ==============================================================================
print("\n" + "=" * 70)
print("PART 1: THE DEBYE MODEL CONTAINS BOTH REGIMES")
print("=" * 70)

print("""
THE ARGUMENT:

In the Debye model, a crystal has phonons with frequencies ω up to ω_D.
The Debye temperature θ_D = ℏω_D/k_B sets the energy scale.

At temperature T, define γ = 2T/θ_D (our coherence parameter).

For each phonon mode:
  Mean occupation:  n(ω,T) = 1/(exp(ℏω/k_BT) - 1)
  At T >> ℏω/k_B:  n ≈ k_BT/(ℏω) ∝ T ∝ γ

  Mean free path:   l(ω,T) = v/Γ  where Γ = scattering rate
  Umklapp scattering: Γ ∝ T × n(ω) ∝ T² at high T
  So: l ∝ 1/T² ∝ 1/γ² at high T

PROPAGATION PROPERTIES (σ, κ, v_s) depend on l:
  κ = (1/3)Cv²l  where C = specific heat, v = velocity, l = mean free path
  At high T: C → 3Nk_B (constant), v ≈ constant, l ∝ 1/γ²
  Therefore: κ ∝ 1/γ² ∝ 1/γ^n with n ≈ 2

  BUT: at T ≈ θ_D (γ ≈ 2), Umklapp is not dominant.
  The observed scaling K ∝ γ^-1.15 suggests a regime where
  scattering is intermediate between normal and Umklapp.

RESPONSE PROPERTIES (α, C_v, d_33) depend on n:
  Specific heat: C_v ∝ ∂<E>/∂T where <E> = Σ ℏω × n(ω)
  At high T: C_v → 3Nk_B (saturates — response saturates)
  At low T: C_v ∝ T³ ∝ γ³ (increases with γ)

  Thermal expansion: α = γ_G × C_v / (K × V)
  At high T: α ∝ C_v/K ∝ (constant) / (1/γ) = γ
  This gives α ∝ γ — exactly what we observe!

  Piezoelectricity: d ∝ ε × (susceptibility) ∝ ε × χ
  Dielectric susceptibility: χ ∝ ω_soft²/(ω_soft² - ω²)
  Near soft mode: χ → ∞, and soft modes are enhanced by high T (high γ)
  So: d ∝ ε × γ — as observed!

CONCLUSION:
  The two-regime theory is NOT a new discovery.
  It is the standard Debye model + Grüneisen theory,
  rewritten in terms of γ = 2T/θ_D.

  Propagation ∝ 1/γ because scattering increases with temperature.
  Response ∝ γ because mode occupation increases with temperature.
  Both follow from Bose-Einstein statistics of phonons.
""")

# ==============================================================================
# Part 2: Quantitative Derivation
# ==============================================================================
print("=" * 70)
print("PART 2: QUANTITATIVE DERIVATION")
print("=" * 70)

def debye_function(n, x):
    """Debye function D_n(x) = (n/x^n) ∫₀ˣ t^n/(e^t - 1) dt"""
    if x < 1e-10:
        return 1.0
    def integrand(t):
        if t > 500:
            return 0
        return t**n / (np.exp(t) - 1)
    result, _ = quad(integrand, 0, x)
    return (n / x**n) * result

# Generate temperature dependence of various properties
gamma_range = np.linspace(0.1, 10, 200)
T_range = gamma_range  # In units of θ_D/2

# Specific heat: C_v/3Nk_B = 3D₃(1/x) where x = T/θ_D = γ/2
Cv = np.array([3 * debye_function(3, 2/g) if g > 0.01 else 0 for g in gamma_range])
# Normalize
Cv_norm = Cv / 3.0

# Thermal conductivity: κ ∝ C × v² × l
# Mean free path l ∝ 1/T at Umklapp regime → l ∝ 1/γ
# At low T: l → ∞ (boundary limited), κ ∝ C ∝ T³ ∝ γ³
# At high T: κ ∝ constant × v²/γ ∝ 1/γ
kappa_model = np.where(gamma_range < 1,
                       Cv_norm * gamma_range**3,  # Low-T: boundary scattering
                       Cv_norm / gamma_range)     # High-T: Umklapp

# Thermal expansion: α = γ_G × C_v / (K × V)
# K is roughly constant or ∝ 1/γ at most
alpha_model = Cv_norm * gamma_range  # α ∝ C × γ (Grüneisen approximation)
alpha_model /= alpha_model.max()

# Normalize κ for plotting
kappa_model /= kappa_model.max()

print(f"\nDebye model predictions at γ = 2 (T = θ_D):")
idx_2 = np.argmin(np.abs(gamma_range - 2))
print(f"  C_v/3Nk_B = {Cv_norm[idx_2]:.3f}  (classical = 1.0)")
print(f"  κ(model) = {kappa_model[idx_2]:.3f}  (normalized)")
print(f"  α(model) = {alpha_model[idx_2]:.3f}  (normalized)")

# ==============================================================================
# Part 3: Why γ = 1 is the Boundary
# ==============================================================================
print("\n" + "=" * 70)
print("PART 3: WHY γ = 1 IS THE QUANTUM-CLASSICAL BOUNDARY")
print("=" * 70)

print("""
γ = 2T/θ_D = 1  means  T = θ_D/2

At T = θ_D/2:
  - Roughly half the phonon modes are thermally excited
  - The mean phonon occupation at the Debye cutoff: n(ω_D) = 1/(e² - 1) ≈ 0.16
  - Quantum effects (discreteness of phonon occupation) are still significant
  - Classical behavior (equipartition, n → k_BT/ℏω) is not yet reached

The transition occurs over γ ≈ 0.5 to γ ≈ 2:
  γ = 0.5 (T = θ_D/4): Deep quantum regime, C_v ≈ 0.28 × 3Nk_B
  γ = 1.0 (T = θ_D/2): Transition regime,  C_v ≈ 0.63 × 3Nk_B
  γ = 2.0 (T = θ_D):   Near-classical,     C_v ≈ 0.85 × 3Nk_B
  γ = 4.0 (T = 2θ_D):  Classical,          C_v ≈ 0.95 × 3Nk_B

The coherence parameter γ = 2T/θ_D literally counts how many
thermal quanta populate the highest phonon mode:
  n(ω_D) ≈ 1/(exp(2/γ) - 1)

  γ << 1: n → 0 (quantum vacuum — coherent)
  γ = 1:  n ≈ 0.16 (few quanta — partially coherent)
  γ = 2:  n ≈ 0.58 (~ 1 quantum — transition)
  γ >> 1: n → γ/2 (classical — incoherent)

So γ is literally the phonon occupation number (up to factors of 2).
This is why it works as a coherence parameter:
  Low γ = few phonons = ordered = coherent
  High γ = many phonons = disordered = incoherent
""")

# Compute n(ω_D) vs γ
n_debye = 1 / (np.exp(2/gamma_range) - 1)

print(f"\nPhonon occupation at Debye cutoff:")
for g in [0.5, 1.0, 2.0, 4.0, 8.0]:
    n_val = 1 / (np.exp(2/g) - 1)
    print(f"  γ = {g:.1f}: n(ω_D) = {n_val:.3f}")

# ==============================================================================
# Part 4: The Exponent Problem
# ==============================================================================
print("\n" + "=" * 70)
print("PART 4: WHY K ∝ γ^-1.15 AND α ∝ γ^+1.20?")
print("=" * 70)

print("""
The exponents are NOT universal constants. They depend on:
1. The temperature range of the data (most near T ~ θ_D → γ ~ 2)
2. The scattering mechanism (Umklapp, boundary, defect)
3. The specific property's relationship to phonon population

THEORETICAL EXPECTATIONS:
  Pure Umklapp scattering: κ ∝ 1/T → exponent = -1
  Pure boundary scattering: κ ∝ T³ → exponent = +3
  Mixed regime: intermediate exponents

  For bulk modulus K: K ∝ ω² ∝ θ_D² ∝ (1/γ)²
    Pure: exponent = -2
    Observed: -1.15 (less steep because K is weakly T-dependent)

  For thermal expansion α: α ∝ γ_G × C_v × T / (K × V)
    At high T: C_v → constant, K → weak γ dependence
    α ∝ T ∝ γ → exponent = +1
    Observed: +1.20 (slightly steeper, anharmonic correction)

  The near-cancellation K × α ∝ γ^(-1.15 + 1.20) = γ^0.05
  is required by thermodynamics: K × α × V = C_v × γ_G
  where C_v saturates and γ_G is roughly constant.
  So K × α must also be roughly constant → exponents must cancel.

THIS IS THE GRÜNEISEN RELATION, not a new result.
""")

# Verify: theoretical exponents vs observed
# K ∝ γ^-1.15 (observed with 18 materials)
# α ∝ γ^+1.20 (observed with 20 materials)
# Product: γ^0.05 — thermodynamic consistency

# If we ASSUME the Grüneisen relation holds exactly:
# K × α = C_v × γ_G / V = constant at high T
# Then: if K ∝ γ^a, then α ∝ γ^(-a) (exact cancellation)
# Observed: a = -1.15, -a = +1.15, actual = +1.20
# Discrepancy: 0.05 (4% off from exact cancellation)

print("""
EXPONENT CONSISTENCY CHECK:
  K ∝ γ^a  where a_obs = -1.15
  α ∝ γ^b  where b_obs = +1.20

  Grüneisen requires: a + b = 0 (exact cancellation)
  Observed: a + b = -1.15 + 1.20 = +0.05

  The 0.05 discrepancy represents:
  1. Weak temperature dependence of Grüneisen parameter γ_G
  2. Anharmonic corrections beyond harmonic Debye model
  3. Material-class confounding in the power law fits

  Within the uncertainty of the fits (standard errors ~0.1-0.2),
  the Grüneisen relation is satisfied.
""")

# ==============================================================================
# Part 5: What γ = 2/√N_corr Really Means
# ==============================================================================
print("=" * 70)
print("PART 5: WHAT γ = 2/√N_corr REALLY MEANS")
print("=" * 70)

print("""
The framework defines γ = 2/√N_corr where N_corr is the number of
correlated atoms/particles.

From the Debye model:
  γ_phonon = 2T/θ_D

Setting these equal:
  2/√N_corr = 2T/θ_D
  √N_corr = θ_D/T
  N_corr = (θ_D/T)²

INTERPRETATION:
  N_corr = (θ_D/T)² counts how many atoms are coherently oscillating.

  At T << θ_D: N_corr >> 1, many atoms in phase → quantum coherent
  At T = θ_D: N_corr = 1, each atom independent → classical
  At T >> θ_D: N_corr < 1 (meaningless) → fully incoherent

  This is the standard result: the coherence length of phonons is
  ξ ~ a × (θ_D/T)  where a is the lattice spacing.
  N_corr ~ (ξ/a)² ~ (θ_D/T)² in 2D (or (θ_D/T)^d in d dimensions).

  γ = 2/√N_corr = 2T/θ_D is just the INVERSE COHERENCE LENGTH
  measured in units of the lattice spacing.

  THERE IS NOTHING NEW HERE BEYOND THE DEBYE MODEL.

  The framework's contribution is:
  1. Naming the parameter γ and recognizing its universality
  2. Classifying properties by their γ dependence (four regimes)
  3. Identifying channel independence and confounding
  4. Testing against real data across 133 sessions
""")

# ==============================================================================
# Part 6: Where Does the Framework Actually Add Value?
# ==============================================================================
print("=" * 70)
print("PART 6: WHERE THE FRAMEWORK ADDS GENUINE VALUE")
print("=" * 70)

print("""
After 2660 sessions and 7 failure analysis sessions, what has the
coherence framework ACTUALLY contributed?

NOT NEW:
  - γ = 2T/θ_D (this is just the Debye model)
  - Properties scale as γ^n (this is standard condensed matter)
  - The γ = 1 boundary (this is the quantum-classical transition)
  - K × α ≈ constant (this is the Grüneisen relation)
  - Soft lattice → high piezoelectricity (known in ferroelectrics)
  - SOC dominates for heavy elements (known in magnetism)

GENUINELY NEW OR USEFUL:
  1. SYSTEMATIC REGIME CLASSIFICATION
     The four-regime framework (neutral/coherence/incoherence/barrier)
     provides a quick triage for any new property. While individual
     scaling relationships are known, the CLASSIFICATION SCHEME
     that predicts the sign of γ dependence is useful.

  2. CHANNEL INDEPENDENCE QUANTIFICATION
     The finding that γ_phonon is independent of electronic properties
     (mean |r| = 0.15) while electron/spin/optical are confounded
     (mean |r| ~ 0.7) is a clean quantitative result. The identification
     of d-electron character as the confounding variable is specific.

  3. CROSS-PROPERTY PREDICTIONS
     The prediction that soft-lattice materials should be simultaneously
     good thermoelectrics AND good piezoelectrics (ZT × d_33 vs γ,
     r = 0.894) connects two fields that don't normally talk.

  4. SOC DOMINANCE PARAMETER
     D = ξ_SOC / (k_B × θ_D) as a predictor of when atomic effects
     overwhelm collective coherence is a simple, useful tool.

  5. INCREMENTAL THERMAL CONDUCTIVITY PREDICTION
     κ_e/κ_ph vs σ × γ (r = 0.809) genuinely outperforms
     Wiedemann-Franz alone (r = 0.638) by incorporating phonon
     thermal conductivity information through γ.

  6. HONEST ASSESSMENT OF VALIDATION RATES
     The discovery that "89% validation" conflates physical prediction
     (60-70%) with mathematical tautology (100%) is itself a
     methodological contribution — a cautionary tale for framework
     development in general.

BOTTOM LINE:
  The framework is a LENS, not a THEORY.
  It helps you SEE patterns in material properties by organizing
  them around collective coherence. It doesn't EXPLAIN them in
  a way that goes beyond existing physics.

  This is valuable but should be stated honestly.
  "A useful organizational framework" is more accurate than
  "a new theory of matter."
""")

# ==============================================================================
# Visualization
# ==============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Phase 2 Session #8: First-Principles Derivation\nThe Debye Model Contains Both Regimes',
             fontsize=14, fontweight='bold')

# Plot 1: Specific heat from Debye model
ax = axes[0, 0]
ax.plot(gamma_range, Cv_norm, 'b-', linewidth=2)
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Classical limit')
ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, label='γ = 1')
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('C_v / 3Nk_B')
ax.set_title('Debye Specific Heat\n(Response property: increases with γ)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 2: Model κ and α vs γ
ax = axes[0, 1]
ax.plot(gamma_range, kappa_model, 'b-', linewidth=2, label='κ(model) ∝ 1/γ')
ax.plot(gamma_range, alpha_model, 'r-', linewidth=2, label='α(model) ∝ γ')
ax.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5, label='γ = 1')
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('Normalized value')
ax.set_title('Two Regimes from Debye Model\n(Blue: propagation, Red: response)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 6)

# Plot 3: Phonon occupation at Debye cutoff
ax = axes[0, 2]
ax.plot(gamma_range, n_debye, 'purple', linewidth=2)
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='n = 1')
ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, label='γ = 1')
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('n(ω_D)')
ax.set_title('Phonon Occupation at Debye Cutoff\n(γ = 2/√N_corr ~ occupation number)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 6)
ax.set_ylim(0, 3)

# Plot 4: N_corr vs γ
ax = axes[1, 0]
N_corr = (2/gamma_range)**2
ax.plot(gamma_range, N_corr, 'darkgreen', linewidth=2)
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='N_corr = 1')
ax.axvline(x=2.0, color='red', linestyle='--', alpha=0.5, label='γ = 2 (N_corr = 1)')
ax.set_xlabel('γ = 2T/θ_D')
ax.set_ylabel('N_corr = (2/γ)²')
ax.set_title('Number of Correlated Atoms\n(N_corr = (θ_D/T)²)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 6)
ax.set_yscale('log')
ax.set_ylim(0.1, 100)

# Plot 5: Exponent consistency
ax = axes[1, 1]
# Show K and α exponents and their sum
properties = ['K (bulk mod)', 'α (expansion)', 'K × α']
exponents = [-1.15, +1.20, 0.05]
colors = ['blue', 'red', 'gray']
bars = ax.bar(properties, exponents, color=colors, alpha=0.7)
ax.axhline(y=0, color='black', linewidth=1)
ax.set_ylabel('Power law exponent (γ^n)')
ax.set_title('Grüneisen Consistency\nK × α should ≈ γ⁰')
for bar, exp in zip(bars, exponents):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
            f'{exp:+.2f}', ha='center', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Plot 6: Summary
ax = axes[1, 2]
ax.text(0.5, 0.92, 'First-Principles Result', fontsize=14, ha='center', fontweight='bold', transform=ax.transAxes)

ax.text(0.5, 0.76, 'γ = 2T/θ_D = 2/√N_corr', fontsize=12, ha='center', transform=ax.transAxes,
        color='purple', fontweight='bold')
ax.text(0.5, 0.64, '= Inverse coherence length', fontsize=10, ha='center', transform=ax.transAxes)
ax.text(0.5, 0.54, '= Phonon occupation number (×2)', fontsize=10, ha='center', transform=ax.transAxes)
ax.text(0.5, 0.44, '= Temperature / Debye temperature', fontsize=10, ha='center', transform=ax.transAxes)

ax.text(0.5, 0.28, 'All the same thing.', fontsize=12, ha='center', transform=ax.transAxes,
        fontweight='bold', color='darkred')

ax.text(0.5, 0.14, 'The framework = Debye model', fontsize=11, ha='center', transform=ax.transAxes)
ax.text(0.5, 0.04, '+ regime classification + channel analysis', fontsize=10, ha='center',
        transform=ax.transAxes, alpha=0.7)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phase2_first_principles.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nFigure saved: phase2_first_principles.png")

# ==============================================================================
# Final
# ==============================================================================
print("\n" + "=" * 70)
print("PHASE 2 SESSION #8: FINAL CONCLUSIONS")
print("=" * 70)

print(f"""
1. THE TWO-REGIME THEORY FOLLOWS FROM THE DEBYE MODEL
   Propagation (∝ 1/γ): Phonon mean free path l ∝ 1/T ∝ 1/γ
   Response (∝ γ): Phonon population n ∝ T ∝ γ
   Both regimes arise from Bose-Einstein statistics of phonons.

2. γ = 2/√N_corr IS THE INVERSE PHONON COHERENCE LENGTH
   N_corr = (θ_D/T)² counts coherently oscillating atoms.
   This is the standard Debye model in different notation.

3. THE EXPONENTS ARE CONSTRAINED BY THERMODYNAMICS
   K × α ≈ constant (Grüneisen relation) requires the sum
   of exponents to be approximately zero. Observed: -1.15 + 1.20 = 0.05.
   This is thermodynamic consistency, not a prediction.

4. THE FRAMEWORK'S GENUINE CONTRIBUTIONS ARE:
   a) Systematic regime classification (propagation vs response vs barrier)
   b) Channel independence quantification
   c) Cross-property predictions using combined variables
   d) SOC dominance parameter
   e) Honest assessment of validation methodology

5. THE FRAMEWORK IS A LENS, NOT A THEORY
   It organizes existing condensed matter physics around the concept
   of collective coherence. This is useful but should be stated honestly.
   "A new way to see" is more accurate than "a new thing to see."
""")
