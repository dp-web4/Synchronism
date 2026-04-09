"""
Phase 4 Session 5: Allotrope Test — Breaking the Bonding-Type Confound

Session #4 found bonding type dominates crystal structure for L_crit (p=0.0002).
The allotrope test uses SAME ELEMENT, DIFFERENT STRUCTURE to deconfound.

Key systems:
- Fe: BCC (alpha) -> FCC (gamma) -> BCC (delta) -> liquid
- Sn: beta (diamond cubic) -> alpha (tetragonal)
- Ti: HCP (alpha) -> BCC (beta)
- Co: HCP (alpha) -> FCC (beta)

If L_crit depends on Voronoi geometry: L should change with structure
If L_crit depends on bonding character: L should be ~constant across allotropes

Also: Cooper pair entity classification in the two-entity framework.

Author: Chemistry Track Phase 4
Date: 2026-04-09
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# PART 1: ALLOTROPE LINDEMANN PARAMETERS FROM LITERATURE
# =============================================================================

print("=" * 70)
print("PART 1: ALLOTROPE TEST — SAME ELEMENT, DIFFERENT STRUCTURE")
print("=" * 70)

# Literature data for Lindemann parameters at melting/transition
# Sources:
#   - Gilvarry (1956): Lindemann parameters for elemental metals
#   - Grimvall et al. (2012): Lattice instabilities in metallic elements
#   - Wallace (1991): Thermodynamics of Crystals, Ch 14
#   - Errandonea (2013): High-pressure melting curves
#   - Lawson et al. (2006): High-T neutron diffraction of Fe

# Note: L = sqrt(<u^2>)/d_nn where d_nn = nearest-neighbor distance
# At transition temperature T_tr, L reaches L_crit for that phase

print("\n--- Iron Allotropes ---")
print("Fe undergoes BCC(alpha) -> FCC(gamma) at 1185K")
print("                FCC(gamma) -> BCC(delta) at 1667K")
print("                BCC(delta) -> liquid     at 1811K")

# Iron lattice parameters and Debye temperatures
# alpha-Fe (BCC): a = 2.87 A, d_nn = a*sqrt(3)/2 = 2.485 A, theta_D = 470K
# gamma-Fe (FCC): a = 3.65 A, d_nn = a/sqrt(2) = 2.581 A, theta_D = 385K
# delta-Fe (BCC): a = 2.93 A, d_nn = a*sqrt(3)/2 = 2.537 A, theta_D = 350K (estimated, softened)

# Mean-square displacements from neutron diffraction (Lawson et al. 2006)
# and Mossbauer spectroscopy (Seto et al. 2003)

# alpha-Fe at 1185K (just before alpha->gamma transition):
# <u^2> from Debye model: <u^2> = (9*hbar^2*T)/(m*k_B*theta_D^2) at high T
# For Fe: m = 55.845 amu = 9.274e-26 kg

hbar = 1.0546e-34  # J*s
k_B = 1.381e-23    # J/K
amu = 1.661e-27     # kg

def u2_debye_highT(T, theta_D, mass_amu):
    """Mean-square displacement in Debye model, high-T limit (T >> theta_D)"""
    m = mass_amu * amu
    return 9 * hbar**2 * T / (m * k_B * theta_D**2)

def lindemann_param(T, theta_D, mass_amu, d_nn_angstrom):
    """Lindemann parameter L = sqrt(<u^2>)/d_nn"""
    u2 = u2_debye_highT(T, theta_D, mass_amu)
    u_rms = np.sqrt(u2)
    d_nn = d_nn_angstrom * 1e-10  # to meters
    return u_rms / d_nn

# Iron phases
Fe_mass = 55.845

# alpha-Fe (BCC)
alpha_Fe = {
    'phase': 'alpha-Fe (BCC)',
    'theta_D': 470,  # K
    'd_nn': 2.485,   # Angstrom
    'T_transition': 1185,  # K (alpha -> gamma)
    'structure': 'BCC'
}

# gamma-Fe (FCC)
gamma_Fe = {
    'phase': 'gamma-Fe (FCC)',
    'theta_D': 385,   # K - lower due to larger volume, magnetic transition
    'd_nn': 2.581,    # Angstrom
    'T_transition': 1667,  # K (gamma -> delta)
    'structure': 'FCC'
}

# delta-Fe (BCC)
delta_Fe = {
    'phase': 'delta-Fe (BCC)',
    'theta_D': 350,   # K - estimated, further softened
    'd_nn': 2.537,    # Angstrom
    'T_transition': 1811,  # K (delta -> liquid, i.e., melting)
    'structure': 'BCC'
}

fe_phases = [alpha_Fe, gamma_Fe, delta_Fe]

print(f"\n{'Phase':<20} {'Structure':<8} {'T_tr (K)':<10} {'theta_D (K)':<12} {'d_nn (A)':<10} {'L at T_tr':<10}")
print("-" * 70)

fe_L_values = {}
for phase in fe_phases:
    L = lindemann_param(phase['T_transition'], phase['theta_D'], Fe_mass, phase['d_nn'])
    fe_L_values[phase['phase']] = L
    print(f"{phase['phase']:<20} {phase['structure']:<8} {phase['T_transition']:<10} {phase['theta_D']:<12} {phase['d_nn']:<10.3f} {L:<10.4f}")

print(f"\n--- Analysis ---")
L_alpha = fe_L_values['alpha-Fe (BCC)']
L_gamma = fe_L_values['gamma-Fe (FCC)']
L_delta = fe_L_values['delta-Fe (BCC)']

print(f"L(alpha-Fe BCC at 1185K) = {L_alpha:.4f}")
print(f"L(gamma-Fe FCC at 1667K) = {L_gamma:.4f}")
print(f"L(delta-Fe BCC at 1811K) = {L_delta:.4f}")

# The Voronoi prediction from Session #4:
# L_crit(BCC) > L_crit(FCC) because BCC has larger inscribed Voronoi radius
# So alpha-Fe and delta-Fe (BCC) should tolerate HIGHER L than gamma-Fe (FCC)

# But these L values are at the TRANSITION temperature, not at melting.
# alpha->gamma is a STRUCTURAL transition, not melting
# gamma->delta is a STRUCTURAL transition, not melting
# delta->liquid IS melting

print(f"\nVoronoi prediction: L_crit(BCC) > L_crit(FCC)")
print(f"  L(BCC, alpha->gamma) = {L_alpha:.4f}")
print(f"  L(FCC, gamma->delta) = {L_gamma:.4f}")
print(f"  L(BCC, delta->liquid) = {L_delta:.4f}")

# CRITICAL DISTINCTION: solid-solid vs solid-liquid transitions
print(f"\n*** CRITICAL: alpha->gamma and gamma->delta are SOLID-SOLID transitions ***")
print(f"*** Only delta->liquid is a true melting transition ***")
print(f"*** Solid-solid transitions occur when FREE ENERGY of new phase is lower ***")
print(f"*** NOT when Lindemann criterion is violated ***")


# =============================================================================
# PART 2: PROPER ALLOTROPE TEST — L AT SAME TEMPERATURE
# =============================================================================

print("\n" + "=" * 70)
print("PART 2: LINDEMANN PARAMETER AT SAME TEMPERATURE (FAIR COMPARISON)")
print("=" * 70)

# Compare L values of different phases AT THE SAME TEMPERATURE
# This isolates the structure effect from the temperature effect

T_compare = 1000  # K — well within alpha-Fe stability range, but gamma-Fe metastable

print(f"\nComparing all Fe phases at T = {T_compare}K:")
print(f"{'Phase':<20} {'Structure':<8} {'theta_D':<12} {'d_nn (A)':<10} {'L(1000K)':<10}")
print("-" * 60)

for phase in fe_phases:
    L = lindemann_param(T_compare, phase['theta_D'], Fe_mass, phase['d_nn'])
    print(f"{phase['phase']:<20} {phase['structure']:<8} {phase['theta_D']:<12} {phase['d_nn']:<10.3f} {L:<10.4f}")

L_alpha_1000 = lindemann_param(T_compare, 470, Fe_mass, 2.485)
L_gamma_1000 = lindemann_param(T_compare, 385, Fe_mass, 2.581)
L_delta_1000 = lindemann_param(T_compare, 350, Fe_mass, 2.537)

print(f"\nAt T=1000K: L(alpha BCC) = {L_alpha_1000:.4f}, L(gamma FCC) = {L_gamma_1000:.4f}")
print(f"Ratio L(FCC)/L(BCC) = {L_gamma_1000/L_alpha_1000:.3f}")

# Voronoi prediction: BCC inscribed radius r_in/d = 0.433, FCC r_in/d = 0.354
# So BCC should TOLERATE more displacement (higher L before escape)
# But here we're comparing L VALUES, not L_crit thresholds

# The question is: does gamma-Fe have LOWER L at the same temperature?
# If L(FCC) < L(BCC) at same T: FCC is stiffer — could explain why BCC transitions first
# If L(FCC) > L(BCC) at same T: the higher L is compensated by geometry

print(f"\n--- What This Means ---")
if L_gamma_1000 > L_alpha_1000:
    print(f"L(FCC) > L(BCC) at same T by factor {L_gamma_1000/L_alpha_1000:.3f}")
    print(f"FCC has LARGER atomic displacements despite smaller Voronoi cell")
    print(f"This is because theta_D(FCC) < theta_D(BCC) — softer lattice")
    print(f"The Voronoi cell size difference is OVERWHELMED by the Debye temperature difference")
else:
    print(f"L(FCC) < L(BCC) at same T — FCC stiffer as expected from larger d_nn")


# =============================================================================
# PART 3: OTHER ALLOTROPIC SYSTEMS
# =============================================================================

print("\n" + "=" * 70)
print("PART 3: OTHER ALLOTROPIC ELEMENTS")
print("=" * 70)

# Titanium: HCP (alpha) -> BCC (beta) at 1155K, melts at 1941K
Ti_mass = 47.867
print("\n--- Titanium ---")
print("Ti: HCP(alpha) -> BCC(beta) at 1155K, melts at 1941K")

ti_alpha = {'phase': 'alpha-Ti (HCP)', 'theta_D': 420, 'd_nn': 2.896, 'T_tr': 1155, 'structure': 'HCP'}
ti_beta  = {'phase': 'beta-Ti (BCC)',  'theta_D': 350, 'd_nn': 2.838, 'T_tr': 1941, 'structure': 'BCC'}

for phase in [ti_alpha, ti_beta]:
    L = lindemann_param(phase['T_tr'], phase['theta_D'], Ti_mass, phase['d_nn'])
    print(f"  {phase['phase']:<20} at T={phase['T_tr']}K: L = {L:.4f}")

# At same temperature
T_comp_Ti = 1000
L_ti_hcp = lindemann_param(T_comp_Ti, 420, Ti_mass, 2.896)
L_ti_bcc = lindemann_param(T_comp_Ti, 350, Ti_mass, 2.838)
print(f"  At T={T_comp_Ti}K: L(HCP) = {L_ti_hcp:.4f}, L(BCC) = {L_ti_bcc:.4f}")
print(f"  Ratio L(BCC)/L(HCP) = {L_ti_bcc/L_ti_hcp:.3f}")

# Cobalt: HCP (alpha) -> FCC (beta) at 695K, melts at 1768K
Co_mass = 58.933
print("\n--- Cobalt ---")
print("Co: HCP(alpha) -> FCC(beta) at 695K, melts at 1768K")

co_alpha = {'phase': 'alpha-Co (HCP)', 'theta_D': 445, 'd_nn': 2.507, 'T_tr': 695, 'structure': 'HCP'}
co_beta  = {'phase': 'beta-Co (FCC)',  'theta_D': 385, 'd_nn': 2.507, 'T_tr': 1768, 'structure': 'FCC'}

for phase in [co_alpha, co_beta]:
    L = lindemann_param(phase['T_tr'], phase['theta_D'], Co_mass, phase['d_nn'])
    print(f"  {phase['phase']:<20} at T={phase['T_tr']}K: L = {L:.4f}")

T_comp_Co = 600
L_co_hcp = lindemann_param(T_comp_Co, 445, Co_mass, 2.507)
L_co_fcc = lindemann_param(T_comp_Co, 385, Co_mass, 2.507)
print(f"  At T={T_comp_Co}K: L(HCP) = {L_co_hcp:.4f}, L(FCC) = {L_co_fcc:.4f}")
print(f"  Ratio L(FCC)/L(HCP) = {L_co_fcc/L_co_hcp:.3f}")

# Tin: beta (tetragonal) -> alpha (diamond cubic) at 286K
# This is the "tin pest" transition — diamond cubic is LESS dense
Sn_mass = 118.71
print("\n--- Tin ---")
print("Sn: beta (tetragonal) stable above 286K, melts at 505K")
print("    alpha (diamond cubic) stable below 286K (tin pest)")

sn_beta  = {'phase': 'beta-Sn (BCT)',   'theta_D': 200, 'd_nn': 3.016, 'T_tr': 505, 'structure': 'BCT'}
sn_alpha = {'phase': 'alpha-Sn (DIA)',   'theta_D': 260, 'd_nn': 2.810, 'T_tr': 286, 'structure': 'DIA'}

# beta-Sn melts at 505K
L_sn_beta_melt = lindemann_param(505, 200, Sn_mass, 3.016)
# alpha-Sn at its transition to beta (286K)
L_sn_alpha_tr = lindemann_param(286, 260, Sn_mass, 2.810)

print(f"  beta-Sn (BCT) at melting 505K:     L = {L_sn_beta_melt:.4f}")
print(f"  alpha-Sn (DIA) at transition 286K: L = {L_sn_alpha_tr:.4f}")

# At same temperature (250K — within alpha stability)
T_comp_Sn = 250
L_sn_beta_250 = lindemann_param(T_comp_Sn, 200, Sn_mass, 3.016)
L_sn_alpha_250 = lindemann_param(T_comp_Sn, 260, Sn_mass, 2.810)
print(f"  At T={T_comp_Sn}K: L(BCT) = {L_sn_beta_250:.4f}, L(DIA) = {L_sn_alpha_250:.4f}")
print(f"  Ratio L(BCT)/L(DIA) = {L_sn_beta_250/L_sn_alpha_250:.3f}")


# =============================================================================
# PART 4: SYNTHESIS — STRUCTURE vs BONDING
# =============================================================================

print("\n" + "=" * 70)
print("PART 4: SYNTHESIS — DOES STRUCTURE OR BONDING CONTROL L?")
print("=" * 70)

# Session #4 Voronoi prediction: L_crit(BCC) > L_crit(FCC) > L_crit(DIA)
# due to inscribed Voronoi radius r_in/d:
#   BCC: r_in/d = 0.433
#   FCC: r_in/d = 0.354
#   HCP: r_in/d = 0.354
#   DIA: r_in/d = 0.217

# At SAME temperature, which phase has higher L?
# The Debye model gives L ∝ sqrt(T) / (theta_D * d_nn)
# So L is controlled by theta_D and d_nn, NOT by Voronoi geometry

print("\n--- The Debye Model Prediction ---")
print("L = sqrt(9*hbar^2*T / (m*k_B*theta_D^2)) / d_nn")
print("At fixed T and m: L ∝ 1/(theta_D * d_nn)")
print()
print("This means L is controlled by:")
print("  1. theta_D (lattice stiffness — BONDING)")
print("  2. d_nn (nearest-neighbor distance — partially STRUCTURE)")
print()
print("The Voronoi escape argument says L_CRIT depends on structure.")
print("The Debye model says L VALUE depends on bonding (theta_D).")
print("These are DIFFERENT questions.")

print("\n--- Separating the Questions ---")
print()
print("Question A: Does L at fixed T depend on crystal structure?")
print("  Answer: YES, through theta_D and d_nn, which change at phase transitions")
print("  But these changes are BONDING-mediated: structure changes -> bonding changes -> theta_D changes")
print()
print("Question B: Does L_CRIT (threshold for structural failure) depend on geometry?")
print("  This requires measuring L at the moment of MELTING for each allotrope")
print("  Problem: only ONE allotrope of each element actually melts (the others")
print("  undergo solid-solid transitions before reaching their L_crit)")

# For Fe, only delta-Fe melts. Alpha-Fe and gamma-Fe undergo structural
# transitions BEFORE they could reach their Lindemann thresholds.
# We can't measure L_crit for alpha-Fe or gamma-Fe because they don't melt.

print("\n--- The Fundamental Problem with the Allotrope Test ---")
print()
print("We CANNOT directly measure L_crit for allotropes that don't melt.")
print("  - alpha-Fe (BCC): transitions to gamma at L = {:.4f} (NOT at melting)".format(L_alpha))
print("  - gamma-Fe (FCC): transitions to delta at L = {:.4f} (NOT at melting)".format(L_gamma))
print("  - delta-Fe (BCC): MELTS at L = {:.4f} (THIS is the true L_crit)".format(L_delta))
print()
print("The solid-solid transition occurs when the FREE ENERGY of the new phase")
print("drops below the current phase — this is thermodynamic, not a Lindemann")
print("instability criterion. The Lindemann criterion only applies to MELTING.")


# =============================================================================
# PART 5: WHAT WE CAN EXTRACT — INDIRECT TEST
# =============================================================================

print("\n" + "=" * 70)
print("PART 5: INDIRECT TEST — LINDEMANN EXTRAPOLATION")
print("=" * 70)

# We CAN ask: at what temperature WOULD each Fe allotrope melt if it were
# stabilized to arbitrarily high temperature?
#
# Using Clausius-Clapeyron / Simon-Glatzel melting curves under pressure:
# At ambient pressure, only delta-Fe melts. But the L_crit should be
# similar (~0.12-0.14 for transition metals) regardless of structure.
#
# Alternative approach: use the Debye model to extrapolate
# T_melt = L_crit^2 * m * k_B * theta_D^2 * d_nn^2 / (9 * hbar^2)

def T_melt_from_Lcrit(L_crit, theta_D, mass_amu, d_nn_angstrom):
    """Predicted melting temperature from Lindemann criterion"""
    m = mass_amu * amu
    d = d_nn_angstrom * 1e-10
    return L_crit**2 * m * k_B * theta_D**2 * d**2 / (9 * hbar**2)

# Use delta-Fe's L as the "true" L_crit for Fe
L_crit_Fe = L_delta
print(f"delta-Fe melting: L_crit = {L_crit_Fe:.4f}")
print()

# Hypothetical melting temperatures if each phase were stable to melting
T_melt_alpha = T_melt_from_Lcrit(L_crit_Fe, 470, Fe_mass, 2.485)
T_melt_gamma = T_melt_from_Lcrit(L_crit_Fe, 385, Fe_mass, 2.581)
T_melt_delta_check = T_melt_from_Lcrit(L_crit_Fe, 350, Fe_mass, 2.537)

print(f"Hypothetical T_melt if L_crit is UNIVERSAL for Fe (structure-independent):")
print(f"  alpha-Fe (BCC, theta_D=470K): T_melt = {T_melt_alpha:.0f}K")
print(f"  gamma-Fe (FCC, theta_D=385K): T_melt = {T_melt_gamma:.0f}K")
print(f"  delta-Fe (BCC, theta_D=350K): T_melt = {T_melt_delta_check:.0f}K (actual: 1811K)")

print(f"\n--- Now assume L_crit IS structure-dependent (Voronoi) ---")
# Voronoi prediction: L_crit(BCC)/L_crit(FCC) = r_in(BCC)/r_in(FCC) = 0.433/0.354 = 1.223
# So L_crit(BCC) ~ 1.22 * L_crit(FCC)
voronoi_ratio = 0.433 / 0.354
L_crit_FCC_voronoi = L_crit_Fe / voronoi_ratio  # Assuming delta-Fe gives BCC L_crit
print(f"Voronoi ratio BCC/FCC = {voronoi_ratio:.3f}")
print(f"If L_crit(BCC) = {L_crit_Fe:.4f}, then L_crit(FCC) = {L_crit_FCC_voronoi:.4f}")

T_melt_gamma_voronoi = T_melt_from_Lcrit(L_crit_FCC_voronoi, 385, Fe_mass, 2.581)
print(f"gamma-Fe T_melt with Voronoi L_crit: {T_melt_gamma_voronoi:.0f}K")
print(f"gamma-Fe T_melt with universal L_crit: {T_melt_gamma:.0f}K")


# =============================================================================
# PART 6: THE KEY INSIGHT — theta_D CONFOUND REDUX
# =============================================================================

print("\n" + "=" * 70)
print("PART 6: THE KEY INSIGHT — theta_D DOMINATES EVERYTHING")
print("=" * 70)

print("""
The Lindemann parameter in the Debye model is:
  L = sqrt(9*hbar^2*T / (m*k_B)) / (theta_D * d_nn)

At the melting point T_m, L = L_crit, so:
  T_m = L_crit^2 * (m*k_B/9*hbar^2) * theta_D^2 * d_nn^2

The Lindemann FORMULA predicts T_m from theta_D and d_nn.
This is known — it's the Lindemann LAW (1910).

The question "does L_crit depend on structure?" is asking:
  "Is there a RESIDUAL structure dependence after theta_D and d_nn?"

For Fe allotropes:
  - theta_D changes from 470 (BCC) to 385 (FCC) to 350 (BCC)
  - d_nn changes from 2.485 to 2.581 to 2.537

theta_D changes by 34%. d_nn changes by 4%.
theta_D dominates by (34/4)^2 ~ 70x in variance.

The "structure effect" on L IS the theta_D effect.
And theta_D is determined by bonding + structure inseparably.
""")

# Quantify: how much of L variation is theta_D vs d_nn?
# L ∝ 1/(theta_D * d_nn) at fixed T
# Variance decomposition

thetas = np.array([470, 385, 350])
d_nns = np.array([2.485, 2.581, 2.537])

# Relative variations
theta_cv = np.std(thetas) / np.mean(thetas)
d_nn_cv = np.std(d_nns) / np.mean(d_nns)

print(f"Coefficient of variation:")
print(f"  theta_D: {theta_cv:.3f} ({theta_cv*100:.1f}%)")
print(f"  d_nn:    {d_nn_cv:.3f} ({d_nn_cv*100:.1f}%)")
print(f"  Ratio of variance contributions: {(theta_cv/d_nn_cv)**2:.1f}x")
print(f"\ntheta_D explains ~{(theta_cv**2/(theta_cv**2+d_nn_cv**2))*100:.0f}% of L variation across Fe allotropes")


# =============================================================================
# PART 7: CROSS-ELEMENT ALLOTROPE COMPARISON
# =============================================================================

print("\n" + "=" * 70)
print("PART 7: CROSS-ELEMENT ALLOTROPE COMPARISON TABLE")
print("=" * 70)

# Collect all allotrope data for comparison
allotrope_data = [
    # (Element, Phase, Structure, theta_D, d_nn, T_transition, transition_type)
    ('Fe', 'alpha', 'BCC', 470, 2.485, 1185, 'solid-solid'),
    ('Fe', 'gamma', 'FCC', 385, 2.581, 1667, 'solid-solid'),
    ('Fe', 'delta', 'BCC', 350, 2.537, 1811, 'melting'),
    ('Ti', 'alpha', 'HCP', 420, 2.896, 1155, 'solid-solid'),
    ('Ti', 'beta',  'BCC', 350, 2.838, 1941, 'melting'),
    ('Co', 'alpha', 'HCP', 445, 2.507, 695,  'solid-solid'),
    ('Co', 'beta',  'FCC', 385, 2.507, 1768, 'melting'),
    ('Sn', 'alpha', 'DIA', 260, 2.810, 286,  'solid-solid'),
    ('Sn', 'beta',  'BCT', 200, 3.016, 505,  'melting'),
]

masses = {'Fe': 55.845, 'Ti': 47.867, 'Co': 58.933, 'Sn': 118.71}

print(f"\n{'Elem':<5} {'Phase':<8} {'Struct':<6} {'theta_D':<8} {'d_nn':<6} {'T_tr':<7} {'L(T_tr)':<8} {'Type':<12}")
print("-" * 62)

for elem, phase, struct, theta, d, T_tr, tr_type in allotrope_data:
    L = lindemann_param(T_tr, theta, masses[elem], d)
    marker = " <-- L_crit" if tr_type == 'melting' else ""
    print(f"{elem:<5} {phase:<8} {struct:<6} {theta:<8} {d:<6.3f} {T_tr:<7} {L:<8.4f} {tr_type}{marker}")

print(f"\n--- Pattern ---")
print("For every element: the L at solid-solid transition < L at melting")
print("This is expected: solid-solid transitions occur at LOWER T than melting")
print("The solid-solid L values are NOT L_crit — they're just L at that temperature")


# =============================================================================
# PART 8: WHAT THE ALLOTROPE TEST ACTUALLY TELLS US
# =============================================================================

print("\n" + "=" * 70)
print("PART 8: VERDICT — WHAT THE ALLOTROPE TEST REVEALS")
print("=" * 70)

print("""
SESSION #4 QUESTION: Does L_crit depend on crystal structure or bonding type?

THE ALLOTROPE TEST CANNOT ANSWER THIS DIRECTLY because:

1. Only ONE allotrope of each element melts (the high-T phase)
2. Other allotropes undergo solid-solid transitions before reaching L_crit
3. Solid-solid transitions are FREE ENERGY crossings, not Lindemann instabilities
4. We cannot measure L_crit for a phase that doesn't melt

WHAT WE CAN EXTRACT:

A. theta_D is the dominant variable controlling L (explains ~98% of variance
   across Fe allotropes), not geometric packing

B. theta_D itself changes across allotropes:
   - Fe: BCC(alpha) 470K -> FCC(gamma) 385K -> BCC(delta) 350K
   - This is a BONDING effect (electron density at crystal sites changes)
   - BCC Fe has stronger directional bonding than FCC Fe (magnetic + d-band)

C. The original Session #4 finding STANDS reinforced:
   Bonding type (controlling theta_D) > crystal structure (controlling geometry)
   This is NOT resolvable by allotrope tests because theta_D and structure
   change TOGETHER at phase transitions — they are thermodynamically coupled

D. The Voronoi prediction (BCC tolerates higher L) is UNTESTABLE in practice:
   - Would need to melt BCC and FCC of the same element at the same theta_D
   - This requires a hypothetical material that doesn't exist
   - Or: measuring L_crit under extreme pressure where phase boundaries shift

CONCLUSION: The allotrope test is a PRODUCTIVE DEAD END.
It clarifies WHY the bonding-structure confound cannot be broken by
elemental allotropes: theta_D and structure are thermodynamically inseparable.
The confound is not an experimental limitation — it's a physical coupling.
""")


# =============================================================================
# PART 9: COOPER PAIRS IN THE TWO-ENTITY FRAMEWORK
# =============================================================================

print("=" * 70)
print("PART 9: COOPER PAIRS — OSCILLATORY OR STRUCTURAL ENTITY?")
print("=" * 70)

print("""
Phase 4 Session #3 established two entity types:
  OSCILLATORY: gamma/f < 1 (wave survives one cycle)
  STRUCTURAL:  L < L_crit  (positional order persists)

Where do Cooper pairs fit?

--- Cooper Pair Properties ---

1. PHASE COHERENCE: The superconducting condensate has a well-defined
   macroscopic phase phi. Phase fluctuations delta_phi are quantum-limited.

   For a BCS superconductor:
     delta_phi ~ sqrt(k_B*T / E_cond)
   where E_cond = N(E_F) * Delta^2 / 2 is the condensation energy.

   At T << Tc: delta_phi -> 0 (STRONGLY oscillatory entity)
   At T -> Tc: delta_phi -> pi (phase coherence lost)

2. SPATIAL EXTENT: Cooper pair wavefunction extends over coherence length xi_0.
   - BCS: xi_0 = hbar*v_F / (pi*Delta) ~ 10-1000 nm
   - This is NOT positional order (pairs are delocalized)
   - It IS correlation length (how far the phase relationship extends)

3. ORDER PARAMETER: Delta(r) = |Delta| * exp(i*phi)
   - Amplitude |Delta|: gap size (energy scale)
   - Phase phi: OSCILLATORY entity character (Josephson frequency = 2eV/hbar)

--- Entity Classification ---

Cooper pairs are PURELY OSCILLATORY entities:
  - gamma_pair ~ 1/tau_GL ~ (T - Tc)/Tc near Tc
  - f_pair = 2*Delta/hbar (Josephson frequency for a pair)
  - gamma/f = hbar / (2*Delta*tau_GL)

  At T << Tc: tau_GL -> infinity, so gamma/f -> 0 (perfect oscillatory entity)
  At T -> Tc: Delta -> 0, tau_GL -> 0, gamma/f -> infinity (entity destroyed)
""")

# Quantitative: gamma/f for some superconductors
print("--- Quantitative gamma/f for Cooper Pairs ---\n")

superconductors = [
    # (name, Tc(K), Delta_0(meV), xi_0(nm))
    ('Al',   1.18,  0.18, 1600),
    ('Sn',   3.72,  0.59, 230),
    ('Pb',   7.19,  1.36, 83),
    ('Nb',  9.25,  1.55, 38),
    ('NbN', 16.0,  2.4,  5),
    ('MgB2', 39.0, 7.1,  5),
    ('YBCO', 92.0, 20.0, 1.5),
]

print(f"{'SC':<8} {'Tc(K)':<8} {'Delta(meV)':<12} {'f_J (THz)':<12} {'gamma/f at T=0.5Tc':<20} {'Entity?':<8}")
print("-" * 68)

for name, Tc, Delta_meV, xi_nm in superconductors:
    Delta_J = Delta_meV * 1.602e-22  # meV to J
    f_J = 2 * Delta_J / (2 * np.pi * hbar)  # Josephson frequency (Hz)
    f_THz = f_J / 1e12

    # Ginzburg-Landau lifetime at T = 0.5*Tc
    # tau_GL ~ hbar / (k_B * (Tc - T)) for T near Tc
    # At T = 0.5*Tc: tau_GL ~ hbar / (0.5 * k_B * Tc)
    T = 0.5 * Tc
    tau_GL = hbar / (k_B * (Tc - T))  # Very rough GL estimate
    gamma_pair = 1 / tau_GL
    gamma_over_f = gamma_pair / (f_J * 2 * np.pi)  # Compare to angular frequency

    entity = "YES" if gamma_over_f < 1 else "NO"
    print(f"{name:<8} {Tc:<8.1f} {Delta_meV:<12.2f} {f_THz:<12.2f} {gamma_over_f:<20.4f} {entity:<8}")

print("""
IMPORTANT CAVEAT: The Ginzburg-Landau estimate above is for T NEAR Tc.
At T << Tc, pair breaking is exponentially suppressed:
  gamma ~ Delta * exp(-Delta/(k_B*T))

At T = 0.1*Tc: gamma/f ~ exp(-10) ~ 5e-5 (deeply oscillatory entity)
At T = 0.9*Tc: gamma/f ~ 1 (marginal entity — approaching transition)
At T = Tc:     gamma/f -> infinity (entity destroyed)

The superconducting transition IS the oscillatory entity criterion:
  gamma/f = 1 at T = Tc (within Synchronism vocabulary)
  |Delta| = 0 at T = Tc (in BCS vocabulary)

These are the SAME statement in different languages.
""")


# =============================================================================
# PART 10: IS THIS TAUTOLOGICAL? (THE HONEST CHECK)
# =============================================================================

print("=" * 70)
print("PART 10: CIRCULARITY CHECK — IS COOPER PAIR = ENTITY JUST VOCABULARY?")
print("=" * 70)

print("""
The claim: "Cooper pair condensate is an oscillatory entity with gamma/f < 1"

Is this predictive or just vocabulary mapping?

TEST 1: Does the entity criterion predict Tc?
  gamma/f = 1 at T = Tc
  gamma = k_B*(Tc-T)/hbar near Tc (Ginzburg-Landau)
  f = 2*Delta/(2*pi*hbar) = Delta/(pi*hbar)

  Setting gamma/f = 1:
    k_B*(Tc-T)/hbar = Delta(T)/(pi*hbar)
    Delta(T) = pi*k_B*(Tc-T)

  BCS near Tc: Delta(T) ~ 3.06*k_B*Tc*sqrt(1-T/Tc)

  These give DIFFERENT temperature dependences:
    Entity criterion: Delta ∝ (Tc - T) [linear]
    BCS: Delta ∝ sqrt(Tc - T) [square root]

  VERDICT: NOT equivalent. The entity criterion gives a WRONG temperature
  dependence for the gap near Tc. BCS is correct (confirmed experimentally).

  This means the entity vocabulary DOES NOT reproduce BCS physics.
  It is a weaker statement that captures the qualitative feature
  (coherence lost at Tc) but gets the quantitative details wrong.

TEST 2: Does the two-entity framework predict anything NEW about SC?

  The structural entity criterion (L < L_crit) applies to the lattice.
  The oscillatory criterion (gamma/f < 1) applies to the condensate.

  PREDICTION: Superconductivity requires BOTH entities to coexist:
    - The lattice must be a structural entity (T < T_m, L < L_crit)
    - The condensate must be an oscillatory entity (T < Tc, gamma/f < 1)

  Is this trivial? Almost. But there's one non-trivial implication:

  NEAR MELTING, lattice anharmonicity increases gamma_phonon.
  If gamma_phonon mediates Cooper pair breaking, then:
    gamma_pair ∝ gamma_phonon (through electron-phonon coupling)

  Prediction: Tc/T_m should correlate with the structural entity's
  "distance from failure" (L/L_crit at Tc).

  For BCS superconductors: Tc << T_m (typically Tc/T_m ~ 0.001-0.01)
  This means L(Tc) << L_crit — the lattice is deeply stable when SC forms.
  No interesting interplay.

  For high-Tc (cuprates): Tc/T_m ~ 0.02-0.04 — still no interplay.

  VERDICT: The two-entity framework doesn't predict anything new for SC
  because Tc and T_m are separated by orders of magnitude. The lattice
  entity and the condensate entity don't "compete" in any real material.
""")


# =============================================================================
# PART 11: PHASE 4 ASSESSMENT
# =============================================================================

print("=" * 70)
print("PART 11: PHASE 4 CUMULATIVE ASSESSMENT (5 SESSIONS)")
print("=" * 70)

print("""
PHASE 4 SCORECARD:

Session | Topic                      | Finding                    | Status
--------|----------------------------|----------------------------|--------
#1      | KSS viscosity bound        | Orders quantum->classical  | CONFIRMED
#2      | Entity<->KSS mapping       | 4pi gap, no quantitative   | FAILED
#3      | Lindemann-KSS orthogonality| Two entity types exist     | GENUINE
#3      | Lindemann more universal    | L clusters 2x tighter      | CONFIRMED
#4      | Derive L_crit from substrate| Three routes, all fail     | PRODUCTIVE FAILURE
#4      | BCC > FCC > DIA ordering   | Observed but confounded    | CONFOUNDED
#4      | Bonding > structure for L   | p=0.0002 within BCC       | EMPIRICAL
#5      | Allotrope deconfounding    | Cannot deconfound (coupled)| DEAD END (productive)
#5      | Cooper pairs = oscillatory | gamma/f wrong T-dependence | VOCABULARY (not predictive)
#5      | Two-entity SC prediction   | Tc << T_m, no interplay    | NO NEW PHYSICS

GENUINE CONTRIBUTIONS (will survive):
  1. Two-entity taxonomy (oscillatory vs structural)
  2. Bonding dominates structure for L_crit
  3. theta_D-structure coupling is thermodynamic, not accidental

VOCABULARY MAPPINGS (organizational only):
  4. KSS bound <-> entity criterion (qualitative only)
  5. Cooper pair <-> oscillatory entity
  6. Crystal <-> structural entity

PRODUCTIVE FAILURES (mapped boundaries):
  7. L_crit cannot be derived parameter-free
  8. Entity <-> KSS has unresolvable 4pi gap
  9. Allotrope test cannot deconfound bonding vs structure
  10. Two-entity framework adds nothing to superconductivity

DIMINISHING RETURNS ASSESSMENT:

  Sessions #1-3 produced the two-entity taxonomy — genuine finding.
  Sessions #4-5 attempted to extend it — hit hard limits in all directions.

  The pattern is clear: Phase 4 found ONE genuine insight (two entity types)
  and has spent two sessions confirming it can't be extended further.

  RECOMMENDATION: CLOSE Phase 4 after this session.

  The chemistry track's cumulative contributions across all phases:
  - Phase 1: Cataloguing (organizational, no novel physics)
  - Phase 2: Failure analysis (four-regime classification, eight combined predictions)
  - Phase 3: N-S circularity proof (definitively rules out CFD as predictive tool)
  - Phase 4: Two-entity taxonomy (oscillatory vs structural entities)

  These are real contributions to the Synchronism framework's self-understanding.
  They are predominantly NEGATIVE results: mapping what the framework CANNOT do.
  This is valuable — it prevents future waste on the same dead ends.
""")

print("\n" + "=" * 70)
print("SESSION #5 COMPLETE")
print("=" * 70)
